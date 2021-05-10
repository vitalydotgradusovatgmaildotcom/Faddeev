#include "FMEq.h"
#include "HSpline5.h"
#include "RadEigenProb.h"
#include "CWFCalculator.h"
#include "HybridBasis.h"
#include "wignerSymbols-cpp.h"
#include "SimplePoly.h"


FMEq::FMEq(const System3Body &sys, const double energy) \
	: config(FaddeevConfigurator::getInstance()), \
		ene(energy), sys(sys), Mmin((1 - sys.tau) / 2), Mmax(sys.J), \
			nM(Mmax - Mmin + 1) {
	nComp = config.nComp;
	solComps.reserve(nComp);
	for (Int alpha = 0; alpha < nComp; alpha++)
		solComps.push_back(alpha);
	//energy conversion factor
	eConv = config.eConv;
	ene *= eConv;

	prepareIdentical();
	makeChannels();
	preparePartial();
	makeDiscretization();
	makePotentials();
	makeOperators();

//===================================================
}

bool FMEq::needsSolution() {
	if (lamtit != lambdat.end()) {
		cout << "Partial incoming wave " << *lamtit << " of [" << \
			lambdat.front() << ", " << lambdat.back() << "]";
		cout << endl << endl;

		solved = false;
		makeRHS();

		//multiply rhs by preconditioner
		const PrecFMOp *precFM = \
			dynamic_cast<const PrecFMOp *>(op);
		if (precFM == nullptr)
			assert(false);
		precFM->apprx->solve(*rhs);
		//rhs->coef.print();

		sol = make_shared<Components>(cdiscr);
		*sol = *rhs;

		return true;
	}
	return false;
}

void FMEq::postprocess() {
	projectBound();
	
	//writeSolution(0, -0.8);
	//writeSolution(0, 0.1);
	//writeSolution(0, 0.9);

	lamtit++;
}

void FMEq::calcCrossSect() {
	Int nCh = 0;
	for (auto alpha : solComps)
		nCh += chans[alpha].size();

	Header("Cross sections");

	double pn0 = sqrt(ene - chan0.E / eConv);
	double xConv, yConv;
	xConv = 1.0; yConv = 1.0;
	scaleCoo(xConv, yConv, chan0.alpha);
	double etan0 = 0.5*sys.q123[chan0.alpha] * eConv * yConv / pn0;

	double etan, pn;
	double sigma;
	string csName;
	Complex a;
	Int lam, lamt;
	Int m0max;
	for (Int alpha : solComps) {
		for (Int iCh = 0; iCh < chans[alpha].size(); iCh++) {
			BinChannel & ch = chans[alpha][iCh];
			m0max = min(sys.J, chan0.l);
			sigma = 0.0;
			for (Int m0 = -m0max; m0 <= m0max; m0++) {
				for (Int ilam = 0; ilam < lambda[alpha][iCh].size(); ilam++) {
					lam = lambda[alpha][iCh][ilam];
					a = zzero;
					for (Int ilamt = 0; ilamt < lambdat.size(); ilamt++) {
						lamt = lambdat[ilamt];
						a += sqrt(2 * lamt + 1)*pow(zi, lamt)*exp(zi * sigmaL(lamt, etan0))* \
							WignerSymbols::clebschGordan(lamt, chan0.l, sys.J, 0, -m0, -m0)* \
							blamlam[alpha][iCh][ilam][ilamt];
					}
					a *= sqrt(2 * chan0.l + 1)/(2.0*PI)/pn0;
					if (ch == chan0)
						a += sqrt(4.0*PI*(2 * lam + 1))*(exp(zi * 2.0* sigmaL(lam, etan0)) - zone) / \
						2.0 / zi / pn0 * WignerSymbols::clebschGordan(lam, chan0.l, sys.J, 0, -m0, -m0);
					sigma += pow(abs(a), 2.0);
				}
			}
			//print cross section
			//sigma /= (2.0*chan0.l+1.0)*2.0*sys.m123[alpha]*PI; //IN UNITS of pi*(a_0)^2
			sigma /= (2.0*chan0.l + 1.0)*2.0*sys.m123[chan0.alpha] * PI;
			csName = std::to_string(chan0.alpha) + \
				std::to_string(chan0.n) + std::to_string(chan0.l) + " -> " + \
				std::to_string(ch.alpha) + std::to_string(ch.n) + std::to_string(ch.l);
			cout << csName << " = " << sigma;

			//if (nCh == 1 && sys.J == 0) { //print phase shift additionally
				//TODO J > 0
				//Complex delta = log(zone \
					+exp(-2.0*zi*sigmaL(chan0.l, etan0))*zi*pn0*a/(sqrt(PI*2.0*sys.m123[chan0.alpha]))) / (2.0*zi);
				//cout << delta << endl;
				//cout << ", Phase shift = " << delta.real() << endl;
			//}

			cout << endl;
		}
	}
}

void FMEq::prepareIdentical() {
	if (sys.id == three_sym || sys.id == three_asym) {
		solComps = { 0 };
	}
	if (sys.id == two_sym || sys.id == two_asym) {
		if (max(nComp - 1, (Int)1) == 2) //2 equations
			solComps = { 0, 2 };
		else
			solComps = { 0 };
	}
	//check Merkuriev cutoff
	if (sys.id != no_ident) {
		assert(config.nu[0] == config.nu[1]);
		assert(config.x0[0] == config.x0[1]);
		assert(config.y0[0] == config.y0[1]);
	}
	if (sys.id == three_sym || sys.id == three_asym) {
		assert(config.nu[0] == config.nu[2]);
		assert(config.x0[0] == config.x0[2]);
		assert(config.y0[0] == config.y0[2]);
	}
}

void FMEq::makeChannels() {
	Header("Channels");
	double enePhys = ene * eConv;
	for (auto alpha : solComps) {
		cout << "Pair " << alpha << " bound states:" << endl;
		Pair pair = sys.getPair(alpha);
		Int l = 0; bool cond = true;
		//std::vector<double> rms;
		while (cond) {
			RadEigenProb tbep(pair, l);
			Vector<Complex> eval = tbep.getEval();
			Int iCh = 0;
			//tbep.writeSolution();
			while (iCh < tbep.numberOfBound() && eval[iCh].real() < enePhys) {
				chans[alpha].push_back(BinChannel(alpha, l + iCh, l, 0, eval[iCh].real(), \
					std::dynamic_pointer_cast<Function<Complex, 1>>(tbep.getEvec(iCh))));
				//rms.push_back(tbep.rms(iCh));
				iCh++;
			}
			cond = iCh > 0;
			l++;
		}
		std::vector<BinChannel> &chalph = chans[alpha];
		std::stable_sort(chans[alpha].begin(), chans[alpha].end(), \
			[&chalph](BinChannel ch1, BinChannel ch2) {return ch1.E < ch2.E; });
		//print channels
		for (Int k = 0; k < chans[alpha].size(); k++) {
			chans[alpha][k].print();
			//cout << "rms radius: " << rms[k] << endl;
		}
		cout << endl;
		//get initial channel
		for (auto ch : chans[alpha])
			if (config.chan0 == ch)
				chan0 = ch;
		//get nondegenerate channels
		shared_ptr<const ZeroPot<Complex>> tmp = \
			dynamic_pointer_cast<const ZeroPot<Complex>>(sys.vs[alpha]);
		if (tmp == nullptr)//Coulomb + short-range
			chans_nond[alpha] = chans[alpha];
		else //only Coulomb potential, degeneracy in l
			for (auto ch : chans[alpha])
				if (ch.l == 0)
					chans_nond[alpha].push_back(ch);
	}
	//print initial channel
	assert(chan0.E != 0.0);
	header("Initial channel");
	chan0.print();
}

void FMEq::preparePartial() {
	//find partial incoming waves indices lambda tilde
	Int lambdat_ = abs(sys.J - chan0.l);
	if ((lambdat_ + chan0.l + sys.J + (1 - sys.tau) / 2) % 2 == 1)
		lambdat_++;
	while (lambdat_ <= sys.J + chan0.l) {
		lambdat.push_back(lambdat_);
		lambdat_ += 2;
	}
	lamtit = lambdat.begin();

	//for each open channel nl, find indices lambda
	//of partial decomposition of scattering amplitude
	Int lambda_, l;
	for (Int alpha : solComps) {
		lambda[alpha].resize(chans[alpha].size());
		blamlam[alpha].reserve(chans[alpha].size());
		for (Int iCh = 0; iCh < chans[alpha].size(); iCh++) {
			l = chans[alpha][iCh].l;
			lambda_ = abs(sys.J - l);
			if ((lambda_ + l + sys.J + (1 - sys.tau) / 2) % 2 == 1)
				lambda_++;
			while (lambda_ <= sys.J + l) {
				lambda[alpha][iCh].push_back(lambda_);
				lambda_ += 2;
			}
			blamlam[alpha].push_back( \
				GenMatrix<Complex>(lambda[alpha][iCh].size(), lambdat.size()));
		}
	}
}

void FMEq::makeDiscretization() {
	Header("Discretization");
	Int nFunc = solComps.size() * nM;

	array<shared_ptr<const ABasis<Complex>>, 3> bases;
	array<shared_ptr<const ABasis<Complex>>, 1> basis;
	shared_ptr<const ADiscretization<Complex, 3>> discr;
	vector<shared_ptr<const ADiscretization<Complex, 3>>> discrs;
	shared_ptr<Collocation<Complex, 3, Complex>> cdiscr;
	discrs.reserve(nFunc);
	discry.resize(nComp);
	discrz.resize(nComp);
	double xmax, ymax;
	double pn0 = sqrt(ene - chan0.E/eConv);
	for (auto alpha : solComps) {
//NB! Same grids for all Mbar components atm
		//scale physical coordinates
		xmax = config.xmax[alpha];
		ymax = config.ymax[alpha];
		scaleCoo(xmax, ymax, alpha);
		Grid grX = makeXGrid(alpha, xmax);
		Grid grY = makeYGrid(alpha, ymax);
		Grid grZ = makeZGrid(alpha);
		discry[alpha].reserve(nM);
		discrz[alpha].reserve(nM);
		for (Int Mbar = Mmin; Mbar <= Mmax; Mbar++) {
#ifndef POLY_BASIS_Z_FMEQ
			bases[0] = make_shared<const HSpline5<Complex>>(grZ, none, none, 0.0, 0.0);
#else
			//NB power of polynomial is chosen so that
			//number of basis functions equals that in the case of spline basis
			//with the same nz specified in input file
			bases[0] = make_shared<const SimplePoly<Complex>>(grZ, 3 * config.nz[alpha] - 1);
#endif
			bases[2] = make_shared<const HSpline5<Complex>>(grX, Dir, Dir, 0.0, 0.0);
			
			/*
			if (alpha == chan0.alpha)
				bases[1] = make_shared<const HSpline5<Complex>>(grY, Dir, Mix, 0.0, zi*pn0);
				//bases[1] = make_shared<const HSpline5<Complex>>(grY, Dir, Dir, 0.0, 0.0);
			else {//alpha != alpha0
				if (chans_nond[alpha].empty())
					bases[1] = make_shared<const HSpline5<Complex>>(grY, Dir, Dir, 0.0, 0.0);
				else {
					double pn = sqrt(ene - chans_nond[alpha][0].E / eConv);
					bases[1] = make_shared<const HSpline5<Complex>>(grY, Dir, Mix, 0.0, zi*pn);
				}
			}
			*/
			

			
			if (chans_nond[alpha].empty()) //no open channels
				bases[1] = make_shared<const HSpline5<Complex>>(grY, Dir, Dir, 0.0, 0.0);
			else { //open channels
				unique_ptr<ABasis<Complex>> bas = \
					make_unique<HSpline5<Complex>>(grY, Dir, none, zzero, zzero);
				vector<shared_ptr<CWFCalculator>> cwfs; cwfs.reserve(chans_nond[alpha].size());
				vector<double> pn; pn.reserve(chans_nond[alpha].size());
				vector<double> rnl; rnl.reserve(chans_nond[alpha].size());
				vector<basset> bs; bs.reserve(chans_nond[alpha].size());
				double etan;
				double xConv, yConv;
				double xR, yR;
				xConv = 1.0; yConv = 1.0;
				scaleCoo(xConv, yConv, alpha);
				for (auto ch : chans_nond[alpha]) {
					pn.push_back(sqrt(ene - ch.E / eConv));
					etan = 0.5*sys.q123[alpha] * eConv * yConv / pn.back();
					cwfs.push_back(make_shared<CWFCalculator>(ch.l, etan));
#ifdef CORRECTION_FMEQ
					//hack
					
					if (ch.alpha == 1) {
					//if (ch.n == 0) {
						if (ch.n == 0)
							bs.push_back(basset{ 1, true });
						else
							bs.push_back(basset{ 1, false });
					}
					else {
						bs.push_back(basset{ 0, false });
					}
					
					/*
					if ((ch.alpha == chan0.alpha) && (ch.n == chan0.n))
						//bs.push_back(basset{ 1, false });
						bs.push_back(basset{ 1, true });
					else
						bs.push_back(basset{1, false});
					*/
#else
					bs.push_back(basset{ 0, false });
#endif
					//e+pe- 1_below
					//xR = 1.0; yR = 15.0;
					//e+pe- 2_Ore_gap
					//if (ch.alpha == 0) {
					//	xR = 1.0; yR = 13.0;
					//}
					//else {
					//	xR = 1.0; yR = 15.0;
					//}
					//e+pe- 3_H_2_Ps_2_H_3
					//if (ch.alpha == 0) {
					//	if (ch.n == 0) {
					//		xR = 1.0; yR = 20.0;
					//	}
					//	else { //n = 1
					//		xR = 1.0; yR = 45.0;
					//	}
					//}
					//else {
					//	xR = 1.0; yR = 15.0;
					//}
					//e+Hepe- 1_below
					//xR = 1.0; yR = 10.0;
					//e+Hepe- 2_Hep_2_Ps_1
					//if (ch.n == 0) {
					//	xR = 1.0; yR = 8.0;
					//}
					//else { //n = 1
					//	xR = 1.0; yR = 20.0;
					//}
					//e+Hepe- 3_Ps_1_Hep_4
					//if (ch.alpha == 0) {
					//	if (ch.n == 0) {
					//		xR = 1.0; yR = 15.0;
					//	}
					//	else { //n = 1
					//		xR = 1.0; yR = 25.0;
					//	}
					//}
					//else {
					//	xR = 1.0; yR = 20.0;
					//}
					//scaleCoo(xR, yR, alpha);
					//rnl.push_back(yR);
					rnl.push_back(0.1*grY.getRightmostPoint()); //TODO better
				}
				bases[1] = make_shared<HybridBasis<AHermitSpline<Complex>>>(\
					std::move(bas), cwfs, pn, rnl, bs);
//==========================
				//hack
				/*
				if (alpha == 1) {
					xConv = 1.0; yConv = 1.0;
					scaleCoo(xConv, yConv, alpha);
					etan = 0.5*sys.q123[alpha] * eConv * yConv / pn.back();
					CWFCalculator cwf(chans_nond[alpha].back().l, etan);
					ofstream f;
					f.open("RNL.dat");
					double a = grY.getLeftmostPoint();
					double b = grY.getRightmostPoint();
					Int n = 1001;
					double h = (b - a) / (n - 1);
					double x;
					//Int fnum = bases[1].get()->getNCoef()-1;
					Int fnum = 207;
					//for (Int i = 400; i < n; i++) {
					for (Int i = 0; i < n; i++) {
						x = a + i * h;
						f << x << "  " << bases[1].get()->f(x, fnum).real();
						f << "  " << bases[1].get()->d(x, fnum).real();
						f << "  " << bases[1].get()->dd(x, fnum).real() << endl;
						//f << x << "  " << cwf.ulp(pn.back()*x).real();
						//f << "  " << cwf.dulp(pn.back()*x).real();
						//f << "  " << cwf.ddulp(pn.back()*x).real() << endl;
					}
					f << endl << endl;
					f.close();
				}
				*/
//==========================
			}

			cdiscr = make_shared<Collocation<Complex, 3, Complex>>(bases);
			discr = cdiscr;
			discrs.push_back(discr);

			basis[0] = bases[0];
			discrz[alpha].push_back(make_shared<Collocation<Complex, 1, Complex>>(basis));
			basis[0] = bases[1];
			discry[alpha].push_back(make_shared<Collocation<Complex, 1, Complex>>(basis));
		}
	}
	this->cdiscr = make_shared<const CompsDiscr>(discrs, nComp, Mmin, Mmax, solComps);

	//get discretization parameters
	getDiscrParams();
}

Grid FMEq::makeXGrid(const Int alpha, const double xmax) {
	Grid x(0.0, xmax, config.nx[alpha]);

	//mapping
	config.xmap(alpha, x);

	cout << "X grid, component " << alpha << ":" << endl;
	x.print();
	cout << endl;

	return x;
}

Grid FMEq::makeYGrid(const Int alpha, const double ymax) {
	
	Grid y(0.0, ymax, config.ny[alpha]);

	//mapping
	config.ymap(alpha, y);

//===========================
	/*
	if (alpha == 0) {
		vector<double> poi{ 0.00,   1.42,   3.05,   4.82,   6.70,   8.69,   10.75,  12.90,  \
			15.11,  17.39,  19.72,  22.11,  24.55,  27.03,  29.56,  32.13,  34.75,  37.40,  \
			40.09,  42.81,  45.57,  48.36,  51.19,  54.04,  56.93,  59.84,  62.79,  65.75,  \
			68.75,  71.77,  74.82,  77.89,  80.99,  84.11,  87.25,  90.42,  93.61,  96.82,  \
			100.05,  103.30,  106.57,  109.86,  113.17,  116.51,  119.86,  123.22,  126.61,  130.02, \
			133.44,  136.88,  140.34,  143.81,  147.30,  150.81,  154.33,  157.87,  161.43,  165.00, \
			168.59,  172.19,  175.81,  179.44,  183.09,  186.75,  190.43,  194.12,  197.82,  201.54};
		Grid y_(poi);
		y = y_;
	}
	else {
		vector<double> poi{ 0.00,   1.63,   3.51,   5.57,   7.76,   10.07,  12.48,  14.98,  17.57, \
			20.22,  22.95,  25.74,  28.59,  31.49,  34.45,  37.46,  40.51,  43.61,  46.76,  49.94, \
			53.17,  56.44,  59.74,  63.08,  66.46,  69.87,  73.31,  76.78,  80.29,  83.82,  87.39, \
			90.99,  94.61,  98.26,  101.94,  105.64,  109.37,  113.13,  116.91,  120.72,  124.55,
			128.40,  132.28,  136.17,  140.09,  144.04,  148.00,  151.99,  155.99,  160.02,  164.07, \
			168.14,  172.22,  176.33,  180.45,  184.60,  188.76,  192.94,  197.14,  201.36,  205.59, \
			209.85,  214.12,  218.40,  222.70,  227.02,  231.36 };
		Grid y_(poi);
		y = y_;
	}
	*/
//===========================

	cout << "Y grid, component " << alpha << ":" << endl;
	y.print();
	cout << endl;

	return y;
}

Grid FMEq::makeZGrid(const Int alpha) {
#ifdef POLY_BASIS_Z_FMEQ
	return Grid(-1.0, 1.0, 2);
#endif

	Grid res(-1.0, 1.0, config.nz[alpha]);
	//if (identical)
	//	res = Grid(0.0, 1.0, config.nz[alpha]);

	//mapping
	config.zmap(alpha, res);

	cout << "Z grid, component " << alpha << ":" << endl;
	res.print();
	cout << endl;

	return res;
}

void FMEq::makePotentials() {
	double x0, y0, xmax, ymax;
	double xConv, yConv;
	for (Int alpha = 0; alpha < 3; alpha++) {
		xConv = 1.0; yConv = 1.0;
		scaleCoo(xConv, yConv, alpha);
		xConv = 1.0 / xConv; yConv = 1.0 / yConv;
		vc[alpha] = make_shared<ScaledPot<Complex>>(\
			sys.vc[alpha], eConv, xConv);
		vs[alpha] = make_shared<ScaledPot<Complex>>(\
			sys.vs[alpha], eConv, xConv);
		//cutoff radius for effective potential
		//is taken to be 5% of ymax
		//it is done to get off problems with singularity at y = 0.0
		xmax = config.xmax[alpha]; ymax = config.ymax[alpha];
		scaleCoo(xmax, ymax, alpha);
		xConv = 1.0; yConv = 1.0;
		scaleCoo(xConv, yConv, alpha);
		vcy[alpha] = \
			new CoulombTail<Complex>(eConv*yConv*sys.q123[alpha], 0.05*ymax);
		vcy_sing[alpha] = \
			new CoulombPot<Complex>(eConv*yConv*sys.q123[alpha]);
		//Merkuriev cut-off
		x0 = config.x0[alpha]; y0 = config.y0[alpha];
		scaleCoo(x0, y0, alpha);
		mc[alpha] = \
			new MerkurievCutoff<Complex>(config.nu[alpha], x0, y0);
	}
}

void FMEq::makeOperators() {
	Header("Calculate FM l.h.s. preconditioner");
	unique_ptr<LHSApprox> lhsApp = \
		make_unique<LHSApprox>(cdiscr, vc, vs, vcy, sys, ene);

	Header("Calculate FM operator matrix");
	unique_ptr<FMOp> fmop = \
		make_unique<FMOp>(cdiscr, vc, vs, mc, sys, ene);

//===================================================
	/*
	ofstream f;
	f.open("FMOP.dat");
	Components vec(cdiscr), vvv(cdiscr);
	cout << "Matrix linear size: " << fmop->getRank() << endl;
	f << fmop->getRank() << endl;
	//for (Int i = 0; i < fmop->getRank(); i++) {
	for (Int i = 0; i < 1; i++) {
		if (i%10 == 0)
			cout << i << endl;
		//for (Int j = 0; j < fmop->getRank(); j++) {
		for (Int j = 0; j < 1; j++) {
			vec.coef.fill(zzero);
			vec[j] = zone;
			fmop->times(vec, vvv);
			vec.coef.fill(zzero);
			vec[i] = zone;
			//f << double2315<double> << (vec.coef.dotp(vvv.coef)).real() << " ";
			f << double2315<Complex> << vec.coef.dotp(vvv.coef) << " ";
		}
		f << endl;
	}
	f.close();
	*/
//===================================================

	op = new PrecFMOp(std::move(fmop), std::move(lhsApp));
	cout << "Operator rank = " << op->getRank() << endl << endl;

//===================================================
	/*
	ofstream f;
	f.open("precFMOP.dat");
	Components vec(cdiscr), vvv(cdiscr);
	cout << "Matrix linear size: " << op->getRank() << endl;
	f << op->getRank() << endl;
	for (Int i = 0; i < op->getRank(); i++) {
	//for (Int i = 0; i < 1; i++) {
		if (i%10 == 0)
			cout << i << endl;
		for (Int j = 0; j < op->getRank(); j++) {
		//for (Int j = 0; j < 1; j++) {
			vec.coef.fill(zzero);
			vec[j] = zone;
			op->times(vec, vvv);
			vec.coef.fill(zzero);
			vec[i] = zone;
			//f << double2315<double> << (vec.coef.dotp(vvv.coef)).real() << " ";
			f << double2315<Complex> << vec.coef.dotp(vvv.coef) << " ";
		}
		f << endl;
	}
	f.close();
	*/
//===================================================
}

void FMEq::makeRHS() {
	Header("Calculate FM equation r.h.s.");
	if (rhs == nullptr)
		rhs = new Components(cdiscr);
	rhs->coef.fill(zzero);

	double xi, yi;
	double zroti, yroti, xroti;
	jacobiCoo xyza, xyza0, xyzb;
	Int alpha0 = chan0.alpha;
	double pn0 = sqrt(ene - chan0.E/eConv);
	double xConv, yConv;
	xConv = 1.0; yConv = 1.0;
	scaleCoo(xConv, yConv, alpha0);
	double etan0 = 0.5*sys.q123[alpha0] * eConv * yConv / pn0;
	CWFCalculator cwf(*lamtit, etan0);
	array<double, 1> xarr;
	Complex vlxy;
	std::vector<double> ylmdiv1mz2;
	std::array<Int, 3> is;
	Int i;
	double wba;
	Complex coef, coef2, coef3, sum, fmm;
	std::vector<double> mmult, mmult2;
	GenMatrix<double> f00pm, f00mm;
	ylmdiv1mz2.resize(nM);
	mmult.reserve(nM); mmult2.reserve(nM);
	f00pm.resize(nM); f00mm.resize(nM);
	double tmp, tmp2, x2;
	double llp1 = *lamtit*(*lamtit + 1) - chan0.l*(chan0.l + 1) + sys.J*(sys.J + 1);
	for (Int im = Mmin; im <= Mmax; im++) {
		for (Int imbar = Mmin; imbar <= Mmax; imbar++) {
			f00pm[im - Mmin][imbar - Mmin] = pow(-1.0, imbar) * \
				4.0*WignerSymbols::clebschGordan(*lamtit, chan0.l, sys.J, imbar, 0, imbar) / \
									( sqrt(2.0 + 2.0*DELTA(im, 0)) * (2.0 + 2.0*DELTA(imbar, 0)));
			f00mm[im - Mmin][imbar - Mmin] = f00pm[im - Mmin][imbar - Mmin] * sys.tau*pow(-1.0, im);
		}
		mmult.push_back(2.0);
		mmult2.push_back( mmult.back() );
		tmp2 = WignerSymbols::clebschGordan(*lamtit, chan0.l, sys.J, im, 0, im) / \
																				sqrt(2.0 + 2.0*DELTA(im, 0));
		mmult.back() *= tmp2;
		tmp = (llp1 - 2.0*im*im)*tmp2;
		if (im != Mmin)
			tmp += -lambdaJM(sys.J, -im)*sqrt(1.0 + DELTA(im, 1))*lambdaJM(*lamtit, -im) * \
			WignerSymbols::clebschGordan(*lamtit, chan0.l, sys.J, im - 1, 0, im - 1) / \
				sqrt(2.0 + 2.0*DELTA(im - 1, 0));
		if (im != Mmax)
			tmp += -lambdaJM(sys.J, im)*sqrt(1.0 + DELTA(im, 0))*lambdaJM(*lamtit, im) * \
			WignerSymbols::clebschGordan(*lamtit, chan0.l, sys.J, im + 1, 0, im + 1) / \
			sqrt(2.0 + 2.0*DELTA(im + 1, 0));
		 mmult2.back() *= tmp;
	}
	for (auto alpha : solComps) {
		if (alpha != alpha0) {
			GenMatrix<Complex> vsdotxy(nx[alpha], ny[alpha]);
			std::vector<Complex> vxi, vsxi;
			vxi.reserve(nx[alpha]); vsxi.reserve(nx[alpha]);
			for (Int ix = 0; ix < nx[alpha]; ix++) {
				xi = xcGrid[alpha][ix];
				vxi.push_back( (*vc[alpha])(xi) );
				vsxi.push_back( (*vs[alpha])(xi) );
			}
			for (Int iy = 0; iy < ny[alpha]; iy++) {
				yi = ycGrid[alpha][iy];
				for (Int ix = 0; ix < nx[alpha]; ix++) {
					xi = xcGrid[alpha][ix];
					vsdotxy[ix][iy] = (vxi[ix]*(*mc[alpha])(xi, yi) + vsxi[ix])*xi*yi;
				}
			}
			std::vector<double> inv1mz2pow;
			inv1mz2pow.resize(nM);
			for (Int iz = 0; iz < nz[alpha]; iz++) {
				xyza.z = zcGrid[alpha][iz];
				is[0] = iz;
				for (Int im = Mmin; im <= Mmax; im++)
					inv1mz2pow[im-Mmin] = 1.0 / pow(1.0 - xyza.z * xyza.z, 0.5*im);
				for (Int iy = 0; iy < ny[alpha]; iy++) {
					is[1] = iy;
					yi = xyza.y = ycGrid[alpha][iy];
					for (Int ix = 0; ix < nx[alpha]; ix++) {
						is[2] = ix;
						xi = xyza.x = xcGrid[alpha][ix];
						sys.rotReduced(xyza, alpha, xyza0, alpha0);
						zroti = xyza0.z; xroti = xyza0.x;  yroti = xyza0.y;
						xarr[0] = xroti;
						coef = -(vsdotxy[ix][iy] / xroti / yroti) * (*chan0.radWF)(xarr) * cwf.fl(pn0*yroti);
						wba = acos((sys.s[alpha0][alpha] * yi*xyza.z + sys.c[alpha0][alpha] * xi) / xroti);
						if ((alpha0 - alpha + 1) % 3 != 0)
							wba = 2 * PI - wba;
						for (Int im = Mmin; im <= Mmax; im++) {
							i = cdiscr->getRaw(alpha, im, is);
							coef2 = coef * pow(-1.0, im)*inv1mz2pow[im - Mmin];
							sum = 0.0;
							for (Int imbar = Mmin; imbar <= Mmax && imbar <= *lamtit; imbar++) {
								fmm = f00pm[im - Mmin][imbar - Mmin] * wignerDSmall(sys.J, imbar, im, wba) + \
									f00mm[im - Mmin][imbar - Mmin] * wignerDSmall(sys.J, imbar, -im, wba);
								sum += pLegendre(*lamtit, imbar, zroti)*fmm;
							}
							rhs->coef[i] = coef2*sum;
						}
					}
				}
			}
		}
		else { //alpha = alpha0
			std::vector<Complex> phin0l0; phin0l0.reserve(nx[alpha]);
			for (Int ix = 0; ix < nx[alpha]; ix++) {
				xarr[0] = xcGrid[alpha][ix];
				phin0l0.push_back((*chan0.radWF)(xarr));
			}
			std::vector<Complex> vhat; vhat.reserve(ny[alpha]);
			std::vector<Complex> fl; fl.reserve(ny[alpha]);
			for (Int iy = 0; iy < ny[alpha]; iy++) {
				yi = ycGrid[alpha][iy];
				vhat.push_back((*vcy_sing[alpha])(yi));
				fl.push_back(cwf.fl(pn0*yi));
			}
//======================
			/*
			ofstream f("FL.dat");
			//Int n = 2001;
			//double a = ybasis[alpha]->getGrid().getLeftmostPoint();
			//double b = ybasis[alpha]->getGrid().getRightmostPoint();
			//double h = (b - a) / (n - 1);
			//double y;
			//for (Int i = 0; i < n; i++) {
			//	y = a + i * h;
			//	f << y << "  " << cwf0.fl(pn0*y).real() << "  " << cwf0.fl(pn0*y).imag() << endl;
			//}
			for (Int iy = 0; iy < ny[alpha]; iy++) {
				yi = ycGrid[alpha][iy];
				f << yi << "  " << fl0[iy].real() << "  " << fl0[iy].imag() << endl;
			}
			f.close();
			*/
/*
			ofstream f("FL.dat");
			for (Int ix = 0; ix < nx[alpha]; ix++) {
				f << xcGrid[alpha][ix] << "  " << phin0l0[ix].real() << "  " << \
					phin0l0[ix].imag() << endl;
			}
			f.close();
*/		
//======================
			for (Int iz = 0; iz < nz[alpha]; iz++) {
				xyza.z = zcGrid[alpha][iz];
				is[0] = iz;
				for (Int im = Mmin; im <= Mmax && im <= *lamtit; im++)
					ylmdiv1mz2[im-Mmin] = \
						pLegendre(*lamtit, im, xyza.z)/pow(1.0- xyza.z*xyza.z, 0.5*im);
				for (Int iy = 0; iy < ny[alpha]; iy++) {
					xyza.y = ycGrid[alpha][iy];
					is[1] = iy;
					for (Int ix = 0; ix < nx[alpha]; ix++) {
						is[2] = ix;
						xyza.x = xcGrid[alpha][ix];
						x2 = xyza.x * xyza.x;
						vlxy = zzero;
						for (Int gamma = 0; gamma < nComp; gamma++) {
							sys.rotReduced(xyza, alpha, xyzb, gamma);
							xroti = xyzb.x; yroti = xyzb.y;
							if (gamma != alpha)
								vlxy += (*vc[gamma])(xroti)*mc[gamma]->tail(xroti, yroti);
						}
						for (Int gamma = nComp; gamma < 3; gamma++) {
							sys.rotReduced(xyza, alpha, xyzb, gamma);
							xroti = xyzb.x;
							vlxy += (*vc[gamma])(xroti) + (*vs[gamma])(xroti);
						}
						for (Int im = Mmin; im <= Mmax && im <= *lamtit; im++) {
							i = cdiscr->getRaw(alpha, im, is);
							coef3 = (vlxy - vhat[iy])*mmult[im - Mmin] + mmult2[im-Mmin]/x2;
							rhs->coef[i] = -coef3 * phin0l0[ix] * fl[iy] * ylmdiv1mz2[im - Mmin];
						}
					}
				}
			}
		}
	}
	//rhs->coef.write("RHS.dat");
}

void FMEq::getDiscrParams() {
	nx.reserve(nComp); ny.reserve(nComp); nz.reserve(nComp);
	xbasis.reserve(nComp); ybasis.reserve(nComp); zbasis.reserve(nComp);
	xcGrid.reserve(nComp); ycGrid.reserve(nComp); zcGrid.reserve(nComp);
	for (Int alpha = 0; alpha < nComp; alpha++) {
		nz.push_back(cdiscr->get(alpha, Mmin).getNi(0));
		ny.push_back(cdiscr->get(alpha, Mmin).getNi(1));
		nx.push_back(cdiscr->get(alpha, Mmin).getNi(2));
		const Collocation<Complex, 3, Complex> &cd = \
			dynamic_cast<const Collocation<Complex, 3, Complex> &> \
			(cdiscr->get(alpha, Mmin));
		zbasis.push_back(cd.bases[0].get());
		ybasis.push_back(cd.bases[1].get());
		xbasis.push_back(cd.bases[2].get());
		zcGrid.push_back(cd.cGrids[0]);
		ycGrid.push_back(cd.cGrids[1]);
		xcGrid.push_back(cd.cGrids[2]);
	}
}

void FMEq::projectBound() {
	
	vector<Function<Complex, 1>> projs;
	GenMatrix<Complex> m, mtilde;
	Vector<Complex> mcnl;
	std::vector<Vector<Complex>> mpu;
	std::array<Int, 3> ijkarr;
	Int ijkm; Complex sum;
	Int iCh;
	for (auto alpha : solComps) {
		m = xbasis[alpha]->getOverlap();
		//make partition of unity
#ifndef POLY_BASIS_Z_FMEQ
		const AHermitSpline<Complex> *zbasis_;
		zbasis_ = dynamic_cast<const AHermitSpline<Complex> *>(zbasis[alpha]);
		if (zbasis_ == nullptr)
			assert(false);
		Int nDer = zbasis_->getNDeriv();
		vector<Complex> val; val.resize(zbasis_->getGrid().getNPoints()*(nDer + 1));
		std::fill(val.begin(), val.end(), zzero);
		for (Int i = 0; i < zbasis_->getGrid().getNPoints(); i++)
			val[i*(nDer + 1)] = zone;
		Function<Complex, 1> pu(discrz[alpha][0]);
		zbasis_->interpH(pu, val);
#else
		const SimplePoly<Complex> *zbasis_;
		zbasis_ = dynamic_cast<const SimplePoly<Complex> *>(zbasis[alpha]);
		if (zbasis_ == nullptr)
			assert(false);
		Function<Complex, 1> pu(discrz[alpha][0]);
		pu.coef.fill(zzero);
		pu.coef[0] = zone / zbasis_->f(0.0, 0);
#endif
		//pu.coef.print();
		iCh = 0;
		for (auto chan : chans[alpha]) {
			mcnl = chan.radWF->coef;
			mcnl.conj();
			mcnl = m * mcnl;

			//for (Int im = Mmin; im <= Mmax; im++) {
			for (Int im = Mmin; im <= Mmin; im++) { //im = Mmin
			//for (Int im = Mmax; im <= Mmax; im++) { //im = Mmax
				projs.clear();
				projs.shrink_to_fit();
				mpu.clear(); mpu.shrink_to_fit(); mpu.reserve(lambda[alpha][iCh].size());

				for (auto lam : lambda[alpha][iCh]) {
					projs.push_back(Function<Complex, 1>(discry[alpha][im - Mmin]));
					projs.back().coef.fill(zzero);
					mtilde = zbasis_->getOverlapW(\
						[&im, &lam](double z) { return pow(1.0-z*z, 0.5*im)*pLegendre(lam, im, z); }, lam+im);
					mpu.push_back( mtilde * pu.coef );
				}

				for (Int j = 0; j < ny[alpha]; j++) {
					ijkarr[1] = j;
					for (Int k = 0; k < nz[alpha]; k++) {
						ijkarr[0] = k;
						sum = zzero;
						for (Int i = 0; i < nx[alpha]; i++) {
							ijkarr[2] = i;
							ijkm = cdiscr->getRaw(alpha, im, ijkarr);
							sum += mcnl[i] * sol->coef[ijkm];
						}
						for (Int ilam = 0; ilam < lambda[alpha][iCh].size(); ilam++) {
							projs[ilam].coef[j] += mpu[ilam][k] * sum;
						}
					}
				}
				for (Int ilam = 0; ilam < lambda[alpha][iCh].size(); ilam++) {
					projs[ilam].coef *= 2.0*PI;
					//projs[ilam].coef.print();
				}
				//writeProjections(alpha, im, iCh, projs);
				getAmpls(alpha, im, iCh, projs);
			}
			iCh++;
		}
	}
}

void FMEq::writeProjections(const Int alpha, const Int mbar, const Int iCh, \
		const vector<Function<Complex, 1>> &projs) {
	string filename;
	ofstream fout;
	Complex val;
	std::array<double, 1> ys;
	double a, b, h, x;
	const Int nPoi = WRITE_PROJ_NPOI_FMEQ;

	double xConv, yConv;
	a = 0.0;
	b = ybasis[alpha]->getGrid().getRightmostPoint();
	h = (b - a) / (nPoi - 1);

	Int lam;
	BinChannel & chan = chans[alpha][iCh];
	for (Int ilam = 0; ilam < lambda[alpha][iCh].size(); ilam++) {
		lam = lambda[alpha][iCh][ilam];
		xConv = 1.0; yConv = 1.0;
		scaleCoo(xConv, yConv, alpha);
		//double pn = sqrt(ene - chan.E / eConv);
		//double etan = 0.5*sys.q123[chan.alpha] * eConv * yConv / pn;
		//CWFCalculator cwf(chan.l, etan);
		filename = "projection_lamtalphMlam_" + \
			std::to_string(*lamtit) + std::to_string(alpha) + std::to_string(mbar) + std::to_string(lam) + \
			"_out_" + std::to_string(chan.alpha) + std::to_string(chan.n) + std::to_string(chan.l) + \
			".dat";
		fout.open(filename);
		//for (Int i = 0; i < nPoi; i++) {
		for (Int i = 1; i < nPoi; i++) {
			ys[0] = a + i * h;
			fout << ys[0] / yConv << "  ";
			val = projs[ilam](ys);
			//val = cwf.ulp(pn*ys[0]);
			//val = cwf.fl(pn*ys[0]);
			//if (chan == chan0)
			//	val += cwf.fl(pn*ys[0]);
			fout << val.real() << "  " << val.imag() << endl;
		}
	}
	fout.close();
}

/*
void FMEq::writeSolution(const Int alpha, const double z) {
	assert(config.J == 0);
	ofstream fout;
	string filename = "solution_alpha_" + \
		std::to_string(alpha) + "_z_" + std::to_string(z) + ".dat";
	fout.open(filename);

	double ax, bx, hx;
	double ay, by, hy;
	const Int nPoi = WRITE_SOL_NPOI_FMEQ;
	ax = 0.0; ay = 0.0;
	bx = xbasis[alpha]->getGrid().getRightmostPoint();
	by = ybasis[alpha]->getGrid().getRightmostPoint();
	hx = (bx - ax) / (nPoi - 1); hy = (by - ay) / (nPoi - 1);
	array<double, 3> zyx;
	zyx[0] = z;
	shared_ptr<Components> sol_ = \
		dynamic_pointer_cast<Components>(sol);
	const Function<Complex, 3> &func = sol_->getF(alpha, 0);
	Complex val;
	double xConv, yConv;
	xConv = 1.0; yConv = 1.0;
	scaleCoo(xConv, yConv, alpha);
	for (Int iy = 0; iy < nPoi; iy++) {
		zyx[1] = ay + hy * iy;
		for (Int ix = 0; ix < nPoi; ix++) {
			zyx[2] = ax + hx * ix;
			fout << zyx[2] / xConv << "  " << zyx[1] / yConv << "  ";
			val = func(zyx);
			fout << val.real() << endl;
		}
		fout << endl;
		fout << endl;
	}

	fout.close();
}
*/

void FMEq::getAmpls(const Int alpha, const Int mbar, const Int iCh, \
	const vector<Function<Complex, 1>> &projs) {
	double pn0 = sqrt(ene - chan0.E / eConv);
	double xConv, yConv;
	xConv = 1.0; yConv = 1.0;
	scaleCoo(xConv, yConv, chan0.alpha);
	double etan0 = 0.5*sys.q123[chan0.alpha] * eConv * yConv / pn0;

	BinChannel  & chan = chans[alpha][iCh];
	double pn, etan;
	double ry;
	ry = ybasis[alpha]->getGrid().getRightmostPoint();
	
	pn = sqrt(ene - chan.E / eConv);
	xConv = 1.0; yConv = 1.0;
	scaleCoo(xConv, yConv, alpha);
	etan = 0.5*sys.q123[alpha] * eConv * yConv / pn;
	CWFCalculator cwf(chan.l, etan);

	Complex coef;
	Int lam;
	for (Int ilam = 0; ilam < lambda[alpha][iCh].size(); ilam++) {
		lam = lambda[alpha][iCh][ilam];
		coef = projs[ilam](array<double, 1>{ry}) / cwf.ulp(pn*ry);
		coef /= pow(zi, chan.l)*sqrt((2*chan.l+1)/(4.0*PI))* \
			WignerSymbols::clebschGordan(lam, chan.l, sys.J, mbar, 0, mbar) * \
			2.0*sqrt(pn0/pn)/(2.0*PI)/sqrt(2+2*DELTA(mbar, 0));
		//cout << "coef = " << coef << endl;
		blamlam[alpha][iCh][ilam][lamtit - lambdat.begin()] = coef;
	}
}

void FMEq::scaleCoo(double &x, double &y, const Int alpha) {
	//physical to reduced coordinates
	sys.reduceCoo(x, y, alpha);
	//scale coordinates due to energy conversion factor
	x /= sqrt(eConv);
	y /= sqrt(eConv);
}

double FMEq::lambdaJM(const Int j, const Int m) const {
	return sqrt(j*(j + 1.0) - m * (m + 1.0));
}

/*
void FMEq::polyCollocGridZ(\
	const Int alpha, Collocation<Complex, 3, Complex> &colloc) {
	//colloc.cGrids[0].print();
	std::vector<double> zeros;
	Int nCollocPoi = 3 * config.nz[alpha];
	pLegZeros(config.J, config.J + nCollocPoi, zeros);
	colloc.setCollocGrid(0, Grid(zeros));
	//colloc.cGrids[0].print();
}
*/

FMEq::~FMEq(void) {
	//for (Int k = 0; k < vc.size(); k++)
	//	if (vc[k] != nullptr)
	//		delete vc[k];
	for (Int k = 0; k < vcy.size(); k++)
		if (vcy[k] != nullptr)
			delete vcy[k];
	for (Int k = 0; k < vcy_sing.size(); k++)
		if (vcy_sing[k] != nullptr)
			delete vcy_sing[k];
	for (Int k = 0; k < mc.size(); k++)
		if (mc[k] != nullptr)
			delete mc[k];
}