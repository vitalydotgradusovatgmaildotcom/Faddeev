#include "FMEigenProb.h"
#include "HSpline5.h"
#include "SimplePoly.h"
#include "GMRESLinearEqSolver.h"
#include "SimpleEq.h"
#include "FMOpPrec.h"
#include "TriDiagM.h"


FMEigenProb::FMEigenProb(const System3Body &sys, \
	const double energy) : ene(energy),\
				config(FaddeevConfigurator::getInstance()), sys(sys) {
	nComp = config.nComp;
	solComps.reserve(nComp);
	for (Int alpha = 0; alpha < nComp; alpha++)
		solComps.push_back(alpha);
	//energy conversion factor
	eConv = config.eConv;
	ene *= eConv;

	Header("Discretization");
	prepareIdentical();
	makeDiscretization();
	makePotentials();
	makeOperators();

	evec.reserve(1);
	evec.push_back(make_unique<Components>(cdiscr));
	evec[0]->coef.setZero();
}

void FMEigenProb::prepareIdentical() {
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

void FMEigenProb::makeDiscretization() {
	Int Mmin = (1 - config.tau) / 2;
	Int Mmax = config.J;
	Int nM = Mmax - Mmin + 1;
	Int nFunc = solComps.size() * nM;

	array<shared_ptr<const ABasis<Complex>>, 3> bases;
	shared_ptr<Collocation<Complex, 3, Complex>> cdiscr;
	shared_ptr<const ADiscretization<Complex, 3>> discr;
	vector<shared_ptr<const ADiscretization<Complex, 3>>> discrs;
	discrs.reserve(nFunc);
	double xmax, ymax;
	for (auto alpha : solComps) {
		//NB! Same grids for all Mbar components atm
		//scale physical coordinates
		xmax = config.xmax[alpha];
		ymax = config.ymax[alpha];
		scaleCoo(xmax, ymax, alpha);
		Grid grX = makeXGrid(alpha, xmax);
		Grid grY = makeYGrid(alpha, ymax);
		Grid grZ = makeZGrid(alpha);
		for (Int Mbar = Mmin; Mbar <= Mmax; Mbar++) {
#ifndef POLY_BASIS_Z_FMEIGEN
			bases[0] = make_shared<const HSpline5<Complex>>(grZ, none, none, 0.0, 0.0);
#else
			//NB power of polynomial is chosen so that
			//number of basis functions equals that in the case of spline basis
			//with the same nz specified in input file
			bases[0] = make_shared<const SimplePoly<Complex>>(grZ, 3*config.nz[alpha]-1);
#endif
			bases[1] = make_shared<const HSpline5<Complex>>(grY, Dir, Dir, 0.0, 0.0);
			bases[2] = make_shared<const HSpline5<Complex>>(grX, Dir, Dir, 0.0, 0.0);

			cdiscr = make_shared<Collocation<Complex, 3, Complex>>(bases);
			discr = cdiscr;
			discrs.push_back(discr);
		}
	}
	this->cdiscr = make_shared<const CompsDiscr>(discrs, nComp, Mmin, Mmax, solComps);
}

Grid FMEigenProb::makeXGrid(const Int alpha, const double xmax) {
	Grid x(0.0, xmax, config.nx[alpha]);

	//mapping
	config.xmap(alpha, x);

	cout << "X grid, component " << alpha << ":" << endl;
	x.print();
	cout << endl;

	return x;
}

Grid FMEigenProb::makeYGrid(const Int alpha, const double ymax) {
	Grid y(0.0, ymax, config.ny[alpha]);

	//mapping
	config.ymap(alpha, y);

	cout << "Y grid, component " << alpha << ":" << endl;
	y.print();
	cout << endl;

	return y;
}

Grid FMEigenProb::makeZGrid(const Int alpha) {
#ifdef POLY_BASIS_Z_FMEIGEN
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

void FMEigenProb::makePotentials() {
	double x0, y0, xmax, ymax;
	double xConv, yConv;
	for (Int alpha = 0; alpha < 3; alpha++) {
		xConv = 1.0; yConv = 1.0;
		scaleCoo(xConv, yConv, alpha);
		xConv = 1.0 / xConv; yConv = 1.0 / yConv;
		vc[alpha] = make_shared<ScaledPot<Complex>>( \
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
		//Merkuriev cut-off
		x0 = config.x0[alpha]; y0 = config.y0[alpha];
		scaleCoo(x0, y0, alpha);
		mc[alpha] = \
			new MerkurievCutoff<Complex>(config.nu[alpha], x0, y0);
	}
}

void FMEigenProb::makeOperators() {

	Header("Calculate Identity operator matrix");
	idop = new FMIdOp(cdiscr, sys);

#ifdef DIRECT_SOLV_FMEIGEN
	unique_ptr<LHSApprox> lhsApp(nullptr);
#else
#ifndef ILU0_PREC_FMEIGEN
	Header("Calculate Preconditioner");
	unique_ptr<LHSApprox> lhsApp = \
		make_unique<LHSApprox>(cdiscr, vc, vs, vcy, sys, ene);
#else
	unique_ptr<LHSApprox> lhsApp = nullptr;
#endif
#endif

	Header("Calculate FM operator matrix");
	unique_ptr<FMOp> fmop = \
		make_unique<FMOp>(cdiscr, vc, vs, mc, sys, ene);
//===================================================
	/*
		ofstream f;
	#ifdef STORE_MATRIX_FMOP
		f.open("FMOP_store.dat");
	#else
		f.open("FMOP_nostore.dat");
	#endif
		Components vec(cdiscr), vvv(cdiscr);
		cout << "Matrix linear size: " << fmop->getRank() << endl;
		f << fmop->getRank() << endl;
		for (Int i = 0; i < fmop->getRank(); i++) {
		//for (Int i = 0; i < 1; i++) {
			if (i%10 == 0)
				cout << i << endl;
			for (Int j = 0; j < fmop->getRank(); j++) {
			//for (Int j = 0; j < 1; j++) {
				vec.coef.fill(zzero);
				vec[j] = zone;
				//lhsApp->solve(vec);
				//vvv = vec;
				//lhsApp->times(vec, vvv);
				fmop->times(vec, vvv);
				vec.coef.fill(zzero);
				vec[i] = zone;
				//f << double2315<double> << (vec.coef.dotp(vvv.coef)).real() << " ";
				f << double2315<Complex> << vec.coef.dotp(vvv.coef) << " ";
				//if (vec.coef.dotp(vvv.coef) != zzero)
				//	f << 1;
				//else
				//	f << 0;
			}
			f << endl;
		}
		f.close();
		*/
//===================================================
//Check preconditioner
	
	/*
	Components sol(cdiscr);
	sol.coef.fill(zone);
	Components rhs(cdiscr);
	rhs = sol;
	cout << "before fmop times norm = " << rhs.coef.norm() << endl;
	fmop->times(sol, rhs);
	double rhsNorm1 = rhs.coef.norm();
	cout << "after fmop times norm = " << rhsNorm1 << endl;
	lhsApp->solve(rhs);
	double rhsNorm2 = rhs.coef.norm();
	cout << "after lhsApp solve norm = " << rhsNorm2 << endl;
	cout << "(after lhsApp norm)/(after fmop norm) = " << \
		rhsNorm2 / rhsNorm1 << endl;
	cout << "Rel. error norm = " << \
		double135<double> << \
			(sol.coef - rhs.coef).norm()/sol.coef.norm() << endl;
		*/
//===================================================
#ifndef ILU0_PREC_FMEIGEN
	precFM = new PrecFMOp(std::move(fmop), std::move(lhsApp));
#else
	Header("Calculate Preconditioner");
	unique_ptr<FMOpPrec> FMPrec = \
		make_unique<FMOpPrec>(*fmop);
	precFM = new PrecFMOp(std::move(fmop), \
		nullptr, std::move(FMPrec));
#endif

	op = new IterativeOp(cdiscr, precFM, idop);
	cout << "Operator rank = " << op->getRank() << endl;
}

void FMEigenProb::scaleCoo(double &x, double &y, const Int alpha) const {
	//physical to reduced coordinates
	sys.reduceCoo(x, y, alpha);
	//scale coordinates due to energy conversion factor
	x /= sqrt(eConv);
	y /= sqrt(eConv);
}

Vector<Complex> FMEigenProb::getEval() const {
	Vector<Complex> res(eval);
	res *= eConv;
	return res;
}

shared_ptr<AFunction<Complex>> FMEigenProb::getEvec(const Int i) {
	return evec[i];
}

FMEigenProb::~FMEigenProb() {
	if (precFM != nullptr)
		delete precFM;
	if (idop != nullptr)
		delete idop;
	for (Int k = 0; k < vcy.size(); k++)
		if (vcy[k] != nullptr)
			delete vcy[k];
	for (Int k = 0; k < mc.size(); k++)
		if (mc[k] != nullptr)
			delete mc[k];
}


IterativeOp::IterativeOp(const shared_ptr<const CompsDiscr> &cdiscr, \
	PrecFMOp *precFM, FMIdOp *idop) \
	: VectOperator<Complex, 3>(cdiscr), precFM(precFM), idop(idop) { }

void IterativeOp::times(const AFunction<Complex> &u, \
	AFunction<Complex> &res) {
	
	//S(H-E^*)^{-1}
	shared_ptr<AFunction<Complex>> tmp \
		= u.clone();
	
#ifdef DIRECT_SOLV_FMEIGEN
	precFM->fmop->solve(*tmp);
#else
#ifndef ILU0_PREC_FMEIGEN
	precFM->apprx->solve(*tmp);
#else
	precFM->prec->times(u, *tmp);
#endif
	res = *tmp; //use res as rhs
	SimpleEq<Complex> eq(precFM, &res, tmp);
	GMRESLinearEqSolver<Complex> solv;
	solv.calculateSolution(eq);
#endif

	idop->times(*tmp, res);

	/*
	//(H-E^*)^{-1}S
	idop->times(u, res);

#ifdef DIRECT_SOLV_FMEIGEN
	precFM->fmop->solve(res);
#else
	precFM->lhsApp->solve(res);
	unique_ptr<AFunction<Complex>> rhs \
		= res.clone();
	SimpleEq<Complex> eq(precFM, rhs.get(), &res);
	GMRESLinearEqSolver<Complex> solv;
	solv.calculateSolution(eq);
#endif
*/
}

void IterativeOp::solve(AFunction<Complex> &rhssol) {
	assert(false);
}

/*
void FMEigenProb::polyCollocGridZ( \
		const Int alpha, Collocation<Complex, 3, Complex> &colloc) {
	//colloc.cGrids[0].print();
	std::vector<double> zeros;
	Int nCollocPoi = 3 * config.nz[alpha];
	pLegZeros(config.J, config.J + nCollocPoi, zeros);
	colloc.setCollocGrid(0, Grid(zeros));
	//colloc.cGrids[0].print();
}
*/

//Choose Gaussian collocation points 
//as close as possible (in the least squares sense)
//to associated Legendre roots
/*
void FMEigenProb::splGridZ(const Int alpha, Grid &gr) {
	std::vector<double> zeros;
	pLegZeros(config.J, config.J + 3 * config.nz[alpha], zeros);
	Int nNodes = config.nz[alpha];
	Int nCollocPoi = 3 * config.nz[alpha];
	vector<vector<double>> ksi;
	ksi.reserve(nNodes);
	vector<double> gaussP, gaussW;
	
	//More nodes at endpoints
	if (nNodes == 2) {
		gaussP.resize(6); gaussW.resize(6);
		gauleg(-1.0, 1.0, gaussP, gaussW);
		ksi.push_back(gaussP);
	}
	else {
		gaussP.resize(5); gaussW.resize(5);
		gauleg(-1.0, 1.0, gaussP, gaussW);
		ksi.push_back(gaussP);
		gaussP.resize(3); gaussW.resize(3);
		gauleg(-1.0, 1.0, gaussP, gaussW);
		for (Int k = 1; k < nNodes - 2; k++)
			ksi.push_back(gaussP);
		gaussP.resize(4); gaussW.resize(4);
		gauleg(-1.0, 1.0, gaussP, gaussW);
		ksi.push_back(gaussP);
	}
	ksi.push_back(vector<double>());
	

	//for (Int k = 0; k < ksi.size(); k++) {
	//	for (Int i = 0; i < ksi[k].size(); i++)
	//		cout << ksi[k][i] << " ";
	//	cout << endl;
	//}
	
	Vector<double> rhs(nNodes - 2);
	if (nNodes > 2) {
		TriDiagM<double> A(nNodes - 2);
		double aij;
		Int sh = 0;

		aij = 0.0;
		for (auto ksii : ksi[0])
			aij += pow(1.0 + ksii, 2.0);
		for (auto ksii : ksi[1])
			aij += pow(1.0 - ksii, 2.0);
		A.set(0, 0, 0.5*aij);

		rhs[0] = 0.0;
		for (Int i = 0; i < ksi[0].size(); i++) {
			double ksii = ksi[0][i];
			rhs[0] += (1.0 + ksii)*(zeros[sh + i] + 0.5*(1.0 - ksii));
			//cout << sh + i << " ";
		}
		//cout << endl;
		sh += ksi[0].size();

		if (nNodes > 3) {

			aij = 0.0;
			for (auto ksii : ksi[1])
				aij += 1.0 - ksii * ksii;
			A.set(0, 1, 0.5*aij);

			for (Int i = 0; i < ksi[1].size(); i++) {
				double ksii = ksi[1][i];
				rhs[0] += (1.0 - ksii)*zeros[sh + i];
				//cout << sh + i << " ";
			}
			//cout << endl;

			for (Int k = 2; k < nNodes - 2; k++) {
				aij = 0.0;
				for (auto ksii : ksi[k - 1])
					aij += 1.0 - ksii * ksii;
				A.set(k - 1, k - 2, 0.5*aij);

				aij = 0.0;
				for (auto ksii : ksi[k - 1])
					aij += pow(1.0 + ksii, 2.0);
				for (auto ksii : ksi[k])
					aij += pow(1.0 - ksii, 2.0);
				A.set(k - 1, k - 1, 0.5*aij);

				aij = 0.0;
				for (auto ksii : ksi[k])
					aij += 1.0 - ksii * ksii;
				A.set(k - 1, k, 0.5*aij);

				rhs[k - 1] = 0.0;
				for (Int i = 0; i < ksi[k - 1].size(); i++) {
					double ksii = ksi[k - 1][i];
					rhs[k - 1] += (1.0 + ksii)*zeros[sh + i];
					//cout << sh + i << " ";
				}
				//cout << endl;
				sh += ksi[k - 1].size();
				for (Int i = 0; i < ksi[k].size(); i++) {
					double ksii = ksi[k][i];
					rhs[k - 1] += (1.0 - ksii)*zeros[sh + i];
					//cout << sh + i << " ";
				}
				//cout << endl;
			}

			aij = 0.0;
			for (auto ksii : ksi[nNodes - 3])
				aij += 1.0 - ksii * ksii;
			A.set(nNodes - 3, nNodes - 4, 0.5*aij);

			aij = 0.0;
			for (auto ksii : ksi[nNodes - 3])
				aij += pow(1.0 + ksii, 2.0);
			for (auto ksii : ksi[nNodes - 2])
				aij += pow(1.0 - ksii, 2.0);
			A.set(nNodes - 3, nNodes - 3, 0.5*aij);

			rhs[nNodes - 3] = 0.0;
			for (Int i = 0; i < ksi[nNodes - 3].size(); i++) {
				double ksii = ksi[nNodes - 3][i];
				rhs[nNodes - 3] += (1.0 + ksii)*zeros[sh + i];
				//cout << sh + i << " ";
			}
			//cout << endl;
			sh += ksi[nNodes - 3].size();

		}

		for (Int i = 0; i < ksi[nNodes - 2].size(); i++) {
			double ksii = ksi[nNodes - 2][i];
			rhs[nNodes - 3] += (1.0 - ksii)*(zeros[sh + i] - 0.5*(1.0 + ksii));
			//cout << sh + i << " ";
		}
		//cout << endl;
		
		A.solve(rhs);
		assert(rhs[0] > -1.0);
		assert(rhs[nNodes - 3] < 1.0);
		for (Int k = 1; k < nNodes - 2; k++)
			assert(rhs[k - 1] < rhs[k]);
	}

	vector<double> poi;
	poi.resize(nNodes);
	poi[0] = -1.0;
	for (Int k = 1; k < nNodes - 1; k++)
		poi[k] = rhs[k-1];
	poi[nNodes - 1] = 1.0;
	gr = Grid(poi);

	
	//cout << "Ass. Legendre zeros:" << endl;
	//for (auto zero : zeros)
	//	cout << double135<double> << zero << " ";
	//cout << endl;
	//cout << "Z spline nodes:" << endl;
	//for (Int k = 0; k < nNodes; k++)
	//	cout << poi[k] << endl;
	//cout << "Collocation points:" << endl;
	//for (Int k = 0; k < nNodes - 1; k++)
	//	for (auto ksii : ksi[k]) {
	//		cout << 0.5*(poi[k] + poi[k + 1]) + \
	//			0.5*(poi[k + 1] - poi[k])*ksii << " ";
	//	}
	//cout << endl;

}
*/
