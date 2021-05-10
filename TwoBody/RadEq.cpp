#include "RadEq.h"
#include "HSpline5.h"
#include "HybridBasis.h"


RadEq::RadEq(const Pair &pair, const Int l, const double ene) : \
config(FaddeevConfigurator::getInstance()), pair(pair), l(l), ene(ene) {
	//energy conversion factor
	eConv = config.eConv;
	k = sqrt(ene * eConv);
	double xConv;
	xConv = 1.0;
	scaleCoo(xConv);
	eta = xConv * eConv * 0.5 * pair.q12 / k;
	coul = new CWFCalculator(l, eta);
	makeDiscretization();
	makePotentials();
	solve();
}

void RadEq::makeDiscretization() {
	double xmax = config.xmax[pair.alpha];

	scaleCoo(xmax);
	Grid gr(0.0, xmax, config.nx[pair.alpha]);
	//mapping
	config.xmap(pair.alpha, gr);
	array<shared_ptr<const ABasis<Complex>>, 1> bases;
	//bases[0] = make_shared<const HSpline5<Complex>>(gr, Dir, Mix, 0.0, zi*k);
	//bases[0] = make_shared<const HSpline5<Complex>>(gr, Dir, Dir, 0.0, 0.0);
	//bases[0] = make_shared<const HSpline5<Complex>>(gr, Dir, none, 0.0, 0.0);
	//bases[0] = make_shared<const HSpline5<Complex>>(gr, none, none, 0.0, 0.0);
	
	
	unique_ptr<ABasis<Complex>> bas = \
		make_unique<HSpline5<Complex>>(gr, Dir, none, zzero, zzero);
	vector<shared_ptr<CWFCalculator>> cwfs; cwfs.reserve(1);
	cwfs.push_back(make_shared<CWFCalculator>(l, eta));
	//cwfs.push_back(make_shared<CWFCalculator>(l+0.5, eta));
	vector<double> pn; pn.resize(1);
	pn[0] = k;
	vector<basset> bs; bs.resize(1);
	bs[0] = basset{0, false};
	vector<double> rnl; rnl.resize(1);
	rnl[0] = 2.5;
	scaleCoo(rnl[0]);
	/*
	double xConv;
	xConv = 1.0;
	scaleCoo(xConv);
	vector<shared_ptr<CWFCalculator>> cwfs; cwfs.reserve(2);
	cwfs.push_back(make_shared<CWFCalculator>(l+0.5, eta));
	cwfs.push_back(make_shared<CWFCalculator>(l, xConv * eConv * 0.5 * pair.q12 / 2.0));
	vector<double> pn; pn.resize(2);
	pn[0] = k; pn[1] = 2.0;
	vector<double> rnl; rnl.resize(2);
	rnl[0] = 2.5; rnl[1] = 2.5;
	scaleCoo(rnl[0]); scaleCoo(rnl[1]);
	*/
	/*
	double xConv;
	xConv = 1.0;
	scaleCoo(xConv);
	vector<shared_ptr<CWFCalculator>> cwfs; cwfs.reserve(4);
	cwfs.push_back(make_shared<CWFCalculator>(l, eta));
	cwfs.push_back(make_shared<CWFCalculator>(l, xConv * eConv * 0.5 * pair.q12 / 2.0));
	cwfs.push_back(make_shared<CWFCalculator>(l, xConv * eConv * 0.5 * pair.q12 / 3.0));
	cwfs.push_back(make_shared<CWFCalculator>(l, xConv * eConv * 0.5 * pair.q12 / 4.0));
	cwfs.push_back(make_shared<CWFCalculator>(l, xConv * eConv * 0.5 * pair.q12 / 5.0));
	cwfs.push_back(make_shared<CWFCalculator>(l, xConv * eConv * 0.5 * pair.q12 / 6.0));
	vector<double> pn; pn.resize(6);
	pn[0] = k; pn[1] = 2.0; pn[2] = 3.0; pn[3] = 4.0; pn[4] = 5.0; pn[5] = 6.0;
	vector<double> rnl; rnl.resize(6);
	rnl[0] = 2.5; rnl[1] = 2.5; rnl[2] = 2.5; rnl[3] = 2.5; rnl[4] = 2.5; rnl[5] = 2.5;
	scaleCoo(rnl[0]); scaleCoo(rnl[1]); scaleCoo(rnl[2]); scaleCoo(rnl[3]);
	scaleCoo(rnl[4]); scaleCoo(rnl[5]);
	*/
	bases[0] = make_shared<HybridBasis<AHermitSpline<Complex>>>(\
			std::move(bas), cwfs, pn, rnl, bs);
	
	/*
	unique_ptr<ABasis<Complex>> bas = \
		make_unique<HSpline5<Complex>>(gr, Dir, none, zzero, zzero);
	vector<unique_ptr<CWFCalculator>> cwfs; cwfs.reserve(1);
	cwfs.push_back(make_unique<CWFCalculator>(l, eta));
	vector<double> pn; pn.resize(1);
	pn[0] = k;
	vector<double> rnl; rnl.resize(1);
	rnl[0] = 2.5;
	scaleCoo(rnl[0]);
	vvv = make_unique<HybridBasis<AHermitSpline<Complex>>>(\
		std::move(bas), std::move(cwfs), pn, rnl);

	bas = \
		make_unique<HSpline5<Complex>>(gr, Dir, none, zzero, zzero);
	cwfs.resize(0);
	pn.resize(0);
	rnl.resize(0);
	bases[0] = make_shared<HybridBasis<AHermitSpline<Complex>>>(\
		std::move(bas), std::move(cwfs), pn, rnl);
	*/

	discr = make_shared<Collocation<Complex, 1, Complex>>(bases);
	basis = discr->bases[0].get();
	cGrid = discr->cGrids[0];
	nx = discr->getNi(0);

	//gr.print();
	//cGrid.print();
}

void RadEq::makePotentials() {
	double xConv;
	xConv = 1.0;
	scaleCoo(xConv);
	xConv = 1.0 / xConv;
	vc = make_shared<ScaledPot<Complex>>(\
		pair.vc, eConv, xConv);
	vs = make_shared<ScaledPot<Complex>>(\
		pair.vs, eConv, xConv);
}

void RadEq::scaleCoo(double &x) const {
	//physical to reduced coordinates
	pair.reduceCoo(x);
	//scale coordinates due to energy conversion factor
	x /= sqrt(eConv);
}

void RadEq::solve() {
	sol = make_shared<Function<Complex, 1>>(discr);

	GenMatrix<Complex> h(nx);
	Vector<Complex> y(nx);
	NumbersList limits;
	NumbersList::const_iterator c;
	double xi; Int j;
	Complex fi;

	h.fill(zzero); y.fill(zzero);
	Int il = l;
	double k2 = k * k;
	for (Int i = 0; i < nx; i++) {
		xi = cGrid[i];
		basis->getNonzero(xi, limits);
		for (c = limits.begin(); c != limits.end(); c++) {
			j = *c;
			fi = basis->f(xi, j);
			h[i][j] = -basis->dd(xi, j) \
				+ fi * (il*(il + 1.0) / (xi*xi) + (*vc)(xi) + (*vs)(xi) - k2);
		}
		y[i] = -(*vs)(xi)*coul->flExp(k*xi);
		//y[i] = -(*vs)(xi)*coul->flExp(k*xi) - \
			abig* \
			(- vvv->dd(xi, vvv->getNCoef()-1) \
			+ vvv->f(xi, vvv->getNCoef() - 1) * \
				(il*(il + 1.0) / (xi*xi) + (*vc)(xi) + (*vs)(xi) - k2) );
	}

	h.solve(y);
	sol->coef = y;
	solved = true;
}

void RadEq::writeSol() const {
	const string filename = "solution.dat";
	shared_ptr<Function<Complex, 1>> solf;
	solf = dynamic_pointer_cast<Function<Complex, 1>>(sol);
	double a = 0.0;
	double b = config.xmax[pair.alpha];
	const Int nPoi = WRITE_SOL_NPOI_RADEQ;
	double h = (b - a) / (nPoi - 1);
	ofstream fout;
	fout.open(filename);
	std::array<double, 1> xs;
	Complex val;
	solf->coef.print();
	//Residue
	//for (Int i = 0; i < cGrid.getNPoints(); i++) {
	//	val = -solf->d2f(array<double, 1>({cGrid[i]})) + \
			(*solf)(array<double, 1>({ cGrid[i] })) * \
			(l*(l + 1.0) / (cGrid[i] * cGrid[i]) + (*vc)(cGrid[i]) + (*vs)(cGrid[i]) - k * k) + \
			(*vs)(cGrid[i])*coul->flExp(k*cGrid[i]) + \
			abig* \
			(-vvv->dd(cGrid[i], vvv->getNCoef() - 1) \
				+ vvv->f(cGrid[i], vvv->getNCoef() - 1) * \
				(l*(l + 1.0) / (cGrid[i] * cGrid[i]) + (*vc)(cGrid[i]) + (*vs)(cGrid[i]) - k*k));
	//	cout << cGrid[i] << "  " << abs(val) << endl;
	//}
	//cout << endl;
	for (Int i = 0; i < nPoi; i++) {
		xs[0] = a + i * h;
		fout << xs[0] << "  ";
		scaleCoo(xs[0]);
		val = (*solf)(xs);
		//val += coul->flExp(k*xs[0]);

		//Residue
		//if (xs[0] > 0.3)
		//	val = -solf->d2f(xs) + \
				(*solf)(xs) * (l*(l + 1.0) / (xs[0] * xs[0]) + (*vc)(xs[0]) + (*vs)(xs[0]) - k*k) + \
					(*vs)(xs[0])*coul->flExp(k*xs[0]);
		//else
		//	val = zzero;

		//Exact solution
		//if (xs[0] >= 2.5)
		//	val = abig*coul->ulp(k*xs[0]);
		//else
		//	val = asmall*coul0.fl(k*xs[0]) - coul->flExp(k*xs[0]);

		//Exact - add. bas. func.
		//if (xs[0] >= 2.5)
		//	val = abig*coul->ulp(k*xs[0]) - abig* vvv->f(xs[0], vvv->getNCoef() - 1);
		//else
		//	val = asmall*coul0.fl(k*xs[0]) - coul->flExp(k*xs[0]) - \
										abig * vvv->f(xs[0], vvv->getNCoef() - 1);

		//RHS
		//val = -(*vs)(xs[0])*coul->flExp(k*xs[0]) - \
			abig* \
			(-vvv->dd(xs[0], vvv->getNCoef() - 1) \
				+ vvv->f(xs[0], vvv->getNCoef() - 1) * \
				(l*(l + 1.0) / (xs[0] * xs[0]) + (*vc)(xs[0]) + (*vs)(xs[0]) - k*k));

		fout << pow(abs(val), 2.0) << endl;
		//fout << val.real() << endl;
	}
	fout << endl;
	fout.close();
}

/*
void RadEq::test() const {
	Function<Complex, 1> func(discr);

	CWFCalculator coul0(l, zzero);
	Complex abig = exp(zi*sigmaL(l, eta)) * \
		(coul->fl(k*2.5)*coul0.dfl(k*2.5) - coul->dfl(k*2.5)*coul0.fl(k*2.5)) / \
		(coul0.fl(k*2.5)*coul->dulp(k*2.5) - coul0.dfl(k*2.5)*coul->ulp(k*2.5));
	Complex asmall = exp(zi*sigmaL(l, eta)) * \
		(coul->fl(k*2.5)*coul->dulp(k*2.5) - coul->dfl(k*2.5)*coul->ulp(k*2.5)) / \
		(coul0.fl(k*2.5)*coul->dulp(k*2.5) - coul0.dfl(k*2.5)*coul->ulp(k*2.5));

	const Grid &gr = discr->bases[0]->getGrid();
	Complex val;
	Int m = 3*( gr.getNPoints() - 1 ) - 1;
	assert(discr->bases[0].get()->getNCoef() == m);
	GenMatrix<Complex> vdm(m);
	for (Int i = 0; i < gr.getNPoints()-1; i++) {
		double poi = gr[i];
		
		if (i != 0) {
			for (Int j = 0; j < m; j++) {
				vdm[3 * i-1][j] = discr->bases[0].get()->f(poi, j);
			}
			if (poi >= 2.5)
				val = abig * coul->ulp(k*poi) - abig * vvv->f(poi, vvv->getNCoef() - 1);
			else
				val = asmall * coul0.fl(k*poi) - coul->flExp(k*poi) - \
				abig * vvv->f(poi, vvv->getNCoef() - 1);
			func.coef[3 * i-1] = val;
		}

		for (Int j = 0; j < m; j++) {
			vdm[3 * i][j] = discr->bases[0].get()->d(poi, j);
		}
		if (poi >= 2.5)
			val = abig * k * coul->dulp(k*poi) - abig * vvv->d(poi, vvv->getNCoef() - 1);
		else
			val = asmall * k * coul0.dfl(k*poi) - k * exp(zi*sigmaL(l, eta)) * coul->dfl(k*poi) - \
			abig * vvv->d(poi, vvv->getNCoef() - 1);
		func.coef[3 * i] = val;

		for (Int j = 0; j < m; j++) {
			vdm[3 * i + 1][j] = discr->bases[0].get()->dd(poi, j);
		}
		if (poi >= 2.5)
			val = abig * k * k *coul->ddulp(k*poi) - abig * vvv->dd(poi, vvv->getNCoef() - 1);
		else
			val = asmall * k * k * coul0.ddfl(k*poi) - k * k * exp(zi*sigmaL(l, eta)) * coul->ddfl(k*poi) - \
			abig * vvv->dd(poi, vvv->getNCoef() - 1);
		func.coef[3 * i + 1] = val;
	}
	//vdm.write("RAD_MATR.dat");
	//func.coef.print();
	//vdm.solve(func.coef);
	//func.coef.print();

	Int fnum = 1;
	ofstream f;
	f.open("RNL.dat");
	double a = discr->bases[0].get()->getGrid().getLeftmostPoint();
	double b = discr->bases[0].get()->getGrid().getRightmostPoint();
	Int n = 1001;
	double h = (b - a) / (n - 1);
	double x;
	for (Int i = 0; i < n; i++) {
		x = a + i * h;
		//f << x << "  " << vvv->f(x, vvv->getNCoef()-1).real();
		//f << "  " << vvv->d(x, vvv->getNCoef() - 1).real();
		//f << "  " << vvv->dd(x, vvv->getNCoef() - 1).real() << endl;
		f << x << "  " << discr->bases[0].get()->f(x, func).real() << endl;
	}
	f << endl << endl;
	f.close();
}
*/

bool RadEq::needsSolution() {
	return !solved;
}

RadEq::~RadEq() {
	if (coul != nullptr)
		delete coul;
}
