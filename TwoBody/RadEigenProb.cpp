#include "RadEigenProb.h"
#include "HSpline5.h"
#include "Function.h"
#include "Algorithms.h"

RadEigenProb::RadEigenProb(const Pair &pair, const Int l) : \
	config(FaddeevConfigurator::getInstance()), pair(pair), l(l) {
	//energy conversion factor
	eConv = config.eConv;
	makeDiscretization();
	makePotentials();
	solve();
	normalizeBound();
}

Vector<Complex> RadEigenProb::getEval() const {
	Vector<Complex> res(eval);
	res *= 1.0/eConv;
	return res;
}

shared_ptr<AFunction<Complex>> RadEigenProb::getEvec(const Int i) {
	return evec[i];
}

Int RadEigenProb::numberOfBound() const {
	return nBound;
}

void RadEigenProb::writeSolution() const {
	double a = 0.0;
	double b = config.xmax[pair.alpha];
	const Int nPoi = WRITE_SOL_NPOI_RADEP;
	double h = (b - a) / (nPoi - 1);

	string filename;
	ofstream fout;
	shared_ptr<Function<Complex, 1>> evecf;
	std::array<double, 1> xs;
	Complex val;
	for (Int k = 0; k < nBound; k++) {
		filename = "solution_alpha_" + std::to_string(pair.alpha) + \
			"_l_" + std::to_string(l) + "_n_" + std::to_string(k) + ".dat";
		fout.open(filename);
		evecf = dynamic_pointer_cast<Function<Complex, 1>>(evec[k]);
		for (Int i = 0; i < nPoi; i++) {
			xs[0] = a + i * h;
			fout << xs[0] << "  ";
			scaleCoo(xs[0]);
			val = (*evecf)(xs);
			fout << val.real() << endl;
		}
		fout.close();
	}
}

double RadEigenProb::rms(const Int i) const {
	GenMatrix<double> mtilde_ = \
		basis->getOverlapW([](double x) { return x*x; }, 2);
	GenMatrix<Complex> mtilde(mtilde_.nrows());
	for (Int i_ = 0; i_ < mtilde.nrows(); i_++)
		for (Int j = 0; j < mtilde.ncols(); j++)
			mtilde[i_][j] = mtilde_[i_][j];
	mtilde_.resize(0);
	double xConv = 1.0;
	scaleCoo(xConv);
	return sqrt( (evec[i]->coef.scal(mtilde*evec[i]->coef)).real() ) / xConv;
}

void RadEigenProb::makeDiscretization() {
	double xmax = config.xmax[pair.alpha];
	scaleCoo(xmax);
	Grid gr(0.0, xmax, config.nx[pair.alpha]);
	//mapping
	config.xmap(pair.alpha, gr);
	//gr.print();
	array<shared_ptr<const ABasis<double>>, 1> bases;
	bases[0] = make_shared<const HSpline5<double>>(gr, Dir, Dir, 0.0, 0.0);
	discr = make_shared<Collocation<Complex, 1, double>>(bases);
	basis = discr->bases[0].get();
	cGrid = discr->cGrids[0];
	nx = discr->getNi(0);
}

void RadEigenProb::makePotentials() {
	double xConv;
	xConv = 1.0;
	scaleCoo(xConv);
	xConv = 1.0 / xConv;
	vc = make_shared<ScaledPot<Complex>>( \
		pair.vc, eConv, xConv);
	vs = make_shared<ScaledPot<Complex>>( \
		pair.vs, eConv, xConv);
}

void RadEigenProb::scaleCoo(double &x) const {
	//physical to reduced coordinates
	pair.reduceCoo(x);
	//scale coordinates due to energy conversion factor
	x /= sqrt(eConv);
}

void RadEigenProb::solve() {
	eval.resize(nx);
	evec.reserve(nx);
	for (Int i = 0; i < nx; i++)
		evec.push_back(make_shared<Function<Complex, 1>>(discr));

	GenMatrix<Complex> s(nx), d(nx);
	NumbersList limits;
	NumbersList::const_iterator c;
	double xi, fi; Int j;

	s.fill(zzero); d.fill(zzero);
	Int il = l;
	for (Int i = 0; i < nx; i++) {
		xi = cGrid[i];
		basis->getNonzero(xi, limits);
		for (c = limits.begin(); c != limits.end(); c++) {
			j = *c;
			fi = basis->f(xi, j);
			s[i][j] = fi;
			d[i][j] = -basis->dd(xi, j) \
				+ fi * (il*(il + 1.0) / (xi*xi) + (*vc)(xi) + (*vs)(xi));
		}
	}

	GenMatrix<Complex> wL, wR;
	genEEV(d, s, wL, eval, wR, false, true);
	for (Int j = 0; j < wR.ncols(); j++)
		for (Int i = 0; i < wR.nrows(); i++)
			//evec[i]->coef[j] = wR[i][j];
			evec[j]->coef[i] = wR[i][j];
}


void RadEigenProb::normalizeBound() {
	//reorder eval (and evec) so that real negative eigenvalues
	//and corresponding eigenvectors are first
	//normalize these eigenvectors

	Complex *p = eval.v;
	sortSimult(eval.v, eval.v + nx, evec.begin(), evec.end(), [&p](Int i, Int j) {return p[i].real() < p[j].real(); });

	GenMatrix<double> matrixOverlap = basis->getOverlap();
	Complex sum;
	nBound = 0;
	while (eval[nBound].real() < 0.0 && \
			//IS_EPS(eval[nBound].imag() / eval[nBound].real())) { //NB! severe condition
			abs(eval[nBound].imag() / eval[nBound].real()) < IM_EV_THRESHOLD ) {
		sum = zzero;
		for (Int i = 0; i < evec[nBound]->size(); i++)
			for (Int j = 0; j < evec[nBound]->size(); j++)
				sum += evec[nBound]->coef[i] * conj(evec[nBound]->coef[j]) * matrixOverlap[i][j];
		evec[nBound]->coef *= 1.0 / sqrt(sum.real());

		nBound++;
	}
		
}