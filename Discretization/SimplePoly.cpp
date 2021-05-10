#include "SimplePoly.h"

template class SimplePoly<double>;
template class SimplePoly<Complex>;

template <class TBas>
SimplePoly<TBas>::SimplePoly(const Grid &gr, const Int deg) \
	: basPerNode(0), ABasis<TBas>(gr), deg(deg) {
	assert(deg <= 100);
	assert(gr.getNPoints() == 2);
	a = grid[0]; b = grid[1];
	assert(a < b);
	twodivbma = 2.0 / (b-a);
	twodivbmapow2 = twodivbma * twodivbma;
	nFree = nBas = basPerNode = deg + 1;
}

template <class TBas>
void SimplePoly<TBas>::getNonzero(const double x, NumbersList &limits) const {
	if (x < a || x > b) {
		limits.resize(0);
		return;
	}
	limits.resize(deg+1);
	std::iota(limits.begin(), limits.end(), 0);
}

template <class TBas>
Grid SimplePoly<TBas>::collocGrid() const {
	vector<double> gaussP, gaussW;
	gaussP.resize(deg+1); gaussW.resize(deg+1);
	gauleg(a, b, gaussP, gaussW);
	return Grid(gaussP);
}

template <class TBas>
GenMatrix<TBas> SimplePoly<TBas>::getOverlap() const {
	GenMatrix<TBas> ovlpM(nFree);
	ovlpM.fill(TBas());
	for (Int k = 0; k <= deg; k++)
			ovlpM[k][k] = (b-a)/(2*k+1);
	return ovlpM;
}

template <class TBas>
GenMatrix<TBas> SimplePoly<TBas>::getOverlapW(\
	std::function<double(double)> rho, const Int fdeg) const {
	GenMatrix<TBas> ovlpM(nFree);
	ovlpM.fill(TBas());
	//
	vector<double> gaussX(deg + 1 + (fdeg + 1) / 2), \
		gaussW(deg + 1 + (fdeg + 1) / 2);
	gauleg(a, b, gaussX, gaussW);
	//
	for (Int k = 0; k <= deg; k++)
		for (Int l = k; l <= deg; l++) {
			for (Int i = 0; i < gaussX.size(); i++)
				ovlpM[k][l] += \
				gaussW[i] * f(gaussX[i], k)*f(gaussX[i], l)*rho(gaussX[i]);
			ovlpM[l][k] = ovlpM[k][l];
		}
	return ovlpM;
}

template <class TBas>
TBas SimplePoly<TBas>::f(const double x, const Int n) const {
	return pLegendre(n, twodivbma*(x-a)-1.0);
}

template <class TBas>
TBas SimplePoly<TBas>::d(const double x, const Int n) const {
	return twodivbma*dpLegendre(n, twodivbma*(x - a) - 1.0);
}

template <class TBas>
TBas SimplePoly<TBas>::dd(const double x, const Int n) const {
	return twodivbmapow2*ddpLegendre(n, twodivbma*(x - a) - 1.0);
}

template <class TBas>
TBas SimplePoly<TBas>::f(const double x, const AFunction<double> &func) const {
	std::vector<double> pvals;
	pLegendre(deg, twodivbma * (x - a) - 1.0, pvals);

	TBas res = TBas();
	for (Int k = deg; k >= 0; k--)
		res += func[k]*pvals[k];

	return res;
}

template <class TBas>
Complex SimplePoly<TBas>::f(const double x, const AFunction<Complex> &func) const {
	std::vector<double> pvals;
	pLegendre(deg, twodivbma * (x - a) - 1.0, pvals);

	Complex res = zzero;
	for (Int k = deg; k >= 0; k--)
		res += func[k] * pvals[k];

	return res;
}

template <class TBas>
TBas SimplePoly<TBas>::d(const double x, const AFunction<double> &func) const {
	std::vector<double> pvals;
	dpLegendre(deg, twodivbma * (x - a) - 1.0, pvals);

	TBas res = TBas();
	for (Int k = deg; k >= 0; k--)
		res += func[k] * pvals[k];

	return twodivbma*res;
}

template <class TBas>
Complex SimplePoly<TBas>::d(const double x, const AFunction<Complex> &func) const {
	std::vector<double> pvals;
	dpLegendre(deg, twodivbma * (x - a) - 1.0, pvals);

	Complex res = zzero;
	for (Int k = deg; k >= 0; k--)
		res += func[k] * pvals[k];

	return twodivbma*res;
}

template <class TBas>
TBas SimplePoly<TBas>::dd(const double x, const AFunction<double> &func) const {
	std::vector<double> pvals;
	ddpLegendre(deg, twodivbma * (x - a) - 1.0, pvals);

	TBas res = TBas();
	for (Int k = deg; k >= 0; k--)
		res += func[k] * pvals[k];

	return twodivbmapow2*res;
}

template <class TBas>
Complex SimplePoly<TBas>::dd(const double x, const AFunction<Complex> &func) const {
	std::vector<double> pvals;
	ddpLegendre(deg, twodivbma * (x - a) - 1.0, pvals);

	Complex res = zzero;
	for (Int k = deg; k >= 0; k--)
		res += func[k] * pvals[k];

	return twodivbmapow2*res;
}

template <class TBas>
void SimplePoly<TBas>::interpH(AFunction<TBas> &func, \
	const vector<double> &poi, const vector<TBas> &val) const {
	//TODO remove - very bad realization!
	Int n = val.size();
	Int m = poi.size();
	assert(n == deg+1);
	Int nDer = n / m - 1;
	assert(nDer <= 2); //only second derivative supported atm
	GenMatrix<TBas> vdm(n);
	Int i = 0;
	for (auto p : poi) {
		for (Int iDer = 0; iDer <= nDer; iDer++) {
			switch (iDer) {
			case 0:
				for (Int j = 0; j <= deg; j++)
					vdm[i][j] = f(p, j);
				break;
			case 1:
				for (Int j = 0; j <= deg; j++)
					vdm[i][j] = d(p, j);
				break;
			case 2:
				for (Int j = 0; j <= deg; j++)
					vdm[i][j] = dd(p, j);
				break;
			}
			func.coef[i] = val[i];
			i++;
		}
	}
	vdm.solve(func.coef);
}