#pragma once

#include "ABasis.h"
#include "SpecialFunctions.h"

//Basis of (non-normalized) Legendre polynomials!

//TBas = double in fact
template <class TBas>
class SimplePoly :
	public ABasis<TBas> {
public:
	SimplePoly(const Grid &gr, const Int deg);
	void getNonzero(const double x, NumbersList &limits) const override;
	inline Int nFuncOnInterval() const override;
	Grid collocGrid() const override;
	GenMatrix<TBas> getOverlap() const override;
	GenMatrix<TBas> getOverlapW(\
		std::function<double(double)> rho, const Int fdeg) const override;
	TBas f(const double x, const Int n) const override;
	TBas d(const double x, const Int n) const override;
	TBas dd(const double x, const Int n) const override;
	TBas f(const double x, const AFunction<double> &func) const override;
	Complex f(const double x, const AFunction<Complex> &func) const override;
	TBas d(const double x, const AFunction<double> &func) const override;
	Complex d(const double x, const AFunction<Complex> &func) const override;
	TBas dd(const double x, const AFunction<double> &func) const override;
	Complex dd(const double x, const AFunction<Complex> &func) const override;
	//val contains values of function and its derivatives at nodes in order
	//value--first derivative--second derivative--... at first node and so on
	void interpH(AFunction<TBas> &func, \
		const vector<double> &poi, const vector<TBas> &val) const;
	~SimplePoly() override = default;
protected:
	using ABasis<TBas>::bcTypeLeft;
	using ABasis<TBas>::bcTypeRight;
	using ABasis<TBas>::nBas;
	using ABasis<TBas>::nFree;
	using ABasis<TBas>::grid;
	Int deg; //degree
	double a, b;
	double twodivbma, twodivbmapow2;
	Int basPerNode;
};

template <class TBas>
inline Int SimplePoly<TBas>::nFuncOnInterval() const {
	return deg + 1;
}