//var 4 --- added basis functions on last interval
#pragma once

#include "AProjDiscr.h"
#include "Function.h"
#include "CWFCalculator.h"

//A hybrid basis for scattering calculations
//Adds functions of the form: R_{nl}(y), when y < y2_{nl},
//  F_{nl}(y, pn), when y > y2_{nl}
//	with F_{nl}(y, pn) = u_l^+(p_n y)
//  and F_{nl}(y, pn) = u_l^{+-}(p_n y) / (p_n y)^k, degree k > 0, to a given basis "bas".
//The Coulomb w.f. u_l^+(p_n y) ~ e^(-i sigma_n) e^( i theta_l(y, p_n) ),
//	with theta_l(y, p_n) = p_n y - eta_n log( 2 p_n y ) - pi l/2 + sigma_n.
//R_{nl}(y) is a linear combination of "bas" functions
//  its value and derivatives equal 0 at the leftmost point
//  and it connects smoothly with F_{nl}(y, pn) at y = y2_{nl};
//	y1_{nl} and y2_{nl} are subsequent nodes of the basis "bas"
//	such that a given y0_{nl} is in [y1_{nl}, y2_{nl})

struct basset {
	Int deg = 0;
	bool pm = false;
};

struct Fnl {
	Fnl(shared_ptr<CWFCalculator> &cwf, const double pn, \
					const Int k, const bool plus) : \
							cwf(cwf), pn(pn), k(k), plus(plus) {};
	shared_ptr<CWFCalculator> cwf;
	double pn;
	Int k; //degree k
	bool plus; //u_l^+ or u_l^-
};

template <class BasT>
class HybridBasis :
	public ABasis<Complex> {
public:
	HybridBasis(unique_ptr<const ABasis<Complex>> &&bas, \
		std::vector<shared_ptr<CWFCalculator>> &cwf, \
		const std::vector<double> &pn, \
		const std::vector<double> &y0nl, \
		const std::vector<basset> &bs);
	Complex f(const double x, const Int n) const override;
	Complex d(const double x, const Int n) const override;
	Complex dd(const double x, const Int n) const override;
	Complex f(const double x, const AFunction<double> &func) const override;
	Complex f(const double x, const AFunction<Complex> &func) const override;
	Complex d(const double x, const AFunction<double> &func) const override;
	Complex d(const double x, const AFunction<Complex> &func) const override;
	Complex dd(const double x, const AFunction<double> &func) const override;
	Complex dd(const double x, const AFunction<Complex> &func) const override;
	void getNonzero(const double x, NumbersList &limits) const override;
	Int nFuncOnInterval() const override;
	Grid collocGrid() const override;
	GenMatrix<Complex> getOverlap() const override;
	GenMatrix<Complex> getOverlapW(\
		std::function<double(double)> rho, const Int fdeg) const override;
	~HybridBasis() override;
protected:
	Int nFree0; //number of used basis functions of "bas"
	const ABasis<Complex> *bas0;
	shared_ptr<AProjDiscr<Complex, 1, Complex>> discr0;
	vector<Int> numBas0; //numbers of "bas" basis functions used (sorted)
	vector<Int> numBas0_inv;
	Int nBas1; //number of additional basis functions
	//Int maxk; //maximum degree of the Coulomb w.f. expansion correction term
	std::vector<Fnl> fnl;
	std::vector<Function<Complex, 1>> rnl;
	Function<Complex, 1> *f0; //for evaluating "bas" part of function
	std::vector<Int> inl; //numbers of intervals I_{nl}
	void makeAddFunc(std::vector<shared_ptr<CWFCalculator>> &cwf, \
		const std::vector<double> &pn, const std::vector<double> &y0nl, \
			const std::vector<basset> &bs);
	void makeRnl();
	Complex vfnl(const Int k, const double y) const;
	Complex dfnl(const Int k, const double y) const;
	Complex ddfnl(const Int k, const double y) const;
};

/* var 3 --- Rnl support on [0, y2_{nl}]
#pragma once

#include "AProjDiscr.h"
#include "Function.h"
#include "CWFCalculator.h"

//A hybrid basis for scattering calculations
//Adds functions of the form: R_{nl}(y), when y < y2_{nl},
//  F_{nl}(y, pn), when y > y2_{nl}
//	with F_{nl}(y, pn) = u_l^+(p_n y), to a given basis "bas".
//The Coulomb w.f. u_l^+(p_n y) ~ e^(-i sigma_n) e^( i theta_l(y, p_n) ),
//	with theta_l(y, p_n) = p_n y - eta_n log( 2 p_n y ) - pi l/2 + sigma_n.
//R_{nl}(y) is a linear combination of "bas" functions
//  its value and derivatives equal 0 at the leftmost point
//  and it connects smoothly with F_{nl}(y, pn) at y = y2_{nl};
//	y1_{nl} and y2_{nl} are subsequent nodes of the basis "bas"
//	such that a given y0_{nl} is in [y1_{nl}, y2_{nl})

struct Fnl {
	Fnl(unique_ptr<CWFCalculator> &&cwf, const double pn) : \
							cwf(std::move(cwf)), pn(pn) {};
	unique_ptr<CWFCalculator> cwf;
	double pn;
};

template <class BasT>
class HybridBasis :
	public ABasis<Complex> {
public:
	HybridBasis(unique_ptr<const ABasis<Complex>> &&bas, \
		std::vector<unique_ptr<CWFCalculator>> &&cwf, \
		const std::vector<double> &pn, \
		const std::vector<double> &y0nl);
	Complex f(const double x, const Int n) const override;
	Complex d(const double x, const Int n) const override;
	Complex dd(const double x, const Int n) const override;
	Complex f(const double x, const AFunction<double> &func) const override;
	Complex f(const double x, const AFunction<Complex> &func) const override;
	Complex d(const double x, const AFunction<double> &func) const override;
	Complex d(const double x, const AFunction<Complex> &func) const override;
	Complex dd(const double x, const AFunction<double> &func) const override;
	Complex dd(const double x, const AFunction<Complex> &func) const override;
	void getNonzero(const double x, NumbersList &limits) const override;
	Int nFuncOnInterval() const override;
	Grid collocGrid() const override;
	GenMatrix<Complex> getOverlap() const override;
	GenMatrix<Complex> getOverlapW(\
		std::function<double(double)> rho, const Int fdeg) const override;
	~HybridBasis() override;
protected:
	Int nFree0; //number of used basis functions of "bas"
	const ABasis<Complex> *bas0;
	shared_ptr<AProjDiscr<Complex, 1, Complex>> discr0;
	Int nBas1; //number of additional basis functions
	std::vector<Fnl> fnl;
	std::vector<Function<Complex, 1>> rnl;
	Function<Complex, 1> *f0; //for evaluating "bas" part of function
	std::vector<Int> inl; //numbers of intervals I_{nl}
	void makeAddFunc(std::vector<unique_ptr<CWFCalculator>> &&cwf, \
		const std::vector<double> &pn, const std::vector<double> &y0nl);
	void makeRnl();
};
*/

//var 2 --- Rnl support on 1 interval
/*
#pragma once

#include "AProjDiscr.h"
#include "Function.h"
#include "CWFCalculator.h"

//A hybrid basis for scattering calculations
//Adds functions of the form R_{nl}(y) F_{nl}(y, pn),
//	with F_{nl}(y, pn) = u_l^+(p_n y), to a given basis "bas".
//The Coulomb w.f. u_l^+(p_n y) ~ e^(-i sigma_n) e^( i theta_l(y, p_n) ),
//	with theta_l(y, p_n) = p_n y - eta_n log( 2 p_n y ) - pi l/2 + sigma_n.
//R_{nl}(y) equals 0 for y <= y1_{nl},
//  equals 1 for y2_{nl} <= y <= R_y
//	and changes from 0 to 1 on interval I_{nl} = [y1_{nl}, y2_{nl}];
//	y1_{nl} and y2_{nl} are subsequent nodes of the basis "bas"
//	such that a given y0_{nl} is in [y1_{nl}, y2_{nl})

struct Fnl {
	Fnl(unique_ptr<CWFCalculator> &&cwf, const double pn) : \
		cwf(std::move(cwf)), pn(pn) {};
	unique_ptr<CWFCalculator> cwf;
	double pn;
};

template <class BasT>
class HybridBasis :
	public ABasis<Complex> {
public:
	HybridBasis(unique_ptr<const ABasis<Complex>> &&bas, \
		std::vector<unique_ptr<CWFCalculator>> &&cwf, \
		const std::vector<double> &pn, \
		const std::vector<double> &y0nl);
	Complex f(const double x, const Int n) const override;
	Complex d(const double x, const Int n) const override;
	Complex dd(const double x, const Int n) const override;
	Complex f(const double x, const AFunction<double> &func) const override;
	Complex f(const double x, const AFunction<Complex> &func) const override;
	Complex d(const double x, const AFunction<double> &func) const override;
	Complex d(const double x, const AFunction<Complex> &func) const override;
	Complex dd(const double x, const AFunction<double> &func) const override;
	Complex dd(const double x, const AFunction<Complex> &func) const override;
	void getNonzero(const double x, NumbersList &limits) const override;
	Int nFuncOnInterval() const override;
	Grid collocGrid() const override;
	GenMatrix<Complex> getOverlap() const override;
	GenMatrix<Complex> getOverlapW(\
		std::function<double(double)> rho, const Int fdeg) const override;
	~HybridBasis() override;
protected:
	Int nFree0; //number of used basis functions of "bas"
	const ABasis<Complex> *bas0;
	shared_ptr<AProjDiscr<Complex, 1, Complex>> discr0;
	Int nBas1; //number of additional basis functions
	std::vector<Fnl> fnl;
	std::vector<Function<Complex, 1>> rnl;
	Function<Complex, 1> *f0; //for evaluating "bas" part of function
	std::vector<Int> inl; //numbers of intervals I_{nl}
	void makeAddFunc(std::vector<unique_ptr<CWFCalculator>> &&cwf, \
		const std::vector<double> &pn, const std::vector<double> &y0nl);
	void makeRnl();
};

*/
