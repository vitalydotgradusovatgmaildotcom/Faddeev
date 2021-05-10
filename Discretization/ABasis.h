#pragma once

#include "cph.h"
#include "Grid.h"
#include "AFunction.h"
//#include <forward_list>
#include "GenMatrix.h"

/*struct SummationLimits {
	Int lLim;
	Int rLim;
	SummationLimits() : lLim(0), rLim(-1) {};
	SummationLimits(const Int leftLimit, const Int rightLimit) : \
		lLim(leftLimit), rLim(rightLimit) {};
};*/

typedef std::vector<Int> NumbersList;

enum bc { Dir, Neu, Mix, none };

//Template parameter defines whether basis functions
//are real or complex valued
template <class TBas>
class ABasis {
public:
	ABasis(const Grid &grid) : nBas(0), nFree(0), \
		bcLeft(NaN), bcRight(NaN), \
		bcTypeLeft(none), bcTypeRight(none), \
		grid(grid) { };
	ABasis (const Grid &gr, \
		const enum bc bcTypeLeft, const enum bc bcTypeRight, \
		const TBas bcLeft, const TBas bcRight) \
		: nBas(0), nFree(0), grid(gr), \
		bcLeft(bcLeft), bcRight(bcRight), \
		bcTypeLeft(bcTypeLeft), bcTypeRight(bcTypeRight) { };
	ABasis(const ABasis &rhs) = delete;//default;
	ABasis(ABasis &&rhs) = delete;
	//virtual ABasis & operator=(const ABasis &rhs) = default;
	ABasis & operator=(const ABasis &rhs) = delete;
	ABasis & operator=(ABasis &&rhs) = delete;
	virtual TBas f(const double x, const Int n) const = 0;
	virtual TBas d(const double x, const Int n) const = 0;
	virtual TBas dd(const double x, const Int n) const = 0;
	virtual TBas f(const double x, const AFunction<double> &func) const = 0;
	virtual Complex f(const double x, const AFunction<Complex> &func) const = 0;
	virtual TBas d(const double x, const AFunction<double> &func) const = 0;
	virtual Complex d(const double x, const AFunction<Complex> &func) const = 0;
	virtual TBas dd(const double x, const AFunction<double> &func) const = 0;
	virtual Complex dd(const double x, const AFunction<Complex> &func) const = 0;
	virtual void getNonzero(const double x, NumbersList &limits) const = 0;
			//limits contains numbers of basis functions
			//that are nonzero on Interval that contains x
			//numbers must be sorted!
	virtual Int nFuncOnInterval() const = 0;
	const Grid & getGrid() const { return grid; };
	virtual Grid collocGrid() const = 0;
	//gives matrix of size nFree which
	//does not include functions which coefficients are fixed
	virtual GenMatrix<TBas> getOverlap() const = 0;
	//gives weighted overlap matrix of size nFree
	//function f has "degree" fdeg
	virtual GenMatrix<TBas> getOverlapW(\
		std::function<double(double)> rho, const Int fdeg) const = 0;
	//returns number of free coefficients
	Int getNCoef() const { return nFree; };
	enum bc getBCTypeLeft() const { return bcTypeLeft; };
	enum bc getBCTypeRight() const { return bcTypeRight; };
	TBas getBCLeft() const { return bcLeft; };
	TBas getBCRight() const { return bcRight; };
	virtual ~ABasis(void) = default;
protected:
	Int nBas;
	Grid grid;
	Int nFree; //number of free coefficients
	TBas bcLeft, bcRight;
	enum bc bcTypeLeft, bcTypeRight;
};