#pragma once

#include "ABasis.h"
#include "SpecialFunctions.h"

//A local basis of Hermite splines
//splines must be enumerated so that
//0 - spline with value 1 at first node
//1 - spline with derivative 1 and so on
//Mixed boundary condition is of the form y'/y = bcLeft (bcRight)
template <class TBas>
class AHermitSpline :
	public ABasis<TBas> {
public:
	AHermitSpline(const Grid &gr);
	AHermitSpline(const Grid &gr, \
		const enum bc bcTypeLeft, const enum bc bcTypeRight, \
		const TBas bcLeft, const TBas bcRight);
	//AHermitSpline(const AHermitSpline &rhs) = default;
	//AHermitSpline & operator=(const AHermitSpline &rhs) = default;
	void getNonzero(const double x, NumbersList &limits) const override;
	inline Int nFuncOnInterval() const override;
	Grid collocGrid() const override;
	GenMatrix<TBas> getOverlap() const override;
	GenMatrix<TBas> getOverlapW(\
		std::function<double(double)> rho, const Int fdeg) const override;
	inline Int getNDeriv() const; //number of derivatives that a spline interpolates
	inline Int getDeg() const;
	//val contains values of function and its derivatives at nodes in order
	//value--first derivative--second derivative--... at first node and so on
	void interpH(AFunction<TBas> &func, \
		const vector<TBas> &val) const;
	~AHermitSpline(void) override = default;
protected:
	using ABasis<TBas>::bcTypeLeft, ABasis<TBas>::bcTypeRight;
	using ABasis<TBas>::bcLeft, ABasis<TBas>::bcRight;
	using ABasis<TBas>::nBas, ABasis<TBas>::nFree;
	using ABasis<TBas>::grid;
	vector<double> lPoi, mPoi, rPoi;
	Int deg; //degree
	Int basPerNode;
	void initRefPoints();
	inline Int toInternalN(const Int n) const;
};

template <class TBas>
inline Int AHermitSpline<TBas>::nFuncOnInterval() const {
	return 2*basPerNode;
}

template <class TBas>
inline Int AHermitSpline<TBas>::getNDeriv() const {
	return deg / 2;
}

template <class TBas>
inline Int AHermitSpline<TBas>::getDeg() const {
	return deg;
}

template <class TBas>
inline Int AHermitSpline<TBas>::toInternalN(const Int n) const {
	Int res;
	if (n == 0) {
		res = (bcTypeLeft == Dir) ? 1 : 0;
		res = (bcTypeLeft == Mix) ? -1 : res;
	} else {
		res = n;
		if (bcTypeLeft != none) res++;
		//if (n >= nFree-basPerNode+1) {
		if (res == nBas - basPerNode) {
			if (bcTypeRight == Dir)
				res++;
			//else if (bcTypeRight == Neu && n >= nFree-basPerNode+2)
			else if (bcTypeRight == Mix)
				res = -2;
		}
		else if (res == nBas - basPerNode + 1 && bcTypeRight != none)
			res++;
	}
	return res;
}

