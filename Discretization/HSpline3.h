#pragma once

#include "AHermitSpline.h"

template <class TBas>
class HSpline3 :
	public AHermitSpline<TBas> {
public:
	HSpline3(const Grid &gr);
	HSpline3(const Grid &gr, \
		const enum bc bcTypeLeft, const enum bc bcTypeRight, \
		const TBas bcLeft, const TBas bcRight);
	//HSpline3(const HSpline3 &rhs) = default;
	//HSpline3 & operator=(const HSpline3 &rhs) = default;
	//HSpline3 & operator=(const ABasis &rhs) override;
	TBas f(const double x, const Int n) const override;
	TBas d(const double x, const Int n) const override;
	TBas dd(const double x, const Int n) const override;
	TBas f(const double x, const AFunction<double> &func) const override;
	Complex f(const double x, const AFunction<Complex> &func) const override;
	TBas d(const double x, const AFunction<double> &func) const override;
	Complex d(const double x, const AFunction<Complex> &func) const override;
	TBas dd(const double x, const AFunction<double> &func) const override;
	Complex dd(const double x, const AFunction<Complex> &func) const override;
	~HSpline3(void) = default;
protected:
	using ABasis<TBas>::bcLeft, ABasis<TBas>::bcRight;
	using ABasis<TBas>::bcTypeLeft, ABasis<TBas>::bcTypeRight;
	using ABasis<TBas>::nBas, ABasis<TBas>::nFree;
	using AHermitSpline<TBas>::deg, AHermitSpline<TBas>::basPerNode;
	using AHermitSpline<TBas>::lPoi, AHermitSpline<TBas>::rPoi, \
																AHermitSpline<TBas>::mPoi;
	inline double phi(const double x, \
		const double xl, const double xm, const double xr) const;
	inline double dphi(const double x, \
		const double xl, const double xm, const double xr) const;
	inline double ddphi(const double x, \
		const double xl, const double xm, const double xr) const;
	inline double psi(const double x, \
		const double xl, const double xm, const double xr) const;
	inline double dpsi(const double x, \
		const double xl, const double xm, const double xr) const;
	inline double ddpsi(const double x, \
		const double xl, const double xm, const double xr) const;
};

template <class TBas>
inline double HSpline3<TBas>::phi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = -2.0 * DL * DL * (DR - 0.5 * HK) / (HK * HK * HK);
	} else if (x == xm) {
		res = 1.0; //for cases like x = xm = xr
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = 2.0 * (DL + 0.5 * HK) * DR * DR / (HK * HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline3<TBas>::dphi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = -2.0 * (DL * (DR - 0.5 * HK) + DL * (DR - 0.5 * HK) + DL * DL) / (HK * HK * HK);
	} else if (x == xm) {
		res = 0.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = 2.0 * (DR * DR + (DL + 0.5 * HK) * DR + (DL + 0.5 * HK) * DR) / (HK * HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline3<TBas>::ddphi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = -4.0 * (DR - 0.5 * HK + 2.0 * DL) / (HK * HK * HK);
	} else if (x == xm) {
		res = 0.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = 4.0 * (DL + 0.5 * HK + 2.0 * DR) / (HK * HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline3<TBas>::psi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = DL * DL * DR / (HK * HK);
	} else if (x == xm) {
		res = 0.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = DL * DR * DR / (HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline3<TBas>::dpsi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = DL * (DR + DR + DL) / (HK * HK);
	} else if (x == xm) {
		res = 1.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = DR * (DR + DL + DL) / (HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline3<TBas>::ddpsi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = 2.0 * (2.0 * DL + DR) / (HK * HK);
	} else if (x == xm) {
		res = 0.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = 2.0 * (DL + 2.0 * DR) / (HK * HK);
	}
	return res;
}

