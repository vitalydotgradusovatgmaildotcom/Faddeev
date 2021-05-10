#pragma once

#include "AHermitSpline.h"

template <class TBas>
class HSpline5 :
	public AHermitSpline<TBas> {
public:
	HSpline5(const Grid &gr);
	HSpline5(const Grid &gr, \
		const enum bc bcTypeLeft, const enum bc bcTypeRight, \
		const TBas bcLeft, const TBas bcRight);
	TBas f(const double x, const Int n) const override;
	TBas d(const double x, const Int n) const override;
	TBas dd(const double x, const Int n) const override;
	TBas f(const double x, const AFunction<double> &func) const override;
	Complex f(const double x, const AFunction<Complex> &func) const override;
	TBas d(const double x, const AFunction<double> &func) const override;
	Complex d(const double x, const AFunction<Complex> &func) const override;
	TBas dd(const double x, const AFunction<double> &func) const override;
	Complex dd(const double x, const AFunction<Complex> &func) const override;
	~HSpline5(void) = default;
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
	inline double chi(const double x, \
		const double xl, const double xm, const double xr) const;
	inline double dchi(const double x, \
		const double xl, const double xm, const double xr) const;
	inline double ddchi(const double x, \
		const double xl, const double xm, const double xr) const;
};

template <class TBas>
inline double HSpline5<TBas>::phi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = DL * DL * DL * (6 * DR * DR - 3 * DR * HK + HK * HK) / (HK * HK * HK * HK * HK);
	} else if (x == xm) {
		res = 1.0; //for cases like x = xm = xr
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = DR * DR * DR * (-6 * DL * DL - 3 * DL * HK - HK * HK) / (HK * HK * HK * HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline5<TBas>::dphi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = (3.0 * DL * DL * (6.0 * DR * DR - 3.0 * DR * HK + HK * HK) + DL * DL * DL * (12.0 * DR - 3.0 * HK)) / (HK * HK * HK * HK * HK);
	} else if (x == xm) {
		res = 0.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = (3.0 * DR * DR * (-6.0 * DL * DL - 3.0 * DL * HK - HK * HK) + DR * DR * DR * (-12.0 * DL - 3.0 * HK)) / (HK * HK * HK * HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline5<TBas>::ddphi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = (6 * DL * (6 * DR * DR - 3 * DR * HK + HK * HK) + 6 * DL * DL * (12 * DR - 3 * HK) + DL * DL * DL * 12) / (HK * HK * HK * HK * HK);
	} else if (x == xm) {
		res = 0.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = (6 * DR * (-6 * DL * DL - 3 * DL * HK - HK * HK) + 6 * DR * DR * (-12 * DL - 3 * HK) + DR * DR * DR * (-12)) / (HK * HK * HK * HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline5<TBas>::psi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = DL * DL * DL * (DL - 4 * DR) * DR / (HK * HK * HK * HK);
	} else if (x == xm) {
		res = 0.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = DR * DR * DR * (DR - 4 * DL) * DL / (HK * HK * HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline5<TBas>::dpsi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = DL * DL * (3.0 * (DL - 4.0 * DR) * DR - 3.0 * DL * DR + DL * (DL - 4.0 * DR)) / (HK * HK * HK * HK);
	} else if (x == xm) {
		res = 1.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = DR * DR * (3.0 * (DR - 4.0 * DL) * DL - 3.0 * DR * DL + DR * (DR - 4.0 * DL)) / (HK * HK * HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline5<TBas>::ddpsi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = -6 * DL * DR * (4 * DR + 6 * DL) / (HK * HK * HK * HK);
	} else if (x == xm) {
		res = 0.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = -6 * DL * DR * (4 * DL + 6 * DR) / (HK * HK * HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline5<TBas>::chi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = DL * DL * DL * DR * DR * 0.5 / (HK * HK * HK);
	} else if (x == xm) {
		res = 0.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = -DR * DR * DR * DL * DL * 0.5 / (HK * HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline5<TBas>::dchi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = DL * DL * DR * (3.0 * DR + 2.0 * DL) * 0.5 / (HK * HK * HK);
	} else if (x == xm) {
		res = 0.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = -DR * DR * DL * (3.0 * DL + 2.0 * DR) * 0.5 / (HK * HK * HK);
	}
	return res;
}

template <class TBas>
inline double HSpline5<TBas>::ddchi(const double x, \
		const double xl, const double xm, const double xr) const {
	
	if (x < xl || x > xr) return 0.0;

	double res, HK, DL, DR;
	if ( x < xm ) { //x in [xl, xm)
		HK = xm - xl; DL = x - xl; DR = x - xm;
		res = DL * (DL * DL + 3.0 * DR * (DR + 2.0 * DL)) / (HK * HK * HK);
	} else if (x == xm) {
		res = 1.0;
	} else { //x in (xm, xr]
		HK = xr - xm; DL = x - xm; DR = x - xr;
		res = -DR * (DR * DR + 3.0 * DL * (DL + 2.0 * DR)) / (HK * HK * HK);
	}
	return res;
}
