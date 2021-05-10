#pragma once

#include "APotential.h"

//He - He TTYPT
// potential HFD-B(He), distance in a.u., Value in a.u.e
// Tang, Toennies and Yiu, PRL 74, 1546 (1995)

template <class T>
class HeTTY :
	public APotential<T> {
public:
	HeTTY() = default;
	T operator()(const T r) const override;
	~HeTTY(void) = default;
protected:
	static const double beta;
	static const double d;
	static const double fact[24];
	static const double c[12];
	T brkExp(const Int n, const T br) const;
};

template <class T>
const double HeTTY<T>::beta = 1.3443;
template <class T>
const double HeTTY<T>::d = 7.449;
template <class T>
const double HeTTY<T>::fact[24] = {
	0.10000000000000000e+01,
	0.20000000000000000e+01,
	0.60000000000000000e+01,
	0.24000000000000000e+02,
	0.12000000000000000e+03,
	0.72000000000000000e+03,
	0.50400000000000000e+04,
	0.40320000000000000e+05,
	0.36288000000000000e+06,
	0.36288000000000000e+07,
	0.39916800000000000e+08,
	0.47900160000000000e+09,
	0.62270208000000000e+10,
	0.87178291200000000e+11,
	0.13076743680000000e+13,
	0.20922789888000000e+14,
	0.35568742809600000e+15,
	0.64023737057279940e+16,
	0.12164510040883208e+18,
	0.24329020081766400e+19,
	0.51090942171709424e+20,
	0.11240007277776115e+22,
	0.25852016738885062e+23,
	0.62044840173324357e+24
};

template <class T>
const double HeTTY<T>::c[12] = {
	0.0e0,
	0.0e0,
	1.461e0,
	14.11e0,
	183.5e0,
	3213.49316440087568e0,
	75779.3628279993718e0,
	2406338.80806937255e0,
	102895138.969347134e0,
	5924684085.54217720e0,
	459375463771.044128e0,
	47962592392511.7656e0
};

template <>
double HeTTY<double>::operator()(const double r) const {
	double xau = r; //*0.5291772108/0.529177; ///0.529177;
	double vex = d * pow(xau, 3.5 / beta - 1.0) * exp(-2 * beta*xau);
	double b = 2 * beta - (3.5 / beta - 1) / xau;
	double vdisp = 0.0;
	double br = xau * b;
	double f2i = brkExp(2, br);
	double f2ir = 0.0;
	for (Int i = 3; i <= 12; i++) {
		f2i = f2i - exp(-br)*pow(br, 2 * i - 1)*(1.0 + br / (i*2.0)) / fact[2 * i - 2];
		f2ir = f2i;
		if (abs(f2ir) < 4.0e-14)
			f2ir = 0.0;
		vdisp = vdisp - c[i - 1] * f2ir / pow(xau, 2 * i);
	}
	double res = (vdisp + vex);
	double shortRangeCut = 8.0e-2;
	if (res > shortRangeCut)
		res = shortRangeCut;
	return res;// /(3.1577464e5*3.1669e-6);
}

template <>
Complex HeTTY<Complex>::operator()(const Complex r) const {
	Complex xau = r; //*0.5291772108/0.529177; ///0.529177;
	Complex vex = d * pow(xau, 3.5 / beta - 1.0) * exp(-2 * beta*xau);
	Complex b = 2 * beta - (3.5 / beta - 1) / xau;
	Complex vdisp = zzero;
	Complex br = xau * b;
	Complex f2i = brkExp(2, br);
	Complex f2ir = zzero;
	for (Int i = 3; i <= 12; i++) {
		f2i = f2i - exp(-br)*pow(br, 2 * i - 1)*(1.0 + br / (i*2.0)) / fact[2 * i - 2];
		f2ir = f2i;
		if (abs(f2ir) < 4.0e-14)
			f2ir = 0.0;
		vdisp = vdisp - c[i - 1] * f2ir / pow(xau, 2 * i);
	}
	Complex res = (vdisp + vex);
	double shortRangeCut = 8.0e-2;
	if (res.real() > shortRangeCut) //TODO check condition
		res = shortRangeCut;
	return res;// /(3.1577464e5*3.1669e-6);
}

template<class T>
T HeTTY<T>::brkExp(const Int n, const T br) const {
	T s = T();
	Int kk;
	for (Int i = 1; i <= 2 * n; i++) {
		kk = 2 * n - i + 1;
		s += pow(br, kk) * exp(-br) / fact[kk - 1];
	}
	s += exp(-br);
	return 1.0 - s;
}
