#pragma once

#include "APotential.h"

template <class T>
class KTTY :
	public APotential<T> {
public:
	KTTY(double a, double b1, double b2, \
		double c6, double c8, double c10);
	T operator()(const T r) const override;
	~KTTY() = default;
protected:
	double a, b1, b2;
	std::vector<double> c = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	static const double fact[24];
	T brkExp(const Int n, const T br) const;
};

template <class T>
const double KTTY<T>::fact[24] = {
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
KTTY<T>::KTTY(double a, double b1, double b2, \
	double c6, double c8, double c10) : \
				a(a), b1(b1), b2(b2) {
	c[2] = c6; c[3] = c8; c[4] = c10;
	for (Int i = 5; i <= 7; i++)
		c[i] = pow(c[i - 1] / c[i - 2], 3.0)*c[i - 3];
}

template <>
double KTTY<double>::operator()(const double r) const {
	double xau = r;
	double vex = a * exp(-(b1 + b2 * xau)*xau);
	double vdisp = 0.0;
	double br = (b1 + 2 * xau*b2)*xau;
	double f2i = brkExp(2, br);
	double f2ir = 0.0;
	for (Int i = 3; i <= 8; i++) {
		f2i = f2i - exp(-br)*pow(br, 2 * i - 1)*(1.0 + br / (i*2.0)) / fact[2 * i - 2];
		f2ir = f2i;
		if (abs(f2ir) < 4.0e-13) f2ir = 0.0;
		vdisp = vdisp - c[i - 1] * f2ir / pow(xau, 2 * i);
	}
	double res = (vdisp + vex);
	double shortRangeCut = 1.0e-1;
	if (res > shortRangeCut)
		res = shortRangeCut;
	return res;
}

template <>
Complex KTTY<Complex>::operator()(const Complex r) const {
	Complex xau = r;
	Complex vex = a * exp(-(b1 + b2 * xau)*xau);
	Complex vdisp = 0.0;
	Complex br = (b1 + 2.0 * xau*b2)*xau;
	Complex f2i = brkExp(2, br);
	Complex f2ir = 0.0;
	for (Int i = 3; i <= 8; i++) {
		f2i = f2i - exp(-br)*pow(br, 2 * i - 1)*(1.0 + br / (i*2.0)) / fact[2 * i - 2];
		f2ir = f2i;
		if (abs(f2ir) < 4.0e-13) f2ir = 0.0;
		vdisp = vdisp - c[i - 1] * f2ir / pow(xau, 2 * i);
	}
	Complex res = (vdisp + vex);
	double shortRangeCut = 1.0e-1;
	if (res.real() > shortRangeCut) //TODO check condition
		res = shortRangeCut;
	return res;
}

template <class T>
T KTTY<T>::brkExp(const Int n, const T br) const {
	T s = 0.0;
	Int kk;
	for (Int i = 1; i <= 2 * n; i++) {
		kk = 2 * n - i + 1;
		s += pow(br, kk)*exp(-br) / fact[kk - 1];
	}
	s += exp(-br);
	return 1.0 - s;
}