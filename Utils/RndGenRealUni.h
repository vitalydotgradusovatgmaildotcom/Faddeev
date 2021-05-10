#pragma once

#include "cph.h"
#include <random>

template <class T>
class RndGenRealUni {
public:
	RndGenRealUni(const double a, const double b);
	RndGenRealUni(const RndGenRealUni &rhs) = delete;
	RndGenRealUni(RndGenRealUni &&rhs) = delete;
	RndGenRealUni & operator=(const RndGenRealUni &rhs) = delete;
	RndGenRealUni & operator=(RndGenRealUni &&rhs) = delete;
	void fillArr(T *arr, const Int n);
	T generate();
	~RndGenRealUni() = default;
protected:
	std::mt19937 gen;
	std::uniform_real_distribution<double> distr;
};

template <class T>
RndGenRealUni<T>::RndGenRealUni(const double a, const double b) \
	: gen(500), distr(a, b) { }

template <>
void RndGenRealUni<double>::fillArr(double *arr, const Int n) {
	for (Int ind = 0; ind < n; ind++)
		arr[ind] = distr(gen);
}

template <>
void RndGenRealUni<Complex>::fillArr(Complex *arr, const Int n) {
	for (Int ind = 0; ind < n; ind++)
		arr[ind] = Complex(distr(gen), distr(gen));
}

template <>
double RndGenRealUni<double>::generate() {
	return distr(gen);
}

template <>
Complex RndGenRealUni<Complex>::generate() {
	return Complex(distr(gen), distr(gen));
}