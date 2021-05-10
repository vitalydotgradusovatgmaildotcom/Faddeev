#pragma once

#include "APotential.h"

template <class T>
class CoulombTail :
	public APotential<T> {
public:
	CoulombTail(const double z1z2, const double r0);
	T operator()(const T r) const override;
	~CoulombTail(void) = default;
protected:
	double z1z2;
	double r0; //splitting parameter
};

template <class T>
CoulombTail<T>::CoulombTail(const double z1z2, const double r0) \
		: z1z2(z1z2), r0(r0) {
	//assert(r0 != 0.0);
}

template <>
double CoulombTail<double>::operator()(const double r) const {
	return ( 1.0 - exp( -pow(r/r0, 2.0) ) ) * z1z2 / r;
}

template <>
Complex CoulombTail<Complex>::operator()(const Complex r) const {
	Complex pwr = pow(r / r0, 2.0);
	Complex chi;
	if (pwr.real() > 700.0)
		chi = zzero;
	else
		chi = exp(-pwr);
	return (1.0 - chi) * z1z2 / r;
}
