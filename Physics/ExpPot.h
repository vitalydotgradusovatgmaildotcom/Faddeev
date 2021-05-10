#pragma once

#include "APotential.h"

//A potential of the form V_0 exp( -r/a )
template <class T>
class ExpPot :
	public APotential<T> {
public:
	ExpPot(const double v0, const double a);
	T operator()(const T r) const override;
	~ExpPot() = default;
protected:
	double v0;
	double a;
};

template <class T>
ExpPot<T>::ExpPot(const double v0, const double a) : v0(v0), a(a) { }

template <class T>
T ExpPot<T>::operator()(const T r) const {
	return v0*exp(-r/a);
}
