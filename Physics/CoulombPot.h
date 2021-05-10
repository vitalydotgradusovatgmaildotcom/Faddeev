#pragma once

#include "APotential.h"

template <class T>
class CoulombPot :
	public APotential<T> {
public:
	CoulombPot(const double z1z2);
	T operator()(const T r) const override;
	~CoulombPot(void) = default;
protected:
	double z1z2;
};

template <class T>
CoulombPot<T>::CoulombPot(const double z1z2) : z1z2(z1z2) { }

template <class T>
T CoulombPot<T>::operator()(const T r) const {
	//cout << "Potential: r = " << r << "  z1z2 = " << z1z2 << endl;
	return z1z2/r;
}