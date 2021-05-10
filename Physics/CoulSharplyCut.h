#pragma once

#include "APotential.h"

template <class T>
class CoulSharplyCut :
	public APotential<T> {
public:
	CoulSharplyCut(const double z1z2, const double R);
	T operator()(const T r) const override;
	~CoulSharplyCut(void) = default;
protected:
	double z1z2;
	double R;
};

template <class T>
CoulSharplyCut<T>::CoulSharplyCut(const double z1z2, const double R) : \
	z1z2(z1z2), R(R) { }

template <class T>
T CoulSharplyCut<T>::operator()(const T r) const {
	if (abs(r) >= R)
		return T();
	return z1z2 / r;
}

