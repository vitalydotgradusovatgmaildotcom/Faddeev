#pragma once

#include "APotential.h"

template <class T>
class ScaledPot :
	public APotential<T> {
public:
	ScaledPot(const shared_ptr<APotential<T>> &pot, \
		const double eConv, const double rConv);
	T operator()(const T r) const override;
	~ScaledPot() = default;
protected:
	shared_ptr<const APotential<T>> pot;
	double eConv;
	double rConv;
};

template <class T>
ScaledPot<T>::ScaledPot( \
	const shared_ptr<APotential<T>> &pot, const double eConv,
		const double rConv) : pot(pot), eConv(eConv), rConv(rConv) { }

template <class T>
T ScaledPot<T>::operator()(const T r) const {
	return eConv * (*pot)(r*rConv);
}