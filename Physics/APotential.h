#pragma once

#include "cph.h"

template <class T>
class APotential {
public:
	APotential() = default;
	APotential(const APotential<T> &rhs) = delete;
	APotential(APotential<T> &&rhs) = delete;
	APotential<T> & operator=(const APotential<T> &rhs) = delete;
	APotential<T> & operator=(APotential<T> &&rhs) = delete;
	virtual T operator()(const T r) const = 0;
	//returns potential value calculated at rConv*r
	virtual ~APotential(void) = default;
};

template <class T>
class ZeroPot :
	public APotential<T> {
public:
	ZeroPot(void) = default;
	inline T operator()(const T r) const override;
	~ZeroPot(void) = default;
};

template <class T>
inline T ZeroPot<T>::operator()(const T r) const {
	return T();
}