#pragma once

#include "cph.h"
#include "ABasis.h"

template <class T>
class AnOperator {
public:
	AnOperator();
	AnOperator(const AnOperator<T> &rhs) = delete;
	AnOperator(AnOperator<T> &&rhs) = delete;
	AnOperator<T> & operator=(const AnOperator<T> &rhs) = delete;
	AnOperator<T> & operator=(AnOperator<T> &&rhs) = delete;
	//virtual AFunction<T> && operator*(const AFunction<T> &u) const = 0;
	virtual void times(const AFunction<T> &u, AFunction<T> &res) = 0;
	virtual void solve(AFunction<Complex> &rhssol) = 0;
	Int getRank() const;
	virtual ~AnOperator(void);
protected:
	bool changed = false;
	Int rank;
};

template <class T>
AnOperator<T>::AnOperator() { }

/*template<class T>
AnOperator<T> & AnOperator<T>::operator=(const AnOperator<T> &rhs) {
	this->changed = rhs.changed;
	this->rank = rhs.rank;
	return *this;
}*/

template <class T>
Int AnOperator<T>::getRank() const {
	return rank;
}

template <class T>
AnOperator<T>::~AnOperator(void) { }
