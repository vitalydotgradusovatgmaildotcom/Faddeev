#pragma once

#include "Vector.h"
template<class T, Int d> class ADiscretization;

template <class T>
class AFunction {
public:
	Vector<T> coef;
	AFunction();
	AFunction(const AFunction<T> &rhs) = default;
	AFunction(AFunction<T> &&rhs) = default;
	virtual AFunction<T> & operator=(const AFunction<T> &rhs); //= default;
	virtual AFunction<T> & operator=(AFunction<T> &&rhs); //= default;
	// =default does not work - bug???
	//virtual AFunction<T> & operator*(const T x) const = 0; //calculates (scalar x)*this
	//virtual AFunction<T> & operator+(const AFunction<T> &u) const = 0; // calculates this + u
	//virtual AFunction<T> & operator-(const AFunction<T> &u) const = 0; // calculates this - u
	AFunction<T> & operator*=(const T x); //multiplies this by x
	AFunction<T> & operator+=(const AFunction<T> &u); //adds u to this
	AFunction<T> & operator-=(const AFunction<T> &u); //subtracts u from this
	void setZero(); //sets all the coefficients eq to zero
	Int size() const;
	inline T & operator[](const Int i);
	inline const T & operator[](const Int i) const;
	virtual unique_ptr<AFunction<T>> clone() const = 0;
	~AFunction(void);
protected:
	AFunction(Vector<T> &&coef);
};

template <class T>
inline T & AFunction<T>::operator[](const Int i) {
	return coef[i];
}

template <class T>
inline const T & AFunction<T>::operator[](const Int i) const {
	return coef[i];
}

template <class T>
AFunction<T>::AFunction() { }

template <class T>
AFunction<T>::AFunction(Vector<T> &&coef) : coef(std::move(coef)) {}

/*
template <class T>
AFunction<T>::AFunction(const AFunction<T> &rhs) : coef(rhs.coef) { }
*/

template <class T>
AFunction<T> & AFunction<T>::operator=(const AFunction<T> &rhs) {
	this->coef = rhs.coef;
	return *this;
}

template <class T>
AFunction<T> & AFunction<T>::operator=(AFunction<T> &&rhs) {
	coef = std::move(rhs.coef);
	return *this;
}

template <class T>
AFunction<T> & AFunction<T>::operator*=(const T x) {
	this->coef *= x;
	return *this;
}

template <class T>
AFunction<T> & AFunction<T>::operator+=(const AFunction<T> &u) {
	this->coef += u.coef;
	return *this;
}

template <class T>
AFunction<T> & AFunction<T>::operator-=(const AFunction<T> &u) {
	this->coef -= u.coef;
	return *this;
}

template <class T>
void AFunction<T>::setZero() {
	coef.setZero();
}

template <class T>
Int AFunction<T>::size() const {
	return coef.size();
}

template <class T>
AFunction<T>::~AFunction(void) { }
