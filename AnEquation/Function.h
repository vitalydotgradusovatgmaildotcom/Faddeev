#pragma once

#include "AFunction.h"

template <class T, Int d> class VectFunction;
template <class T, Int d> class Operator;

template <class T, Int d>
class Function :
	public AFunction<T> {
public:
	using AFunction<T>::coef;
	Function(const shared_ptr<const ADiscretization<T, d>> &discr);
	Function(const Function<T, d> &rhs) = default;
	Function(Function<T, d> &&rhs) = default;
	Function<T, d> & operator=(const Function<T, d> &rhs) = default;
	Function<T, d> & operator=(const AFunction<T> &rhs) override;
	Function<T, d> & operator=(Function<T, d> &&rhs) = default;
	Function<T, d> & operator=(AFunction<T> &&rhs) override;
	Function<T, d> operator*(const T x) const; //calculates (scalar x)*this
	Function<T, d> operator+(const Function<T, d> &u) const; // calculates this + u
	Function<T, d> operator-(const Function<T, d> &u) const; // calculates this - u
	T operator()(const std::array<double, d> &xs) const; //calculates value at x
	T df(const std::array<double, d> &xs) const; //calculates derivative at x
	T d2f(const std::array<double, d> &xs) const; //calculates second derivative at x
	//shared_ptr<const ADiscretization<T, d>> getDiscr() const;
	unique_ptr<AFunction<T>> clone() const override;
	~Function(void);
	friend class Operator<T, d>;
protected:
	Function();
	Function(const shared_ptr<const ADiscretization<T, d>> &discr, Vector<T> &&coef);
	shared_ptr<const ADiscretization<T, d>> discr;
	friend class VectFunction<T, d>;
};

template <class T, Int d>
Function<T, d>::Function() {}

template <class T, Int d>
Function<T, d>::Function( \
	const shared_ptr<const ADiscretization<T, d>> &discr) \
				: discr(discr) {
	coef.resize(discr->getN());
}

template <class T, Int d>
Function<T, d>::Function(const shared_ptr<const ADiscretization<T, d>> &discr, Vector<T> &&coef) \
	: AFunction<T>(std::move(coef)), discr(discr) {}

/*
template <class T, Int d>
Function<T, d>::Function(const Function<T, d> &rhs) \
	: AFunction<T>(rhs), discr(rhs.discr) { }

template <class T, Int d>
Function<T, d>::Function(Function<T, d> &&rhs) \
	: AFunction<T>(std::move(rhs)), discr(std::move(rhs.discr)) {}

template <class T, Int d>
Function<T, d> & Function<T, d>::operator=(const Function<T, d> &rhs) {
	AFunction::operator=(rhs);
	this->discr = rhs.discr;
	return *this;
}*/

template <class T, Int d>
Function<T, d> & Function<T, d>::operator=(const AFunction<T> &rhs) {
	const Function<T, d> &rhs_ = \
		dynamic_cast<const Function<T, d> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

/*template <class T, Int d>
Function<T, d> & Function<T, d>::operator=(Function<T, d> &&rhs) {
	AFunction<T>::operator=(std::move(rhs));
	discr = rhs.discr; rhs.discr.reset();
	return *this;
}*/

template <class T, Int d>
Function<T, d> & Function<T, d>::operator=(AFunction<T> &&rhs) {
	Function<T, d> &&rhs_ = \
		dynamic_cast<Function<T, d> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T, Int d>
Function<T, d> Function<T, d>::operator*(const T x) const {
	Function<T, d> res(*this);
	res *= x;
	return res;
}

template <class T, Int d>
Function<T, d> Function<T, d>::operator+(const Function<T, d> &u) const {
	Function<T, d> res(*this);
	res += u;
	return res;
}

template <class T, Int d>
Function<T, d> Function<T, d>::operator-(const Function<T, d> &u) const {
	Function<T, d> res(*this);
	res -= u;
	return res;
}

template <class T, Int d>
T Function<T, d>::operator()(const array<double, d> &xs) const {
	return discr->f(xs, *this);
}

template <class T, Int d>
T Function<T, d>::df(const array<double, d> &xs) const {
	return discr->df(xs, *this);
}

template <class T, Int d>
T Function<T, d>::d2f(const array<double, d> &xs) const {
	return discr->d2f(xs, *this);
}

/*
template <class T, Int d>
shared_ptr<const ADiscretization<T, d>> Function<T, d>::getDiscr() const {
	return discr;
}*/

template <class T, Int d>
unique_ptr<AFunction<T>> Function<T, d>::clone() const {
	return make_unique<Function<T, d>>(*this);
}

template <class T, Int d>
Function<T, d>::~Function(void) { }