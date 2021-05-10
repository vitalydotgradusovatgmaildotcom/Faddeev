#pragma once

#include "AnOperator.h"
#include "Function.h"
#include "ADiscretization.h"

template <class T, Int d>
class Operator :
	public AnOperator<T> {
public:
	Operator(const shared_ptr<const ADiscretization<T, d>> &discr);
	//Function<T, d> & operator*(const AFunction<T> &u) const override = 0;
	//shared_ptr<const ADiscretization<T, d>> getDiscr() const;
	~Operator(void);
protected:
	using AnOperator<T>::rank;
	shared_ptr<const ADiscretization<T, d>> discr;
//	Function<T, d> && result;
};

template <class T, Int d>
Operator<T, d>::Operator(const shared_ptr<const ADiscretization<T, d>> &discr) \
				: discr(discr) {
	rank = discr->getN();
}

/*
template <class T, Int d>
shared_ptr<const ADiscretization<T, d>> Operator<T, d>::getDiscr() const {
	return discr;
}*/

template <class T, Int d>
Operator<T, d>::~Operator(void) { }
