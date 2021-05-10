#pragma once

#include "AnOperator.h"
#include "VectDiscr.h"
#include "VectFunction.h"

template <class T, Int d>
class VectOperator :
	public AnOperator<T> {
public:
	VectOperator(const shared_ptr<const VectDiscr<T, d>> &vdiscr);
	shared_ptr<const VectDiscr<T, d>> getDiscr() const;
	~VectOperator(void) = default;
protected:
	using AnOperator<T>::rank;
	shared_ptr<const VectDiscr<T, d>> vdiscr;
};

template <class T, Int d>
VectOperator<T, d>::VectOperator(\
	const shared_ptr<const VectDiscr<T, d>> &vdiscr) \
		: vdiscr(vdiscr) {
	rank = vdiscr->getN();
}


template <class T, Int d>
shared_ptr<const VectDiscr<T, d>> VectOperator<T, d>::getDiscr() const {
	return vdiscr;
}

//template <class T, Int d>
//VectOperator<T, d>::~VectOperator(void) { }
