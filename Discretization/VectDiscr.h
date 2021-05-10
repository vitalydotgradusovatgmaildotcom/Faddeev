#pragma once

#include "ADiscretization.h"

template <class T, Int d> class VectFunction;

template <class T, Int d>
class VectDiscr {
public:
	VectDiscr(const vector<shared_ptr<const ADiscretization<T, d>>> &discrs);
	Int getNFunc() const;
	const ADiscretization<T, d> & operator[](const Int i) const;
	Int getRaw(const Int iFunc, const std::array<Int, d> &is) const;
	void getInds(const Int i, Int &iFunc, std::array<Int, d> &is) const;
	Int getN() const;
	~VectDiscr(void);
protected:
	vector<shared_ptr<const ADiscretization<T, d>>> discrs;
	Int nFunc;
	Int n = 0;
	vector<Int> iStart;
	friend class VectFunction<T, d>;
};

template <class T, Int d>
VectDiscr<T, d>::VectDiscr(const vector<shared_ptr<const ADiscretization<T, d>>> &discrs) \
	: nFunc(discrs.size()) {
	this->discrs.reserve(nFunc);
	iStart.reserve(nFunc);
	for (Int k = 0; k < nFunc; k++) {
		this->discrs.push_back(discrs[k]);
		iStart.push_back(n);
		n += discrs[k]->getN();
	}
}

template <class T, Int d>
Int VectDiscr<T, d>::getNFunc() const {
	return nFunc;
}

template <class T, Int d>
const ADiscretization<T, d> & VectDiscr<T, d>::operator[](const Int i) const {
	return *(discrs[i]);
}

template <class T, Int d>
Int VectDiscr<T, d>::getRaw(const Int iFunc, const std::array<Int, d> &is) const {
	assert(iFunc < nFunc);
	return iStart[iFunc] + discrs[iFunc]->getRaw(is);
}

template <class T, Int d>
void VectDiscr<T, d>::getInds(const Int i, Int &iFunc, std::array<Int, d> &is) const {
	assert(i < n);
	static vector<Int>::const_iterator up;
	up = std::upper_bound(iStart.begin(), iStart.end(), i);
	iFunc = --up - iStart.begin();
	discrs[iFunc]->getInds(i-*up, is);
}

template <class T, Int d>
Int VectDiscr<T, d>::getN() const {
	return n;
}

template <class T, Int d>
VectDiscr<T, d>::~VectDiscr(void) { }
