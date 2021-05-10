#pragma once

#include "Function.h"
#include "VectDiscr.h"

template <class T, Int d> class VectOperator;

template <class T, Int d>
class VectFunction :
	public AFunction<T> {
public:
	using AFunction<T>::coef;
	VectFunction(const shared_ptr<const VectDiscr<T, d>> &vd);
	VectFunction(const VectFunction<T, d> &rhs);
	VectFunction(VectFunction<T, d> &&rhs);
	VectFunction<T, d> & operator=(const VectFunction<T, d> &rhs);
	VectFunction<T, d> & operator=(const AFunction<T> &rhs) override;
	VectFunction<T, d> & operator=(VectFunction<T, d> &&rhs);
	VectFunction<T, d> & operator=(AFunction<T> &&rhs) override;
	VectFunction<T, d> operator*(const T x) const;
	VectFunction<T, d> operator+(const VectFunction<T, d> &u) const;
	VectFunction<T, d> operator-(const VectFunction<T, d> &u) const;
	inline const Function<T, d> & getF(const Int i) const;
	unique_ptr<AFunction<T>> clone() const override;
	//shared_ptr<const VectDiscr<T, d>> getDiscr() const;
	~VectFunction(void);
	friend class VectOperator<T, d>;
protected:
	VectFunction();
	shared_ptr<const VectDiscr<T, d>> vd;
	std::vector<Function<T, d> *> fs;
	Int nFunc;
};

template <class T, Int d>
inline const Function<T, d> & VectFunction<T, d>::getF(const Int i) const {
	return *(fs[i]);
}

template <class T, Int d>
VectFunction<T, d>::VectFunction() {}

template <class T, Int d>
VectFunction<T, d>::VectFunction( \
	const shared_ptr<const VectDiscr<T, d>> &vd) \
		: vd(vd), nFunc(vd->getNFunc()) {
	coef.resize(vd->getN());

	vector<Int> inds; inds.reserve(nFunc);
	Int ind = 0; inds.push_back(ind);
	for (Int k = 1; k < nFunc; k++) {
		ind += (*vd)[k-1].getN();
		inds.push_back(ind);
	}
	vector<Vector<T>> coefs = coef.makePartition(inds);

	fs.resize(nFunc);
	for (Int k = 0; k < nFunc; k++)
		fs[k] = new Function<T, d>(vd->discrs[k], std::move(coefs[k]));
		//fs[k] = make_unique<Function<T, d>>((*vd)[k], std::move(coefs[k]));
		//fs[k] = unique_ptr<Function<T, d>>(new Function<T, d>((*vd)[k], std::move(coefs[k])) );
}

template <class T, Int d>
VectFunction<T, d>::VectFunction<T, d>(const VectFunction<T, d> &rhs) \
		: vd(rhs.vd), nFunc(rhs.nFunc) {
	coef = rhs.coef;
			
	vector<Int> inds; inds.reserve(nFunc);
	Int ind = 0; inds.push_back(ind);
	for (Int k = 1; k < nFunc; k++) {
		ind += (*vd)[k-1].getN();
		inds.push_back(ind);
	}
	vector<Vector<T>> coefs = coef.makePartition(inds);

	fs.resize(nFunc);
	for (Int k = 0; k < nFunc; k++)
		fs[k] = new Function<T, d>(vd->discrs[k], std::move(coefs[k]));
		//fs[k] = make_unique<Function<T, d>>((*vd)[k], std::move(coefs[k]));
		//fs[k] = unique_ptr<Function<T, d>>( new Function<T, d>((*vd)[k], std::move(coefs[k])) );
}

template <class T, Int d>
VectFunction<T, d>::VectFunction<T, d>(VectFunction<T, d> &&rhs) \
		: vd(std::move(rhs.vd)), nFunc(rhs.nFunc) {
	coef = std::move(rhs.coef);
	rhs.nFunc = 0;

	fs.resize(nFunc);
	for (Int k = 0; k < nFunc; k++)
		fs[k] = std::move(rhs.fs[k]);

	rhs.fs.clear(); rhs.fs.shrink_to_fit();
}

template<class T, Int d>
VectFunction<T, d> & VectFunction<T, d>::operator=(const VectFunction<T, d> &rhs) {
	AFunction<T>::operator=(rhs);
	vd = rhs.vd;
	nFunc = rhs.nFunc;

	for (Int k = 0; k < fs.size(); k++)
		if (fs[k] != nullptr)
			delete fs[k];
		//fs[k].reset();

	vector<Int> inds; inds.reserve(nFunc);
	Int ind = 0; inds.push_back(ind);
	for (Int k = 1; k < nFunc; k++) {
		ind += (*vd)[k-1].getN();
		inds.push_back(ind);
	}
	vector<Vector<T>> coefs = coef.makePartition(inds);

	fs.resize(nFunc);
	for (Int k = 0; k < nFunc; k++)
		fs[k] = new Function<T, d>(vd->discrs[k], std::move(coefs[k]));
		//fs[k] = make_unique<Function<T, d>>((*vd)[k], std::move(coefs[k]));
		//fs[k] = unique_ptr<Function<T, d>>( new Function<T, d>((*vd)[k], std::move(coefs[k])) );

	return *this;
}

template<class T, Int d>
VectFunction<T, d> & VectFunction<T, d>::operator=(const AFunction<T> &rhs) {
	const VectFunction<T, d> &rhs_ = \
		dynamic_cast<const VectFunction<T, d> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template<class T, Int d>
VectFunction<T, d> & VectFunction<T, d>::operator=(VectFunction<T, d> &&rhs) {
	AFunction<T>::operator=(std::move(rhs));
	vd = rhs.vd; rhs.vd.reset();
	nFunc = rhs.nFunc; rhs.nFunc = 0;

	for (Int k = 0; k < fs.size(); k++)
		if (fs[k] != nullptr)
			delete fs[k];
		//fs[k].reset();

	fs.resize(nFunc);
	for (Int k = 0; k < fs.size(); k++)
		fs[k] = std::move(rhs.fs[k]);
	rhs.fs.clear(); rhs.fs.shrink_to_fit();

	return *this;
}

template <class T, Int d>
VectFunction<T, d> & VectFunction<T, d>::operator=(AFunction<T> &&rhs) {
	VectFunction<T, d> &&rhs_ = \
		dynamic_cast<VectFunction<T, d> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T, Int d>
VectFunction<T, d> VectFunction<T, d>::operator*(const T x) const {
	VectFunction<T, d> res(*this);
	res *= x;
	return res;
}

template <class T, Int d>
VectFunction<T, d> VectFunction<T, d>::operator+(const VectFunction<T, d> &u) const {
	VectFunction<T, d> res(*this);
	res += u;
	return res;
}

template <class T, Int d>
VectFunction<T, d> VectFunction<T, d>::operator-(const VectFunction<T, d> &u) const {
	VectFunction<T, d> res(*this);
	res -= u;
	return res;
}

template <class T, Int d>
unique_ptr<AFunction<T>> VectFunction<T, d>::clone() const {
	return make_unique<VectFunction<T, d>>(*this);
}

/*
template <class T, Int d>
shared_ptr<const VectDiscr<T, d>> VectFunction<T, d>::getDiscr() const {
	return vd;
}*/

template <class T, Int d>
VectFunction<T, d>::~VectFunction(void) {
	for (Int k = 0; k < nFunc; k++)
		if (fs[k] != nullptr)
			delete fs[k];
}
