#pragma once

#include "AMatrix.h"
#include "GenMatrix.h"

//Tensor product of non square matrices
template<class T, Int d>
class TensorProd :
	public AMatrix<T> {
public:
	std::array<std::shared_ptr<AMatrix<T>>, d> matrs;
	std::array<Int, d> ns;
	std::array<Int, d> ms;
	TensorProd(const std::array<std::shared_ptr<AMatrix<T>>, d> &matrs);
	TensorProd(const TensorProd<T, d> &rhs) = default;
	TensorProd(TensorProd<T, d> &&rhs) = default;
	TensorProd<T, d> & operator=(const TensorProd<T, d> &rhs);
	TensorProd<T, d> & operator=(const AMatrix<T> &rhs) override;
	TensorProd<T, d> & operator=(TensorProd<T, d> &&rhs);
	TensorProd<T, d> & operator=(AMatrix<T> &&rhs) override;
	void fill(const T val) override;
	void resize(const Int nnew) override;
	void solve(Vector<T> & rhssol) override;
	void print() const override;
	void write(const string &filename) const override;
	TensorProd<T, d> & operator*=(const T c) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	double sizeGb() const override;
	GenMatrix<T> toGenMatrix() const; //auxiliary test function,
	//computational cost is BIG!
	~TensorProd(void) = default;
protected:
	using AMatrix<T>::n;
	using AMatrix<T>::m;
	using AMatrix<T>::garbage;
};

template <class T, Int d>
TensorProd<T, d>::TensorProd(const std::array<std::shared_ptr<AMatrix<T>>, d> &matrs) : matrs(matrs) {
	assert( std::find(matrs.begin(), matrs.end(), nullptr) == matrs.end() );
	n = 1; m = 1;
	for (Int k = 0; k < matrs.size(); k++) {
		ns[k] = matrs[k]->nrows();
		ms[k] = matrs[k]->ncols();
		n *= ns[k]; m *= ms[k];
	}
}

template <class T, Int d>
TensorProd<T, d> & TensorProd<T, d>::operator=(const TensorProd<T, d> &rhs) {
	if (this != &rhs) {
		n = rhs.n; m = rhs.m;
		matrs = rhs.matrs;
		ns = rhs.ns; ms = rhs.ms;
	}
	return *this;
}

template <class T, Int d>
TensorProd<T, d> & TensorProd<T, d>::operator=(const AMatrix<T> &rhs) {
	const TensorProd<T, d> &rhs_ = dynamic_cast<const TensorProd<T, d> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T, Int d>
TensorProd<T, d> & TensorProd<T, d>::operator=(TensorProd<T, d> &&rhs) {
	n = rhs.n; rhs.n = 0;
	m = rhs.m; rhs.m = 0;
	matrs = std::move(rhs.matrs);
	ns = std::move(rhs.ns);
	ms = std::move(rhs.ms);
	return *this;
}

template <class T, Int d>
TensorProd<T, d> & TensorProd<T, d>::operator=(AMatrix<T> &&rhs) {
	TensorProd<T, d> && rhs_ = dynamic_cast<TensorProd<T, d> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T, Int d>
void TensorProd<T, d>::fill(const T val) {
	assert(false);
}

template <class T, Int d>
void TensorProd<T, d>::resize(const Int nnew) {
	assert(false);
}

template <class T, Int d>
void TensorProd<T, d>::solve(Vector<T> &rhssol) {
	assert(n == m);
	//works only for square multipliers!
	for (Int k = 0; k < matrs.size(); k++)
		assert(ns[k] == ms[k]);
	Int size = rhssol.size();
	assert(size == n);
	Int i, i1block;
	Int nblock = n, sblock = 1;
	Int nsubblock, ssubblock;
	Vector<T> tmp;
	for (Int k = d-1; k >= 0; k--) {
		tmp.resize(ns[k]);
		nsubblock = ns[k]; ssubblock = sblock;
		nblock /= ns[k]; sblock *= ns[k];
		i1block = 0;
		for (Int kblock = 0; kblock < nblock; kblock++) {
			for (Int isub = 0; isub < ssubblock; isub++) {
				i = i1block + isub;
				for ( Int ksubbl = 0; ksubbl < nsubblock; ksubbl++ ) {
					tmp[ksubbl] = rhssol[i];
					i += ssubblock;
				}
				
				matrs[k]->solve(tmp);
				
				i = i1block + isub;
				for ( Int ksubbl = 0; ksubbl < nsubblock; ksubbl++ ) {
					rhssol[i] = tmp[ksubbl];
					i += ssubblock;
				}
			}
			i1block += sblock;
		}
	}
}

template <class T, Int d>
void TensorProd<T, d>::print() const {
	cout << "Tensor product of matrices:" << endl;
	for (Int k = 0; k < matrs.size(); k++) {
		cout << "Matrix " << k << ":" << endl;
		matrs[k]->print();
	}
}

template <class T, Int d>
void TensorProd<T, d>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	f << "Tensor product of matrices:" << endl;
	for (Int k = 0; k < matrs.size(); k++) {
		f << "Matrix " << k << ":" << endl;
		matrs[k]->write(filename);
	}
	f.close();
}

template <class T, Int d>
TensorProd<T, d> & TensorProd<T, d>::operator*=(const T c) {
	*(matrs[0]) *= c;
	return *this;
}

template <class T, Int d>
Vector<T> TensorProd<T, d>::operator*(const Vector<T> &vec) const {
	Int size = vec.size();
	assert(size == m);
	Int i_, i, i1block_, i1block;
	Int nblock = m, sblock;
	Int nsubblock, ssubblock;
	Int sizeNew = m, sSubBlNew = 1;
	Vector<T> tmp;
	Vector<T> res(vec);
	for (Int k = d-1; k >= 0; k--) {
		size = sizeNew;
		nblock /= ms[k]; nsubblock = ms[k];
		sblock = size/nblock; ssubblock = sSubBlNew;
		sizeNew = (long long int)size * ns[k] / ms[k];
		sSubBlNew = sizeNew / nblock;
		Vector<T> res_(std::move(res));
		res.resize(sizeNew);
		i1block_ = 0; i1block = 0;
		for (Int kblock = 0; kblock < nblock; kblock++) {
			for (Int isub = 0; isub < ssubblock; isub++) {
				tmp.resize(ms[k]); //no resize if ns[k] = ms[k]
				i_ = i1block_ + isub;
				for ( Int ksubbl = 0; ksubbl < ms[k]; ksubbl++ ) {
					tmp[ksubbl] = res_[i_];
					i_ += ssubblock;
				}
				
				tmp = *(matrs[k]) * tmp; //tmp changes size to ns[k]
				
				i = i1block + isub;
				for ( Int ksubbl = 0; ksubbl < ns[k]; ksubbl++ ) {
					res[i] = tmp[ksubbl];
					i += ssubblock;
				}
			}
			i1block_ += sblock;
			i1block += sSubBlNew;
		}
	}

	return res;
}

template <class T, Int d>
double TensorProd<T, d>::sizeGb() const {
	assert(false);
}

template <class T, Int d>
GenMatrix<T> TensorProd<T, d>::toGenMatrix() const {
	GenMatrix<T> res(n, m);

	Vector<T> tmp(m);
	for (Int j = 0; j < m; j++) {
		tmp.fill(T()); tmp[j] = 1.0;
		Vector<T> tmp_(*this * tmp);
		for (Int i = 0; i < n; i++)
			res[i][j] = tmp_[i];
	}

	return res;
}
