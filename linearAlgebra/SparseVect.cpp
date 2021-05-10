#include "SparseVect.h"

template class SparseVect<double>;
template class SparseVect<Complex>;

template <class T>
SparseVect<T>::SparseVect(const Int n, const Int nz) \
		: n(n), nz(nz) {
	assert(nz >= 0); assert(n >= nz);
	v = nz > 0 ? (T *)MALLOC(nz*sizeof(T)) : nullptr;
	inds = nz > 0 ? (Int *)MALLOC(nz*sizeof(Int)) : nullptr;
}

template <>
SparseVect<double>::SparseVect(const Vector<double> &u, \
	const std::vector<Int> &inds_vec) \
		: n(u.size()), nz(inds_vec.size()) {
	assert(n >= nz);
	inds = nz > 0 ? (Int *)MALLOC(nz*sizeof(Int)) : nullptr;
	for (Int i = 0; i < nz; i++)
		inds[i] = inds_vec[i];
	v = nz > 0 ? (double *)MALLOC(nz*sizeof(double)) : nullptr;
	cblas_dgthr(nz, u.v, v, inds);
}

template <>
SparseVect<Complex>::SparseVect(const Vector<Complex> &u, \
	const std::vector<Int> &inds_vec) \
		: n(u.size()), nz(inds_vec.size()) {
	assert(n >= nz);
	inds = nz > 0 ? (Int *)MALLOC(nz*sizeof(Int)) : nullptr;
	for (Int i = 0; i < nz; i++)
		inds[i] = inds_vec[i];
	v = nz > 0 ? (Complex *)MALLOC(nz*sizeof(Complex)) : nullptr;
	cblas_zgthr(nz, u.v, v, inds);
}

template <class T>
SparseVect<T>::SparseVect(const SparseVect<T> &rhs) \
								: n(rhs.n), nz(rhs.nz) {
	v = nz > 0 ? (T *)MALLOC(nz*sizeof(T)) : nullptr;
	inds = nz > 0 ? (Int *)MALLOC(nz*sizeof(Int)) : nullptr;
	for (Int i = 0; i < nz; i++) {
		v[i] = rhs.v[i];
		inds[i] = rhs.inds[i];
	}
}

template <class T>
SparseVect<T>::SparseVect(SparseVect<T> &&rhs) \
						: n(rhs.n), nz(rhs.nz) {
	rhs.n = 0; rhs.nz = 0;
	v = rhs.v; rhs.v = nullptr;
	inds = rhs.inds; rhs.inds = nullptr;
}

template <class T>
SparseVect<T> & SparseVect<T>::operator=(const SparseVect<T> &rhs) {
	this->n = rhs.n;
	if (this->nz != rhs.nz) {
		if (v != nullptr) FREE(v);
		if (inds != nullptr) FREE(inds);
		this->nz = rhs.nz;
		v = nz > 0 ? (T *)MALLOC(nz*sizeof(T)) : nullptr;
		inds = nz > 0 ? (Int *)MALLOC(nz*sizeof(Int)) : nullptr;
	}
	for (Int i = 0; i < nz; i++) {
		this->v[i] = rhs.v[i];
		this->inds[i] = rhs.inds[i];
	}
	return *this;
}

template <class T>
SparseVect<T> & SparseVect<T>::operator=(SparseVect<T> &&rhs) {
	if (v != nullptr) FREE(v);
	if (inds != nullptr) FREE(inds);
	v = rhs.v; rhs.v = nullptr;
	inds = rhs.inds; rhs.inds = nullptr;
	n = rhs.n; rhs.n = 0;
	nz = rhs.nz; rhs.nz = 0;
	return *this;
}

template <>
double SparseVect<double>::scal(const Vector<double> &u) const {
	return cblas_ddoti(nz, v, inds, u.v);
}

template <>
Complex SparseVect<Complex>::scal(const Vector<Complex> &u) const {
	Complex res;
	cblas_zdotci_sub(nz, v, inds, u.v, &res);
	return res;
}

template <>
double SparseVect<double>::dotp(const Vector<double> &u) const {
	return scal(u);
}

template <>
Complex SparseVect<Complex>::dotp(const Vector<Complex> &u) const {
	Complex res;
	cblas_zdotui_sub(nz, v, inds, u.v, &res);
	return res;
}

template <>
double SparseVect<double>::norm() const {
	return cblas_dnrm2(nz, v, 1);
}

template <>
double SparseVect<Complex>::norm() const {
	return cblas_dznrm2(nz, v, 1);
}

template <class T>
void SparseVect<T>::normalize() {
	double norma = norm();
	assert(norma != 0.0);
	*this *= 1.0 / norma;
}

template <class T>
SparseVect<T> SparseVect<T>::operator*(const T x) const {
	SparseVect<T> res(n, nz);
	for (Int i = 0; i < nz; i++) {
		res.v[i] = x*this->v[i];
		res.inds[i] = inds[i];
	}
	return res;
}

template <>
Vector<double> SparseVect<double>::operator+(const Vector<double> &u) const {
	assert(n == u.size());
	Vector<double> res(u);
	cblas_daxpyi(nz, 1.0, v, inds, u.v);
	return res;
}

template <>
Vector<Complex> SparseVect<Complex>::operator+(const Vector<Complex> &u) const {
	assert(n == u.size());
	Vector<Complex> res(u);
	cblas_zaxpyi(nz, &zone, v, inds, u.v);
	return res;
}

template <>
SparseVect<double> & SparseVect<double>::operator*=(const double x) {
	cblas_dscal(nz, x, this->v, 1);
	return *this;
}

template <>
SparseVect<Complex> & SparseVect<Complex>::operator*=(const Complex x) {
	cblas_zscal(nz, &x, this->v, 1);
	return *this;
}

template <class T>
void SparseVect<T>::resize(const Int nnew, const Int nznew) {
	assert(nznew >= 0); assert(nnew >= nznew);
	n = nnew;
	if (nz != nznew) {
		if (v != nullptr) FREE(v);
		if (inds != nullptr) FREE(inds);
		nz = nznew;
		v = nz > 0 ? (T *)MALLOC(nz*sizeof(T)) : nullptr;
		inds = nz > 0 ? (Int *)MALLOC(nz*sizeof(Int)) : nullptr;
	}
}

template <class T>
void SparseVect<T>::fill(const T val) {
	for (Int i = 0; i < nz; i++)
		v[i] = val;
}

template <class T>
Int SparseVect<T>::size() const {
	return n;
}

template <class T>
Int SparseVect<T>::nnonz() const {
	return nz;
}

template <class T>
void SparseVect<T>::print() const {
	for (Int i = 0; i < nz; i++)
		cout << inds[i] << " " << double135<T> << v[i] << endl;
	cout << endl;
}

template <class T>
void SparseVect<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	
	for (Int i = 0; i < nz; i++)
		f << inds[i] << " " << double2315<T> << v[i] << endl;
	f << endl;

	f.close();
}

template <>
Vector<double> SparseVect<double>::toFullForm() const {
	Vector<double> res(n, 0.0);
	cblas_dsctr(nz, v, inds, res.v);
	return res;
}

template <>
Vector<Complex> SparseVect<Complex>::toFullForm() const {
	Vector<Complex> res(n, zzero);
	cblas_zsctr(nz, v, inds, res.v);
	return res;
}

template <class T>
SparseVect<T>::~SparseVect(void) {
	if (v != nullptr) FREE(v);
	if (inds != nullptr) FREE(inds);
}
