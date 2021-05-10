#include "Vector.h"

template class Vector<double>;
template class Vector<Complex>;

template <class T>
Vector<T>::Vector() : n(0), v(nullptr) { }

template <class T>
Vector<T>::Vector(const Int n) : n(n) {
	assert(n >= 0);
	v = n > 0 ? (T *)MALLOC(n*sizeof(T)) : nullptr;
}

template <class T>
Vector<T>::Vector(const Int n, const T val) : n(n) {
	assert(n >= 0);
	v = n > 0 ? (T *)MALLOC(n*sizeof(T)) : nullptr;
	for (Int i = 0; i < n; i++) v[i] = val;
}

template <class T>
Vector<T>::Vector(const Int n, T *vals) \
	: n(n), v(vals), owner(false) {}

template <class T>
Vector<T>::Vector(const Vector<T> &rhs) : n(rhs.n) {
	assert(n >= 0);
	v = n > 0 ? (T *)MALLOC(n*sizeof(T)) : nullptr;
	for (Int i = 0; i < n; i++) v[i] = rhs.v[i];
}

template <class T>
Vector<T>::Vector(Vector<T> &&rhs) : n(rhs.n) {
	rhs.n = 0;
	v = rhs.v; rhs.v = nullptr;
	owner = rhs.owner; rhs.owner = true;
}

template <class T>
Vector<T> & Vector<T>::operator=(const Vector<T> &rhs) {
	if (this->n != rhs.n) {
		assert(owner);
		//if (v != nullptr && owner) FREE(v);
		if (v != nullptr) FREE(v);
		this->n = rhs.n;
		v = n > 0 ? (T *)MALLOC(n*sizeof(T)) : nullptr;
		//owner = true;
	}
	for (Int i = 0; i < n; i++) this->v[i] = rhs.v[i];
	return *this;
}

template <class T>
Vector<T> & Vector<T>::operator=(Vector<T> &&rhs) {
	if (owner) {
		//if (v != nullptr && owner) FREE(v);
		if (v != nullptr) FREE(v);
		v = rhs.v; rhs.v = nullptr;
		n = rhs.n; rhs.n = 0;
		owner = rhs.owner; rhs.owner = true;
	} else { //not owner, copy
		assert(this->n == rhs.n);
		for (Int i = 0; i < n; i++) this->v[i] = rhs.v[i];
	}
	return *this;
}

template <>
double Vector<double>::scal(const Vector<double> &u) const {
	return cblas_ddot(n, v, 1, u.v, 1);
}

template <>
Complex Vector<Complex>::scal(const Vector<Complex> &u) const {
	Complex res;
	cblas_zdotc_sub(n, v, 1, u.v, 1, &res);
	return res;
}

template <>
double Vector<double>::dotp(const Vector<double> &u) const {
	return scal(u);
}

template <>
Complex Vector<Complex>::dotp(const Vector<Complex> &u) const {
	Complex res;
	cblas_zdotu_sub(n, v, 1, u.v, 1, &res);
	return res;
}

template <>
double Vector<double>::norm2() const {
	return scal(*this);
}

template <>
double Vector<Complex>::norm2() const {
	return scal(*this).real();
}

template <>
double Vector<double>::norm() const {
	return cblas_dnrm2(n, v, 1);
}

template <>
double Vector<Complex>::norm() const {
	return cblas_dznrm2(n, v, 1);
}

template <class T>
void Vector<T>::normalize() {
	double norma = norm();
	assert(norma != 0.0);
	*this *= 1.0 / norma;
}

template <class T>
Vector<T> Vector<T>::operator*(const T x) const {
	Vector<T> res(n);
	for (Int i = 0; i < n; i++) res.v[i] = x*this->v[i];
	return res;
}

template <class T>
Vector<T> Vector<T>::operator+(const Vector<T> &u) const {
	assert(n == u.n);
	Vector<T> res(n);
	for (Int i = 0; i < n; i++) res.v[i] = this->v[i] + u.v[i];
	return res;
}

template <class T>
Vector<T> Vector<T>::operator-(const Vector<T> &u) const {
	assert(n == u.n);
	Vector<T> res(n);
	for (Int i = 0; i < n; i++) res.v[i] = this->v[i] - u.v[i];
	return res;
}

template <>
Vector<double> & Vector<double>::operator*=(const double x) {
	cblas_dscal(n, x, this->v, 1);
	return *this;
}

template <>
Vector<Complex> & Vector<Complex>::operator*=(const Complex x) {
	cblas_zscal(n, &x, this->v, 1);
	return *this;
}

template <class T>
Vector<T> & Vector<T>::operator+=(const Vector<T> &u) {
	assert(n == u.n);
	for (Int i = 0; i < n; i++) this->v[i] += u.v[i];
	return *this;
}

template <class T>
Vector<T> & Vector<T>::operator-=(const Vector<T> &u) {
	assert(n == u.n);
	for (Int i = 0; i < n; i++) this->v[i] -= u.v[i];
	return *this;
}

template <class T>
Int Vector<T>::size() const {
	return n;
}

template <class T>
void Vector<T>::print() const {
	for (Int i = 0; i < n; i++)
		cout << i << " " << double135<T> << v[i] << endl;
	cout << endl;
}

template <class T>
void Vector<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	
	for (Int i = 0; i < n; i++)
		f << i << " " << double2315<T> << v[i] << endl;
	f << endl;

	f.close();
}

template <class T>
void Vector<T>::resize(const Int nnew) {
	assert(nnew >= 0);
	assert(owner);
	if (n != nnew) {
		//if (v != nullptr && owner) FREE(v);
		if (v != nullptr) FREE(v);
		n = nnew;
		v = n > 0 ? (T *)MALLOC(n*sizeof(T)) : nullptr;
		//owner = true;
	}
}

template <class T>
void Vector<T>::setZero() {
	T zero = T();
	for (Int i = 0; i < n; i++) v[i] = zero;
}

template <class T>
void Vector<T>::fill(const T val) {
	for (Int i = 0; i < n; i++)
		v[i] = val;
}

template <class T>
vector<Vector<T>> Vector<T>::makePartition(const vector<Int> &inds) const {
//inds contains numbers of first elements of partition vectors
//inds[0] should be 0
	assert(!inds.empty());
	assert(std::is_sorted(inds.begin(), inds.end()));
	assert( ( inds[0] == 0 ) && ( inds.back() < n ) );
	vector<Vector<T>> res;
	res.reserve(inds.size());
	for (Int k = 0; k < inds.size()-1; k++) {
		res.push_back(Vector<T>(inds[k+1]-inds[k], v + inds[k]));
	}
	res.push_back(Vector<T>(n-inds.back(), v + inds.back()));
	return res;
}

template <class T>
void Vector<T>::conj() {
	for (Int i = 0; i < n; i++)
		v[i] = CONJ(v[i]);
}

template <class T>
Vector<T>::~Vector(void) {
	if (v != nullptr && owner) FREE(v);
}
