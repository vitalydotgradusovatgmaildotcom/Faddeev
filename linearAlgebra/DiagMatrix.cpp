#include "DiagMatrix.h"

template class DiagMatrix<double>;
template class DiagMatrix<Complex>;

template <class T>
DiagMatrix<T>::DiagMatrix(void) { }

template <class T>
DiagMatrix<T>::DiagMatrix(const Int n) : d(n) {
	this->n = n;
	this->m = n;
}

template <class T>
DiagMatrix<T> & DiagMatrix<T>::operator=(const DiagMatrix<T> &rhs) {
	if (this != &rhs) {
		n = rhs.n; m = rhs.m;
		garbage = rhs.garbage;
		d = rhs.d;
	}
	return *this;
}

template <class T>
DiagMatrix<T> & DiagMatrix<T>::operator=(const AMatrix<T> &rhs) {
	const DiagMatrix<T> & rhs_ = dynamic_cast<const DiagMatrix<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
DiagMatrix<T> & DiagMatrix<T>::operator=(DiagMatrix<T> &&rhs) {
	n = rhs.n; rhs.n = 0;
	m = rhs.m; rhs.m = 0;
	garbage = rhs.garbage; rhs.garbage = false;
	d = std::move(rhs.d);
	return *this;
}

template <class T>
DiagMatrix<T> & DiagMatrix<T>::operator=(AMatrix<T> &&rhs) {
	DiagMatrix<T> &&rhs_ = dynamic_cast<DiagMatrix<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void DiagMatrix<T>::fill(const T val) {
	garbage = false;
	d.fill(val);
}

template <class T>
void DiagMatrix<T>::resize(const Int nnew) {
	assert(nnew >= 0);
	n = nnew;
	garbage = false;
	d.resize(nnew);
}

template <class T>
void DiagMatrix<T>::solve(Vector<T> & rhssol) {
	assert(!garbage);
	Int size = rhssol.size();
	assert(size <= n);
	for (Int i = 0; i < size; i++)
		rhssol[i] /= d[i];
}

template <class T>
void DiagMatrix<T>::print() const {
	d.print();
}

template <class T>
void DiagMatrix<T>::write(const string &filename) const {
	d.write(filename);
}

template <class T>
DiagMatrix<T> & DiagMatrix<T>::operator*=(const T c) {
	assert(!garbage);
	d *= c;
	return *this;
}

/*template <class T>
DiagMatrix<T> & DiagMatrix<T>::operator+=(const AMatrix<T> &b) {
	assert(!garbage);
	const DiagMatrix<T> & b_ = \
		dynamic_cast<const DiagMatrix<T> &>(b);
	assert(n == b_.n);
	d += b_.d;
	return *this;
}*/

template <class T>
Vector<T> DiagMatrix<T>::operator*(const Vector<T> &vec) const {
	assert(!garbage);
	Int size = vec.size();
	assert(size <= n);
	Vector<T> y(vec);
	for (Int i = 0; i < size; i++)
		y[i] *= d[i];
	return y;
}

template <class T>
double DiagMatrix<T>::sizeGb() const {
	return n * sizeof(T)*bytes2Gbytes;
}
