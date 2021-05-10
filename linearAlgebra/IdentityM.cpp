#include "IdentityM.h"

template class IdentityM<double>;
template class IdentityM<Complex>;

template <class T>
IdentityM<T>::IdentityM(const Int n) {
	this->n = n;
	this->m = n;
}

template <class T>
IdentityM<T> & IdentityM<T>::operator=(const IdentityM<T> &rhs) {
	if (this != &rhs) {
		n = rhs.n;
		m = rhs.m;
	}
	return *this;
}

template <class T>
IdentityM<T> & IdentityM<T>::operator=(const AMatrix<T> &rhs) {
	const IdentityM<T> & rhs_ = dynamic_cast<const IdentityM<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
IdentityM<T> & IdentityM<T>::operator=(IdentityM<T> &&rhs) {
	n = rhs.n; rhs.n = 0;
	m = rhs.m; rhs.m = 0;
	return *this;
}

template <class T>
IdentityM<T> & IdentityM<T>::operator=(AMatrix<T> &&rhs) {
	IdentityM<T> &&rhs_ = dynamic_cast<IdentityM<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void IdentityM<T>::fill(const T val) {
	assert(false);
}

template <class T>
void IdentityM<T>::resize(const Int nnew) {
	assert(nnew >= 0);
	n = nnew;
}

template <class T>
void IdentityM<T>::solve(Vector<T> & rhssol) {
	Int size = rhssol.size();
	assert(size <= n);
}

template <class T>
void IdentityM<T>::print() const {
	cout << "Identity " << n << " by " << n << " matrix" << endl;
}

template <class T>
void IdentityM<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	f << "Identity " << n << " by " << n << " matrix" << endl;
	f.close();
}

template <class T>
IdentityM<T> & IdentityM<T>::operator*=(const T c) {
	assert(false);
	return *this;
}

template <class T>
Vector<T> IdentityM<T>::operator*(const Vector<T> &vec) const {
	Int size = vec.size();
	assert(size <= n);
	Vector<T> y(vec);
	return y;
}

template <class T>
double IdentityM<T>::sizeGb() const {
	return 0.0;
}