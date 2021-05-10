#include "TriDiagM.h"

template class TriDiagM<double>;
template class TriDiagM<Complex>;

template <class T>
TriDiagM<T>::TriDiagM(void) { }

template <class T>
TriDiagM<T>::TriDiagM(const Int n) : d(n), u(n-1), l(n-1) {
	this->n = n;
	this->m = n;
}

template <class T>
TriDiagM<T>::TriDiagM(TriDiagM<T> &&rhs) \
	: AMatrix<T>(std::move(rhs)), LUFactorized(rhs.LUFactorized), \
	u2(std::move(rhs.u2)), ipiv(std::move(rhs.ipiv)), \
	d(std::move(rhs.d)), u(std::move(rhs.u)), l(std::move(rhs.l)) {
	rhs.LUFactorized = false;
}

template <class T>
TriDiagM<T> & TriDiagM<T>::operator=(const TriDiagM<T> &rhs) {
	if (this != &rhs) {
		n = rhs.n; m = rhs.m;
		garbage = rhs.garbage;
		d = rhs.d; u = rhs.u; l = rhs.l;
		ipiv = rhs.ipiv; u2 = rhs.u2;
		LUFactorized = rhs.LUFactorized;
	}
	return *this;
}

template <class T>
TriDiagM<T> & TriDiagM<T>::operator=(const AMatrix<T> &rhs) {
	const TriDiagM<T> & rhs_ = dynamic_cast<const TriDiagM<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
TriDiagM<T> & TriDiagM<T>::operator=(TriDiagM<T> &&rhs) {
	n = rhs.n; rhs.n = 0;
	m = rhs.m; rhs.m = 0;
	garbage = rhs.garbage; rhs.garbage = false;
	d = std::move(rhs.d); u = std::move(rhs.u); l = std::move(rhs.l);
	ipiv = std::move(ipiv); u2 = std::move(u2);
	LUFactorized = rhs.LUFactorized; rhs.LUFactorized = false;
	return *this;
}

template <class T>
TriDiagM<T> & TriDiagM<T>::operator=(AMatrix<T> &&rhs) {
	TriDiagM<T> &&rhs_ = dynamic_cast<TriDiagM<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void TriDiagM<T>::fill(const T val) {
	garbage = false; LUFactorized = false;
	ipiv.clear(); ipiv.shrink_to_fit(); u2.resize(0);
	d.fill(val); u.fill(val); l.fill(val);
}

template <class T>
void TriDiagM<T>::resize(const Int nnew) {
	assert(nnew >= 0);
	n = nnew;
	garbage = false; LUFactorized = false;
	ipiv.clear(); ipiv.shrink_to_fit(); u2.resize(0);
	d.resize(nnew); u.resize(nnew-1); l.resize(nnew-1);
}

template <>
void TriDiagM<double>::solve(Vector<double> & rhssol) {
	Int size = rhssol.size();
	assert(size <= n);
	assert(!garbage || (LUFactorized && (size == n)));
	
	Int info;
	if (!LUFactorized) {
		ipiv.resize(size); u2.resize(max((Int)0, size-2));
		info = LAPACKE_dgttrf(size, l.v, d.v, u.v, u2.v, ipiv.data());
		assert(info == 0);
		garbage = true;
		if (size == n) LUFactorized = true;
	}

	info = LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', size, \
		1, l.v, d.v, u.v, u2.v, ipiv.data(), rhssol.v, 1);
	assert(info == 0);
}

template <>
void TriDiagM<Complex>::solve(Vector<Complex> & rhssol) {
	Int size = rhssol.size();
	assert(size <= n);
	assert(!garbage || (LUFactorized && (size == n)));
	
	Int info;
	if (!LUFactorized) {
		ipiv.resize(size); u2.resize(max((Int)0, size-2));
		info = LAPACKE_zgttrf(size, l.v, d.v, u.v, u2.v, ipiv.data());
		assert(info == 0);
		garbage = true;
		if (size == n) LUFactorized = true;
	}

	info = LAPACKE_zgttrs(LAPACK_ROW_MAJOR, 'N', size, \
		1, l.v, d.v, u.v, u2.v, ipiv.data(), rhssol.v, 1);
	assert(info == 0);
}

template <class T>
void TriDiagM<T>::print() const {
	cout << "Superdiagonal, diagonal, subdiagonal elements:" << endl;
	u.print(); d.print(); l.print();
}

template <class T>
void TriDiagM<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	f << "Superdiagonal, diagonal, subdiagonal elements:" << endl;
	f.close();
	u.write(filename); d.write(filename); l.write(filename);
}

template <class T>
TriDiagM<T> & TriDiagM<T>::operator*=(const T c) {
	assert(!garbage);
	d *= c; u *= c; l *= c;
	return *this;
}

template <class T>
Vector<T> TriDiagM<T>::operator*(const Vector<T> &vec) const {
	assert(!garbage);
	Int size = vec.size();
	assert(size <= n);
	Vector<T> y(vec.size());
	if (size >= 2) {
		y[0] = d[0]*vec[0]+u[0]*vec[1];
		y[size-1] = l[size-2]*vec[size-2]+d[size-1]*vec[size-1];
	} else //size = 1
		y[0] = d[0]*vec[0];
	for (Int i = 1; i < size-1; i++)
		y[i] = l[i-1]*vec[i-1]+d[i]*vec[i]+u[i]*vec[i+1];

	//TODO assert size != 0

	return y;
}

template <class T>
double TriDiagM<T>::sizeGb() const {
	return (3*n - 2) * sizeof(T)*bytes2Gbytes;
}