#include "GenMatrix.h"

template class GenMatrix<double>;
template class GenMatrix<Complex>;

template <class T>
GenMatrix<T>::GenMatrix() : matr(nullptr) {}

template <class T>
GenMatrix<T>::GenMatrix(const Int n, const Int m) {
	this->n = n;
	this->m = m;
	assert((n >= 0) && (m >= 0));
	matr = n > 0 ? (T **)MALLOC(n*sizeof(T *)) : nullptr;
	if (matr)
		matr[0] = m > 0 ? (T *)MALLOC(n*m*sizeof(T)) : nullptr;
	for (Int i = 1; i < n; i++)
			matr[i] = matr[i-1] + m;
}

template <class T>
GenMatrix<T>::GenMatrix(const Int n) {
	this->n = n;
	this->m = n;
	assert(n >= 0);
	matr = n > 0 ? (T **)MALLOC(n*sizeof(T *)) : nullptr;
	if (matr) matr[0] = (T *)MALLOC(n*m*sizeof(T));
	for (Int i = 1; i < n; i++)
			matr[i] = matr[i-1] + m;
}

template <class T>
GenMatrix<T>::GenMatrix(const GenMatrix<T> &rhs) \
	: AMatrix<T>(rhs), LUFactorized(rhs.LUFactorized), ipiv(rhs.ipiv) {
	matr = (n > 0) ? (T **)MALLOC(n*sizeof(T *)) : nullptr;
	if (matr)
		matr[0] = m > 0 ? (T *)MALLOC(n*m*sizeof(T)) : nullptr;
	for (Int i = 1; i < n; i++)
		matr[i] = matr[i-1] + m;
	for (Int i = 0; i < n; i++)
		for (Int j = 0; j < m; j++)
			matr[i][j] = rhs.matr[i][j];	
}

template <class T>
GenMatrix<T>::GenMatrix(GenMatrix<T> &&rhs) \
	: AMatrix<T>(std::move(rhs)), \
	LUFactorized(rhs.LUFactorized), ipiv(std::move(rhs.ipiv)) {
	rhs.LUFactorized = false;
	matr = rhs.matr;
	rhs.matr = nullptr;
}

template <class T>
GenMatrix<T>::GenMatrix(const GenMatrix<T> &rhs, const Int i1, \
	const Int i2, const Int j1, const Int j2) { //LUFactorized = false
	n = i2 - i1 + 1; m = j2 - j1 + 1;
	assert(i1 >= 0 && j1 >= 0 && i2 < rhs.n && j2 < rhs.m);
	garbage = rhs.garbage;
	if (n == rhs.n && m == rhs.m) LUFactorized = rhs.LUFactorized;
	if (LUFactorized) ipiv = rhs.ipiv;
	matr = (n > 0) ? (T **)MALLOC(n*sizeof(T *)) : nullptr;
	if (matr)
		matr[0] = m > 0 ? (T *)MALLOC(n*m*sizeof(T)) : nullptr;
	for (Int i = 1; i < n; i++)
		matr[i] = matr[i-1] + m;
	for (Int i = 0; i < n; i++)
		for (Int j = 0; j < m; j++)
			matr[i][j] = rhs.matr[i1+i][j1+j];
}

template <class T>
GenMatrix<T> & GenMatrix<T>::operator=(const GenMatrix<T> &rhs) {
	if (this != &rhs) {
		garbage = rhs.garbage;
		LUFactorized = rhs.LUFactorized;
		ipiv = rhs.ipiv;
		if ( (n != rhs.n) || (m != rhs.m) ) {
			if (matr != nullptr) {
				FREE(matr[0]); FREE(matr);
			}
			n = rhs.n; m = rhs.m;
			matr = (n > 0) ? (T **)MALLOC(n*sizeof(T *)) : nullptr;
			if (matr)
				matr[0] = m > 0 ? (T *)MALLOC(n*m*sizeof(T)) : nullptr;
			for (Int i = 1; i < n; i++)
				matr[i] = matr[i-1] + m;	
		}
		for (Int i = 0; i < n; i++)
				for (Int j = 0; j < m; j++)
					matr[i][j] = rhs.matr[i][j];
	}
	return *this;
}

template <class T>
GenMatrix<T> & GenMatrix<T>::operator=(const AMatrix<T> &rhs) {
	const GenMatrix<T> & rhs_ = dynamic_cast<const GenMatrix<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
GenMatrix<T> & GenMatrix<T>::operator=(GenMatrix<T> &&rhs) {
	garbage = rhs.garbage; rhs.garbage = false;
	LUFactorized = rhs.LUFactorized; rhs.LUFactorized = false;
	ipiv = std::move(rhs.ipiv);
	n = rhs.n; m = rhs.m; rhs.n = 0; rhs.m = 0;
	if (matr != nullptr) {
		if (matr[0] != nullptr) FREE(matr[0]);
		FREE(matr);
	}
	matr = rhs.matr; rhs.matr = nullptr;
	return *this;
}

template <class T>
GenMatrix<T> & GenMatrix<T>::operator=(AMatrix<T> &&rhs) {
	GenMatrix<T> && rhs_ = dynamic_cast<GenMatrix<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void GenMatrix<T>::resize(const Int nnew, const Int mnew) {
	assert((nnew >= 0) && (mnew >= 0));
	garbage = false; LUFactorized = false;
	ipiv.clear(); ipiv.shrink_to_fit();
	if ( (n != nnew) || (m != mnew) ) {
		if (matr != nullptr) {
			FREE(matr[0]); FREE(matr);
		}
		n = nnew; m = mnew;
		matr = n > 0 ? (T **)MALLOC(n*sizeof(T *)) : nullptr;
		if (matr) matr[0] = m > 0 ? (T *)MALLOC(n*m*sizeof(T)) : nullptr;
		for (Int i = 1; i < n; i++)
			matr[i] = matr[i-1] + m;
	}
}

template <class T>
void GenMatrix<T>::resize(const Int nnew) {
	resize(nnew, nnew);
}

template <class T>
void GenMatrix<T>::fill(const T val) {
	garbage = false; LUFactorized = false;
	ipiv.clear(); ipiv.shrink_to_fit();
	for (Int i = 0; i < n; i++)
		for (Int j = 0; j < m; j++)
			matr[i][j] = val;
}

template <>
void GenMatrix<double>::solve(Vector<double> & rhssol) {
	assert(n == m);
	Int size = rhssol.size();
	assert(size <= n);
	assert(!garbage || (LUFactorized && (size == n)));
	Int info;
	if (!LUFactorized) {
		ipiv.resize(size);
		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, size, size, matr[0], n, ipiv.data());
		assert(info == 0);
		garbage = true;
		if (size == n) LUFactorized = true;
	}

	info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', size, \
		1, matr[0], n, ipiv.data(), rhssol.v, 1);
	assert(info == 0);
}

template <>
void GenMatrix<Complex>::solve(Vector<Complex> & rhssol) {
	assert(n == m);
	Int size = rhssol.size();
	assert(size <= n);
	assert(!garbage || (LUFactorized && (size == n)));
	Int info;
	if (!LUFactorized) {
		ipiv.resize(size);
		info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, size, size, matr[0], n, ipiv.data());
		assert(info == 0);
		garbage = true;
		if (size == n) LUFactorized = true;
	}

	info = LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N', size, \
		1, matr[0], n, ipiv.data(), rhssol.v, 1);
	assert(info == 0);
}

template <class T>
void GenMatrix<T>::print() const {
	for (Int i = 0; i < n; i++) {
		for (Int j = 0; j < m; j++) {
			cout << double135<T> << matr[i][j] << "  ";
		}
		cout << endl;
	}
}

template <class T>
void GenMatrix<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	for (Int i = 0; i < n; i++) {
		for (Int j = 0; j < m; j++) {
			f << double2315<T> << matr[i][j] << "  ";
		}
		f << endl;
	}
	f.close();
}

template <>
GenMatrix<double> & GenMatrix<double>::operator*=(const double c) {
	assert(!garbage);
	cblas_dscal(n*m, c, matr[0], 1);
	return *this;
}

template <>
GenMatrix<Complex> & GenMatrix<Complex>::operator*=(const Complex c) {
	assert(!garbage);
	cblas_zscal(n*m, &c, matr[0], 1);
	return *this;
}

/*template <class T>
GenMatrix<T> & GenMatrix<T>::operator+=(const AMatrix<T> &b) {
	assert(!garbage);
	const GenMatrix<T> & b_ = \
		dynamic_cast<const GenMatrix<T> &>(b);
	assert(n == b_.n || m == b_.m);
	for (Int i = 0; i < n; i++) {
		for (Int j = 0; j < m; j++) {
			matr[i][j] += b_.matr[i][j];
		}
	}
	return *this;
}*/

template <>
Vector<double> GenMatrix<double>::operator*(const Vector<double> &vec) const {
	assert(!garbage);
	Int size = vec.size();
	assert(size <= m);
	Vector<double> y(n, 0.0);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, n, size, 1.0, \
		matr[0], m, vec.v, 1, 0.0, y.v, 1);
	return y;
}

template <>
Vector<Complex> GenMatrix<Complex>::operator*(const Vector<Complex> &vec) const {
	assert(!garbage);
	Int size = vec.size();
	assert(size <= m);
	Vector<Complex> y(n, zzero);
	cblas_zgemv(CblasRowMajor, CblasNoTrans, n, size, &zone, \
		matr[0], m, vec.v, 1, &zzero, y.v, 1);
	return y;
}

template <class T>
double GenMatrix<T>::sizeGb() const {
	return n * m * sizeof(T)*bytes2Gbytes;
}

template <>
void GenMatrix<double>::getEEV(Vector<Complex> &ev, GenMatrix<Complex> &wl, \
				GenMatrix<Complex> &wr, const bool computeLeft, \
				const bool computeRight) {
//NB! returns right eigenvectors in COLUMNS of wr
// left eigenvectors in COLUMNS of wl
	assert(!garbage);
	assert(n == m);
	ev.resize(0); wl.resize(0); wr.resize(0);
	Vector<double> evRe(n), evIm(n);
	GenMatrix<double> wl_(0), wr_(0);
	Int ldwl = computeLeft ? n : 1;
	wl_.resize(n, ldwl);
	Int ldwr = computeRight ? n : 1;
	wr_.resize(n, ldwr);
	ldwl = n; ldwr = n;
	Int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, computeLeft ? 'V' : 'N', \
				computeRight ? 'V' : 'N', n, matr[0], n, evRe.v, evIm.v, \
				wl_.matr[0], ldwl, wr_.matr[0], ldwr);
	assert(info == 0);

	ev.resize(n);
	for (Int i = 0; i < n; i++)
		ev[i] = Complex(evRe[i], evIm[i]);
	evRe.resize(0); evIm.resize(0);

	if (computeLeft) {
		wl.resize(n);
		for (Int j = 0; j < n; j++) {
			if ( IS_EPS(ev[j].imag() / ev[j].real()) ) {
				for (Int i = 0; i < n; i++)
					wl[i][j] = wl_[i][j];
			} else { //complex pair
				for (Int i = 0; i < n; i++) {
					wl[i][j] = Complex(wl_[i][j], wl_[i][j+1]);
					wl[i][j+1] = Complex(wl_[i][j], -wl_[i][j+1]);
				}
				j++;
			}
		}
		wl_.resize(0);
	}
	if (computeRight) {
		wr.resize(n);
		for (Int j = 0; j < n; j++) {
			if ( IS_EPS(ev[j].imag() / ev[j].real()) ) {
				for (Int i = 0; i < n; i++)
					wr[i][j] = wr_[i][j];
			} else { //complex pair
				for (Int i = 0; i < n; i++) {
					wr[i][j] = Complex(wr_[i][j], wr_[i][j+1]);
					wr[i][j+1] = Complex(wr_[i][j], -wr_[i][j+1]);
				}
				j++;
			}
		}
		wr_.resize(0);
	}
	garbage = true;
}

template <>
void GenMatrix<Complex>::getEEV(Vector<Complex> &ev, GenMatrix<Complex> &wl, \
				GenMatrix<Complex> &wr, const bool computeLeft, \
				const bool computeRight) {
//NB! returns right eigenvectors in COLUMNS of wr
// left eigenvectors in COLUMNS of wl
	assert(!garbage);
	assert(n == m);
	ev.resize(n);
	Int ldwl = computeLeft ? n : 1;
	wl.resize(n, ldwl);
	Int ldwr = computeRight ? n : 1;
	wr.resize(n, ldwr);
	ldwl = n; ldwr = n;
	Int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, computeLeft ? 'V' : 'N', \
				computeRight ? 'V' : 'N', n, matr[0], n, ev.v, \
				wl.matr[0], ldwl, wr.matr[0], ldwr);
	assert(info == 0);
	if (!computeLeft) wl.resize(0);
	if (!computeRight) wr.resize(0);
	garbage = true;
}

template<class T>
void GenMatrix<T>::cconj() {
	trans();
	for (Int i = 0; i < n; i++)
		for (Int j = 0; j < n; j++)
			matr[i][j] = CONJ(matr[i][j]);
}

template<class T>
void GenMatrix<T>::trans() {
	assert(!garbage);
	assert(n == m); //not realized for non square matrices
	for (Int i = 0; i < n; i++)
		for (Int j = i+1; j < n; j++)
			SWAP(matr[i][j], matr[j][i]);
}

template <>
void GenMatrix<double>::inv() {
	assert(!garbage || LUFactorized);
	assert(n == m);
	//Int *ipiv = (Int *)MALLOC( n * sizeof(Int) );
	Int info;
	if (!LUFactorized) {
		ipiv.resize(n);
		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, matr[0], n, ipiv.data());
		assert(info == 0);
	}
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, matr[0], n, ipiv.data());
	assert(info == 0);
	//FREE(ipiv);
	garbage = false; LUFactorized = false;
}

template <>
void GenMatrix<Complex>::inv() {
	assert(!garbage || LUFactorized);
	assert(n == m);
	//Int *ipiv = (Int *)MALLOC( n * sizeof(Int) );
	Int info;
	if (!LUFactorized) {
		ipiv.resize(n);
		info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, matr[0], n, ipiv.data());
		assert(info == 0);
	}
	info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, matr[0], n, ipiv.data());
	assert(info == 0);
	//FREE(ipiv);
	garbage = false; LUFactorized = false;
}

template <class T>
GenMatrix<T>::~GenMatrix(void) {
	if (matr != nullptr) {
		if (matr[0] != nullptr) FREE(matr[0]);
		FREE(matr);
	}
}
