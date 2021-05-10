#include "UTMatrix.h"

template class UTMatrix<double>;
template class UTMatrix<Complex>;

template <class T>
UTMatrix<T>::UTMatrix(void) : matr(nullptr) { }

template <class T>
UTMatrix<T>::UTMatrix(const Int n) {
	this->n = n;
	this->m = n;
	assert(n >= 0);
	nnonz = n*(n+1)/2;
	matr = n > 0 ? (T **)MALLOC(n*sizeof(T *)) : nullptr;
	if (matr) matr[0] = (T *)MALLOC(nnonz*sizeof(T));
	for (Int i = 1; i < n; i++)
		matr[i] = matr[i-1] + (n-i+1);
}

template <class T>
UTMatrix<T>::UTMatrix(const UTMatrix<T> &rhs) : AMatrix<T>(rhs) {
	nnonz = n*(n+1)/2;
	matr = (n > 0) ? (T **)MALLOC(n*sizeof(T *)) : nullptr;
	if (matr)
		matr[0] = (T *)MALLOC(nnonz*sizeof(T));
	for (Int i = 1; i < n; i++)
		matr[i] = matr[i-1] + (n-i+1);
	for (Int i = 0; i < n; i++)
		for (Int j = 0; j < n-i; j++)
			matr[i][j] = rhs.matr[i][j];
}


template <class T>
UTMatrix<T>::UTMatrix(UTMatrix<T> &&rhs) \
		: AMatrix<T>(std::move(rhs)), nnonz(rhs.nnonz) {
	rhs.nnonz = 0;
	matr = rhs.matr; rhs.matr = nullptr;
}

template <class T>
UTMatrix<T> & UTMatrix<T>::operator=(const UTMatrix<T> &rhs) {
	if (this != &rhs) {
		garbage = rhs.garbage;
		if (n != rhs.n) {
			if (matr != nullptr) {
				FREE(matr[0]); FREE(matr);
			}
			n = rhs.n; m = rhs.m;
			nnonz = n*(n+1)/2;
			matr = (n > 0) ? (T **)MALLOC(n*sizeof(T *)) : nullptr;
			if (matr)
				matr[0] = (T *)MALLOC(nnonz*sizeof(T));
			for (Int i = 1; i < n; i++)
				matr[i] = matr[i-1] + (n-i+1);
		}
		for (Int i = 0; i < n; i++)
			for (Int j = 0; j < n-i; j++)
				matr[i][j] = rhs.matr[i][j];
	}
	return *this;
}

template <class T>
UTMatrix<T> & UTMatrix<T>::operator=(const AMatrix<T> &rhs) {
	const UTMatrix<T> & rhs_ = dynamic_cast<const UTMatrix<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
UTMatrix<T> & UTMatrix<T>::operator=(UTMatrix<T> &&rhs) {
	garbage = rhs.garbage; rhs.garbage = false;
	n = rhs.n; rhs.n = 0;
	m = rhs.m; rhs.m = 0;
	nnonz = rhs.nnonz; rhs.nnonz = 0;
	if (matr != nullptr) {
		FREE(matr[0]); FREE(matr);
	}
	matr = rhs.matr; rhs.matr = nullptr;
	return *this;
}

template <class T>
UTMatrix<T> & UTMatrix<T>::operator=(AMatrix<T> &&rhs) {
	UTMatrix<T> && rhs_ = dynamic_cast<UTMatrix<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void UTMatrix<T>::resize(const Int nnew) {
	assert(nnew >= 0);
	garbage = false;
	if (n != nnew) {
		if (matr != nullptr) {
			FREE(matr[0]); FREE(matr);
		}
		n = nnew;
		nnonz = n*(n+1)/2;
		matr = n > 0 ? (T **)MALLOC(n*sizeof(T *)) : nullptr;
		if (matr) matr[0] = (T *)MALLOC(nnonz*sizeof(T));
		for (Int i = 1; i < n; i++)
				matr[i] = matr[i-1] + (n-i+1);
	}
}

template <class T>
void UTMatrix<T>::fill(const T val) {
	garbage = false;
	for (Int i = 0; i < n; i++)
		for (Int j = 0; j < n-i; j++)
			matr[i][j] = val;
}

template <>
void UTMatrix<double>::solve(Vector<double> & rhssol) /*const*/ {
	assert(!garbage);
	Int size = rhssol.size();
	assert(size <= n);
	UTMatrix<double> *m_tmp = this;
	if (size != n) {
		m_tmp = new UTMatrix<double>(size);
		for (Int i = 0; i < m_tmp->n; i++)
			for (Int j = 0; j < m_tmp->n-i; j++)
				m_tmp->matr[i][j] = matr[i][j];
	}
	cblas_dtpsv(CblasRowMajor, CblasUpper, CblasNoTrans, \
		CblasNonUnit, size, m_tmp->matr[0], rhssol.v, 1);
	if (size != n) delete m_tmp;
}

template <>
void UTMatrix<Complex>::solve(Vector<Complex> & rhssol)  /*const*/ {
	assert(!garbage);
	Int size = rhssol.size();
	assert(size <= n);
	UTMatrix<Complex> *m_tmp = this;
	if (size != n) {
		m_tmp = new UTMatrix<Complex>(size);
		for (Int i = 0; i < m_tmp->n; i++)
			for (Int j = 0; j < m_tmp->n-i; j++)
				m_tmp->matr[i][j] = matr[i][j];
	}
	cblas_ztpsv(CblasRowMajor, CblasUpper, CblasNoTrans, \
		CblasNonUnit, size, m_tmp->matr[0], rhssol.v, 1);
	if (size != n) delete m_tmp;
}

template <>
void UTMatrix<double>::print() const {
	for (Int i = 0; i < n; i++) {
		for (Int j = 0; j < i; j++)
			cout << setw(6) << ' ';
		for (Int j = 0; j < n-i; j++) {
			cout << double135<double> << matr[i][j] << "  ";
		}
		cout << endl;
	}
}

template <>
void UTMatrix<Complex>::print() const {
	for (Int i = 0; i < n; i++) {
		for (Int j = 0; j < i; j++)
			cout << setw(13) << ' ';
		for (Int j = 0; j < n-i; j++) {
			cout << double135<Complex> << matr[i][j] << "  ";
		}
		cout << endl;
	}
}

template <>
void UTMatrix<double>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	for (Int i = 0; i < n; i++) {
		for (Int j = 0; j < i; j++)
			f << setw(14) << ' ';
		for (Int j = 0; j < n-i; j++) {
			f << double2315<double> << matr[i][j] << "  ";
		}
		f << endl;
	}
	f.close();
}

template <>
void UTMatrix<Complex>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	for (Int i = 0; i < n; i++) {
		for (Int j = 0; j < i; j++)
			f << setw(29) << ' ';
		for (Int j = 0; j < n-i; j++) {
			f << double2315<Complex> << matr[i][j] << "  ";
		}
		f << endl;
	}
	f.close();
}

template <>
UTMatrix<double> & UTMatrix<double>::operator*=(const double c) {
	assert(!garbage);
	cblas_dscal(nnonz, c, matr[0], 1);
	return *this;
}

template <>
UTMatrix<Complex> & UTMatrix<Complex>::operator*=(const Complex c) {
	assert(!garbage);
	cblas_zscal(nnonz, &c, matr[0], 1);
	return *this;
}

/*template <class T>
UTMatrix<T> & UTMatrix<T>::operator+=(const AMatrix<T> &b) {
	assert(!garbage);
	const UTMatrix<T> & b_ = \
		dynamic_cast<const UTMatrix<T> &>(b);
	assert(n == b_.n);
	for (Int i = 0; i < n; i++) {
		for (Int j = 0; j < n-i; j++) {
			matr[i][j] += b_.matr[i][j];
		}
	}
	return *this;
}*/

template <>
Vector<double> UTMatrix<double>::operator*(const Vector<double> &vec) const {
	assert(!garbage);
	Int size = vec.size();
	assert(size <= n);
	const UTMatrix<double> *m_tmp;
	if (size != n) {
		m_tmp = new UTMatrix<double>(size);
		for (Int i = 0; i < m_tmp->n; i++)
			for (Int j = 0; j < m_tmp->n-i; j++)
				m_tmp->matr[i][j] = matr[i][j];
	} else
		 m_tmp = this;
	Vector<double> y(vec);
	cblas_dtpmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, \
		size, m_tmp->matr[0], y.v, 1);
	if (size != n) delete m_tmp;
	return y;
}

template <>
Vector<Complex> UTMatrix<Complex>::operator*(const Vector<Complex> &vec) const {
	assert(!garbage);
	Int size = vec.size();
	assert(size <= n);
	const UTMatrix<Complex> *m_tmp;
	if (size != n) {
		m_tmp = new UTMatrix<Complex>(size);
		for (Int i = 0; i < m_tmp->n; i++)
			for (Int j = 0; j < m_tmp->n-i; j++)
				m_tmp->matr[i][j] = matr[i][j];
	} else
		 m_tmp = this;
	Vector<Complex> y(vec);
	cblas_ztpmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, \
		size, m_tmp->matr[0], y.v, 1);
	if (size != n) delete m_tmp;
	return y;
}

template <class T>
double UTMatrix<T>::sizeGb() const {
	return nnonz * sizeof(T)*bytes2Gbytes;
}

template <class T>
UTMatrix<T>::~UTMatrix(void) {
	if (matr != nullptr) {
		FREE(matr[0]); FREE(matr);
	}
}
