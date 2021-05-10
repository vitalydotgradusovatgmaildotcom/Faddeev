#include "MatrixProd.h"

template class MatrixProd<double>;
template class MatrixProd<Complex>;

template <class T>
MatrixProd<T>::MatrixProd(const std::vector<std::shared_ptr<AMatrix<T>>> &matrs) \
			: matrs(matrs) {
	assert(!matrs.empty());
	assert( std::find(matrs.begin(), matrs.end(), nullptr) == matrs.end() );
	n = matrs[0]->nrows();
	m = matrs.back()->ncols();
	for (Int k = 1; k < matrs.size(); k++)
		assert( matrs[k]->nrows() == matrs[k-1]->ncols() );
	inverse.resize(matrs.size());
	std::fill( inverse.begin(), inverse.end(), false );
}

template <class T>
MatrixProd<T>::MatrixProd(const std::vector<std::shared_ptr<AMatrix<T>>> &matrs, \
			const std::vector<bool> &inverse) \
			: matrs(matrs), inverse(inverse) {
	assert(!matrs.empty());
	assert( std::find(matrs.begin(), matrs.end(), nullptr) == matrs.end() );
	n = matrs[0]->nrows();
	m = matrs.back()->ncols();
	for (Int k = 1; k < matrs.size(); k++)
		assert( matrs[k]->nrows() == matrs[k-1]->ncols() );
	assert(matrs.size() == inverse.size());
	for (Int k = 0; k < matrs.size(); k++)
		assert( !(inverse[k] && \
			( matrs[k]->nrows() !=  matrs[k]->ncols())) );
}

template <class T>
MatrixProd<T> & MatrixProd<T>::operator=(const MatrixProd<T> &rhs) {
	if (this != &rhs) {
		n = rhs.n; m = rhs.m;
		matrs = rhs.matrs;
		inverse = rhs.inverse;
	}
	return *this;
}

template <class T>
MatrixProd<T> & MatrixProd<T>::operator=(const AMatrix<T> &rhs) {
	const MatrixProd<T> &rhs_ = dynamic_cast<const MatrixProd<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
MatrixProd<T> & MatrixProd<T>::operator=(MatrixProd<T> &&rhs) {
	n = rhs.n; rhs.n = 0;
	m = rhs.m; rhs.m = 0;
	matrs = std::move(rhs.matrs);
	inverse = std::move(rhs.inverse);
	return *this;
}

template <class T>
MatrixProd<T> & MatrixProd<T>::operator=(AMatrix<T> &&rhs) {
	MatrixProd<T> && rhs_ = dynamic_cast<MatrixProd<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void MatrixProd<T>::fill(const T val) {
	assert(false);
}

template <class T>
void MatrixProd<T>::resize(const Int nnew) {
	assert(false);
}

template <class T>
void MatrixProd<T>::solve(Vector<T> &rhssol) {
	assert(n == m);
	Int size = rhssol.size();
	assert(size == n);
	for (Int k = 0; k < matrs.size(); k++) {
		if (inverse[k])
			rhssol = *(matrs[k])*rhssol;
		else
			matrs[k]->solve(rhssol);
	}
}

template <class T>
void MatrixProd<T>::print() const {
	cout << "Product of matrices:" << endl;
	for (Int k = 0; k < matrs.size(); k++) {
		if (inverse[k])
			cout << "Inverse of matrix " << k << ":" << endl;
		else
			cout << "Matrix " << k << ":" << endl;
		matrs[k]->print();
	}
}

template <class T>
void MatrixProd<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	f << "Product of matrices:" << endl;
	for (Int k = 0; k < matrs.size(); k++) {
		if (inverse[k])
			f << "Inverse of matrix " << k << ":" << endl;
		else
			f << "Matrix " << k << ":" << endl;
		matrs[k]->write(filename);
	}
	f.close();
}

template <class T>
MatrixProd<T> & MatrixProd<T>::operator*=(const T c) {
	std::vector<bool>::iterator it = \
		std::find(inverse.begin(), inverse.end(), false);
	if (it != inverse.end())
		*(matrs[it-inverse.begin()]) *= c;
	else {
		assert(c != T());
		*(matrs[0]) *= 1.0/c;
	}
	return *this;
}

template <class T>
Vector<T> MatrixProd<T>::operator*(const Vector<T> &vec) const {
	Int size = vec.size();
	assert(size == m);
	Vector<T> res(vec);
	for (Int k = matrs.size()-1; k >= 0 ; k--) {
		if (inverse[k])
			matrs[k]->solve(res);
		else
			res = *matrs[k]*res;
	}
	return res;
}

template <class T>
double MatrixProd<T>::sizeGb() const {
	assert(false);
}