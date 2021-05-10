#include "MatrixSum.h"

template class MatrixSum<double>;
template class MatrixSum<Complex>;

template <class T>
MatrixSum<T>::MatrixSum(const std::vector<std::shared_ptr<AMatrix<T>>> &matrs) \
	: matrs(matrs) {
	assert(!matrs.empty());
	assert(std::find(matrs.begin(), matrs.end(), nullptr) == matrs.end());
	n = matrs[0]->nrows();
	m = matrs[0]->ncols();
	for (Int k = 1; k < matrs.size(); k++) {
		assert(matrs[k]->nrows() == n);
		assert(matrs[k]->ncols() == m);
	}
	inverse.resize(matrs.size());
	std::fill(inverse.begin(), inverse.end(), false);
}

template <class T>
MatrixSum<T>::MatrixSum(const std::vector<std::shared_ptr<AMatrix<T>>> &matrs, \
	const std::vector<bool> &inverse) \
	: matrs(matrs), inverse(inverse) {
	assert(!matrs.empty());
	assert(std::find(matrs.begin(), matrs.end(), nullptr) == matrs.end());
	n = matrs[0]->nrows();
	m = matrs[0]->ncols();
	for (Int k = 1; k < matrs.size(); k++) {
		assert(matrs[k]->nrows() == n);
		assert(matrs[k]->ncols() == m);
	}
	assert(matrs.size() == inverse.size());
	for (Int k = 0; k < matrs.size(); k++)
		assert( !(inverse[k] && \
			( matrs[k]->nrows() != matrs[k]->ncols() )) );
}

template <class T>
MatrixSum<T> & MatrixSum<T>::operator=(const MatrixSum<T> &rhs) {
	if (this != &rhs) {
		n = rhs.n; m = rhs.m;
		matrs = rhs.matrs;
		inverse = rhs.inverse;
	}
	return *this;
}

template <class T>
MatrixSum<T> & MatrixSum<T>::operator=(const AMatrix<T> &rhs) {
	const MatrixSum<T> &rhs_ = dynamic_cast<const MatrixSum<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
MatrixSum<T> & MatrixSum<T>::operator=(MatrixSum<T> &&rhs) {
	n = rhs.n; rhs.n = 0;
	m = rhs.m; rhs.m = 0;
	matrs = std::move(rhs.matrs);
	inverse = std::move(rhs.inverse);
	return *this;
}

template <class T>
MatrixSum<T> & MatrixSum<T>::operator=(AMatrix<T> &&rhs) {
	MatrixSum<T> && rhs_ = dynamic_cast<MatrixSum<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void MatrixSum<T>::fill(const T val) {
	assert(false);
}

template <class T>
void MatrixSum<T>::resize(const Int nnew) {
	assert(false);
}

template <class T>
void MatrixSum<T>::solve(Vector<T> &rhssol) {
	assert(false);
}

template <class T>
void MatrixSum<T>::print() const {
	cout << "Sum of matrices:" << endl;
	for (Int k = 0; k < matrs.size(); k++) {
		if (inverse[k])
			cout << "Inverse of matrix " << k << ":" << endl;
		else
			cout << "Matrix " << k << ":" << endl;
		matrs[k]->print();
	}
}

template <class T>
void MatrixSum<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	f << "Sum of matrices:" << endl;
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
MatrixSum<T> & MatrixSum<T>::operator*=(const T c) {
	for (Int k = 0; k < matrs.size(); k++) {
		if (inverse[k])
			*matrs[k] *= 1.0 / c;
		else
			*matrs[k] *= c;
	}
}

template <class T>
Vector<T> MatrixSum<T>::operator*(const Vector<T> &vec) const {
	assert(vec.size() == m);
	Vector<T> res(n);
	res.fill(T());
	Vector<T> tmp;
	for (Int k = 0; k < matrs.size(); k++)
		if (inverse[k]) {
			tmp = vec;
			matrs[k]->solve(tmp);
			res += tmp;
		} else {
			res += *matrs[k] * vec;
		}
	return res;
}

template <class T>
double MatrixSum<T>::sizeGb() const {
	assert(false);
}