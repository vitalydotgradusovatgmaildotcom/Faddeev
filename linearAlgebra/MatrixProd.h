#pragma once

#include "AMatrix.h"

//Product of non square matrices
template <class T>
class MatrixProd :
	public AMatrix<T> {
public:
	std::vector<std::shared_ptr<AMatrix<T>>> matrs;
	std::vector<bool> inverse; //true if inverse matrix in product
	using AMatrix<T>::n;
	using AMatrix<T>::m;
	MatrixProd(const std::vector<std::shared_ptr<AMatrix<T>>> &matrs);
	MatrixProd(const std::vector<std::shared_ptr<AMatrix<T>>> &matrs, \
				const std::vector<bool> &inverse);
	MatrixProd(const MatrixProd<T> &rhs) = default;
	MatrixProd(MatrixProd<T> &&rhs) = default;
	MatrixProd<T> & operator=(const MatrixProd<T> &rhs);
	MatrixProd<T> & operator=(const AMatrix<T> &rhs) override;
	MatrixProd<T> & operator=(MatrixProd<T> &&rhs);
	MatrixProd<T> & operator=(AMatrix<T> &&rhs) override;
	void fill(const T val) override;
	void resize(const Int nnew) override;
	void solve(Vector<T> & rhssol) override;
	void print() const override;
	void write(const string &filename) const override;
	MatrixProd<T> & operator*=(const T c) override;
	//AMatrix<T> & operator+=(const AMatrix<T> &b) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	double sizeGb() const override;
	~MatrixProd(void) = default;
};