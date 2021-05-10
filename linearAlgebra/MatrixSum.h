#pragma once

#include "AMatrix.h"

//Sum of non square matrices
template <class T>
class MatrixSum :
	public AMatrix<T> {
public:
	std::vector<std::shared_ptr<AMatrix<T>>> matrs;
	std::vector<bool> inverse; //true if inverse matrix in sum
	using AMatrix<T>::n;
	using AMatrix<T>::m;
	MatrixSum(const std::vector<std::shared_ptr<AMatrix<T>>> &matrs);
	MatrixSum(const std::vector<std::shared_ptr<AMatrix<T>>> &matrs, \
		const std::vector<bool> &inverse);
	MatrixSum(const MatrixSum<T> &rhs) = default;
	MatrixSum(MatrixSum<T> &&rhs) = default;
	MatrixSum<T> & operator=(const MatrixSum<T> &rhs);
	MatrixSum<T> & operator=(const AMatrix<T> &rhs) override;
	MatrixSum<T> & operator=(MatrixSum<T> &&rhs);
	MatrixSum<T> & operator=(AMatrix<T> &&rhs) override;
	void fill(const T val) override;
	void resize(const Int nnew) override;
	void solve(Vector<T> & rhssol) override;
	void print() const override;
	void write(const string &filename) const override;
	MatrixSum<T> & operator*=(const T c) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	double sizeGb() const override;
	~MatrixSum() = default;
};

