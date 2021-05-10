#pragma once

#include "AMatrix.h"

//Square with diagonal blocks of different size
template <class T>
class BlockDiagMatr :
	public AMatrix<T> {
public:
	std::vector<std::shared_ptr<AMatrix<T>>> blocks;
	std::vector<bool> inverse; //true if block is an inverse matrix
	std::vector<Int> ns;
	using AMatrix<T>::n;
	using AMatrix<T>::m;
	BlockDiagMatr(const std::vector<std::shared_ptr<AMatrix<T>>> &blocks);
	BlockDiagMatr(const std::vector<std::shared_ptr<AMatrix<T>>> &blocks, \
		const std::vector<bool> &inverse);
	BlockDiagMatr(const BlockDiagMatr<T> &rhs) = default;
	BlockDiagMatr(BlockDiagMatr<T> &&rhs) = default;
	BlockDiagMatr<T> & operator=(const BlockDiagMatr<T> &rhs);
	BlockDiagMatr<T> & operator=(const AMatrix<T> &rhs) override;
	BlockDiagMatr<T> & operator=(BlockDiagMatr<T> &&rhs);
	BlockDiagMatr<T> & operator=(AMatrix<T> &&rhs) override;
	void fill(const T val) override;
	void resize(const Int nnew) override;
	void solve(Vector<T> & rhssol) override;
	void print() const override;
	void write(const string &filename) const override;
	BlockDiagMatr<T> & operator*=(const T c) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	double sizeGb() const override;
	~BlockDiagMatr(void) = default;
};

