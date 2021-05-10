#pragma once

#include "AMatrix.h"

//Square with diagonal blocks of different size 
template <class T>
class BlockMatr :
	public AMatrix<T> {
public:
	vector<vector<shared_ptr<AMatrix<T>>>> blocks;
	vector<vector<bool>> inverse;
	using AMatrix<T>::n;
	using AMatrix<T>::m;
	std::vector<Int> ns;
	Int nBlock;
	BlockMatr( \
		const vector<vector<shared_ptr<AMatrix<T>>>> &blocks);
	BlockMatr(\
		const vector<vector<shared_ptr<AMatrix<T>>>> &blocks, \
			const vector<vector<bool>> inverse);
	BlockMatr(const BlockMatr<T> &rhs) = default;
	BlockMatr(BlockMatr<T> &&rhs) = default;
	BlockMatr<T> & operator=(const BlockMatr<T> &rhs);
	BlockMatr<T> & operator=(const AMatrix<T> &rhs) override;
	BlockMatr<T> & operator=(BlockMatr<T> &&rhs);
	BlockMatr<T> & operator=(AMatrix<T> &&rhs) override;
	void fill(const T val) override;
	void resize(const Int nnew) override;
	void solve(Vector<T> & rhssol) override;
	void print() const override;
	void write(const string &filename) const override;
	BlockMatr<T> & operator*=(const T c) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	double sizeGb() const override;
	~BlockMatr() = default;
};

