#pragma once

#include "AMatrix.h"

//The permutation is specified by its full form
//(1 2 ... n)
//(   perm  )
template <class T>
class PermutationM :
	public AMatrix<T> {
public:
	PermutationM(void);
	PermutationM(std::vector<Int> &&perm);
	PermutationM(const PermutationM<T> &rhs) = default;
	PermutationM(PermutationM<T> &&rhs) = default;
	PermutationM<T> & operator=(const PermutationM<T> &rhs);
	PermutationM<T> & operator=(const AMatrix<T> &rhs) override;
	PermutationM<T> & operator=(PermutationM<T> &&rhs);
	PermutationM<T> & operator=(AMatrix<T> &&rhs) override;
	void fill(const T val) override;
	void resize(const Int nnew) override;
	void solve(Vector<T> & rhssol) override;
	void print() const override;
	void write(const string &filename) const override;
	PermutationM<T> & operator*=(const T c) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	double sizeGb() const override;
	~PermutationM(void) = default;
protected:
	using AMatrix<T>::n;
	using AMatrix<T>::m;
	std::vector<Int> p;
};

