#pragma once

#include "AMatrix.h"

//Square matrix
template <class T>
class TriDiagM :
	public AMatrix<T> {
public:
	Vector<T> d, u, l;
	TriDiagM(void);
	TriDiagM(const Int n);
	TriDiagM(const TriDiagM<T> &rhs) = default;
	TriDiagM(TriDiagM<T> &&rhs);
	TriDiagM<T> & operator=(const TriDiagM<T> &rhs);
	TriDiagM<T> & operator=(const AMatrix<T> &rhs) override;
	TriDiagM<T> & operator=(TriDiagM<T> &&rhs);
	TriDiagM<T> & operator=(AMatrix<T> &&rhs) override;
	void fill(const T val) override;
	void resize(const Int nnew) override;
	void solve(Vector<T> & rhssol) override;
	void print() const override;
	void write(const string &filename) const override;
	TriDiagM<T> & operator*=(const T c) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	inline const T & get(Int i, Int j) const;
	inline void set(Int i, Int j, T val);
	double sizeGb() const override;
	~TriDiagM(void) = default;
protected:
	using AMatrix<T>::n;
	using AMatrix<T>::m;
	using AMatrix<T>::garbage;
	bool LUFactorized = false;
	std::vector<Int> ipiv;
	Vector<T> u2; //for LU factorization
};

template <class T>
inline const T & TriDiagM<T>::get(Int i, Int j) const {
	assert(i <= j+1 && j-1 <= i);
	if (i == j)
		return d[i];
	else if (i == j+1)
		return l[i-1];
	else //i == j-1
		return u[i];
}

template <class T>
inline void TriDiagM<T>::set(Int i, Int j, T val) {
	assert(i <= j+1 && j-1 <= i);
	if (i == j)
		d[i] = val;
	else if (i == j+1)
		l[i-1] = val;
	else //i == j-1
		u[i] = val;
}