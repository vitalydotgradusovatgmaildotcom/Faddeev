#pragma once

#include "AMatrix.h"

//Square matrix
template <class T>
class IdentityM :
	public AMatrix<T> {
public:
	IdentityM(const Int n);
	IdentityM(const IdentityM &rhs) = default;
	IdentityM(IdentityM &&rhs) = default;
	IdentityM<T> & operator=(const IdentityM<T> &rhs);
	IdentityM<T> & operator=(const AMatrix<T> &rhs) override;
	IdentityM<T> & operator=(IdentityM<T> &&rhs);
	IdentityM<T> & operator=(AMatrix<T> &&rhs) override;
	void fill(const T val) override;
	void resize(const Int nnew) override;
	void solve(Vector<T> & rhssol) override;
	void print() const override;
	void write(const string &filename) const override;
	IdentityM<T> & operator*=(const T c) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	double sizeGb() const override;
	~IdentityM(void) = default;
protected:
	using AMatrix<T>::n;
	using AMatrix<T>::m;
	//using AMatrix<T>::garbage;
};

