#pragma once

#include "cph.h"
#include "Vector.h"

template <class T>
class AMatrix {
public:
	AMatrix(void);
	AMatrix(const AMatrix<T> &) = default;
	AMatrix(AMatrix<T> &&rhs);
	virtual AMatrix<T> & operator=(const AMatrix<T> &rhs) = 0;
	virtual AMatrix<T> & operator=(AMatrix<T> &&rhs) = 0;
	virtual void fill(const T val) = 0;
	virtual void resize(const Int nnew) = 0;
	//TODO virtual void resize(const Int nnew, const Int mnew) = 0;
	Int nrows() const;
	Int ncols() const;
	virtual void solve(Vector<T> & rhssol) = 0;
	virtual void print() const = 0;
	virtual void write(const string &filename) const = 0;
	virtual AMatrix<T> & operator*=(const T c) = 0;
	//virtual AMatrix<T> & operator+=(const AMatrix<T> &b) = 0;
	virtual Vector<T> operator*(const Vector<T> &vec) const = 0;
	virtual double sizeGb() const = 0;
	virtual ~AMatrix(void) = default;
protected:
	Int n;
	Int m; //n, m - linear sizes
	bool garbage = false;
};

template <class T>
AMatrix<T>::AMatrix(void) : n(0), m(0) { }

template <class T>
AMatrix<T>::AMatrix(AMatrix<T> &&rhs) \
	: n(rhs.n), m(rhs.m), garbage(rhs.garbage) {
	rhs.n = 0; rhs.m = 0;
	rhs.garbage = false;
}

template <class T>
Int AMatrix<T>::nrows() const {
	return n;
}

template <class T>
Int AMatrix<T>::ncols() const {
	return m;
}