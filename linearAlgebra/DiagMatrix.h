#pragma once

#include "AMatrix.h"

//Square matrix
template <class T>
class DiagMatrix :
	public AMatrix<T> {
public:
	Vector<T> d;
	DiagMatrix(void);
	DiagMatrix(const Int n);
	DiagMatrix(const DiagMatrix<T> &rhs) = default;
	DiagMatrix(DiagMatrix<T> &&rhs) = default;
	DiagMatrix<T> & operator=(const DiagMatrix<T> &rhs);
	DiagMatrix<T> & operator=(const AMatrix<T> &rhs) override;
	DiagMatrix<T> & operator=(DiagMatrix<T> &&rhs);
	DiagMatrix<T> & operator=(AMatrix<T> &&rhs) override;
	void fill(const T val) override;
	void resize(const Int nnew) override;
	void solve(Vector<T> & rhssol) override;
	void print() const override;
	void write(const string &filename) const override;
	DiagMatrix<T> & operator*=(const T c) override;
	//DiagMatrix<T> & operator+=(const AMatrix<T> &b) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	inline T & operator[](const Int i); //i'th diagonal element
	inline const T & operator[](const Int i) const;
	double sizeGb() const override;
	~DiagMatrix(void) = default;
protected:
	using AMatrix<T>::n;
	using AMatrix<T>::m;
	using AMatrix<T>::garbage;
};

template <class T>
inline T & DiagMatrix<T>::operator[](const Int i) {
	return d[i];
}

template <class T>
inline const T & DiagMatrix<T>::operator[](const Int i) const {
	return d[i];
}