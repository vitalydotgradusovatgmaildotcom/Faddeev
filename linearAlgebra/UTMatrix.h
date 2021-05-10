#pragma once

#include "AMatrix.h"

//Square matrix
template <class T>
class UTMatrix :
	public AMatrix<T> {
public:
	T **matr; //NB row-major layout
	UTMatrix(void);
	UTMatrix(const Int n);
	UTMatrix(const UTMatrix<T> &rhs);
	UTMatrix(UTMatrix<T> &&rhs);
	UTMatrix<T> & operator=(const UTMatrix<T> &rhs);
	UTMatrix<T> & operator=(const AMatrix<T> &rhs) override;
	UTMatrix<T> & operator=(UTMatrix<T> &&rhs);
	UTMatrix<T> & operator=(AMatrix<T> &&rhs) override;
	void resize(const Int nnew) override;
	void fill(const T val) override;
	void solve(Vector<T> & rhssol) /*const*/ override;
	//inline T * operator[](const Int i); //i'th row
	//inline const T * operator[](const Int i) const;
	inline const T & get(Int i, Int j) const;
	inline void set(Int i, Int j, T val);
	void print() const override;
	void write(const string &filename) const override;
	UTMatrix<T> & operator*=(const T c) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	double sizeGb() const override;
	~UTMatrix(void);
protected:
	using AMatrix<T>::n;
	using AMatrix<T>::m;
	using AMatrix<T>::garbage;
	Int nnonz;
};

template <class T>
inline const T & UTMatrix<T>::get(Int i, Int j) const {
	assert(i <= j);
	return matr[i][j-i];
}

template <class T>
inline void UTMatrix<T>::set(Int i, Int j, T val) {
	assert(i <= j);
	matr[i][j-i] = val;
}