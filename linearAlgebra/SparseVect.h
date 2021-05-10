#pragma once

#include "Vector.h"

template <class T>
class SparseVect {
public:
	T *v;
	Int *inds;
	SparseVect(const Int n, const Int nz);
	SparseVect(const Vector<T> &u, \
		const std::vector<Int> &inds_vec);
	SparseVect(const SparseVect &rhs);
	SparseVect(SparseVect &&rhs);
	SparseVect<T> & operator=(const SparseVect<T> &rhs);
	SparseVect<T> & operator=(SparseVect<T> &&rhs);
	T scal(const Vector<T> &u) const; //returns the scalar product (this,u)
	T dotp(const Vector<T> &u) const; //returns the vector-vector dot product this*u
	double norm() const; // calculates Euclidean norm
	void normalize(); // normalization of the vector
	SparseVect<T> operator*(const T x) const; //calculates (scalar x)*this
	Vector<T> operator+(const Vector<T> &u) const; // calculates this + u
	SparseVect<T> & operator*=(const T x); //multiplies this by x
	void resize(const Int nnew, const Int nznew);
	void fill(const T val);
	Int size() const;
	Int nnonz() const; //number of nonzero elements
	void print() const;
	void write(const string &filename) const;
	Vector<T> toFullForm() const;
	~SparseVect(void);
protected:
	Int n, nz;
};

