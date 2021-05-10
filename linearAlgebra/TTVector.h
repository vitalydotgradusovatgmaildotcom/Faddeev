#pragma once

#include "cph.h"
#include "MatrixOp.h"

template <class T, Int d>
class TTVector {
public:
	array<vector<GenMatrix<T>>, d> us;
	TTVector();
	TTVector(const array<Int, d> &ns, const array<Int, d-1> &rs);
	TTVector(const TTVector<T, d> &rhs);
	TTVector(TTVector<T, d> &&rhs);
	TTVector<T, d> & operator=(const TTVector<T, d> &rhs);
	TTVector<T, d> & operator=(TTVector<T, d> &&rhs);
	T scal(const TTVector<T, d> &u) const; //returns the scalar product (this,u)
	T dotp(const TTVector<T, d> &u) const; //returns the vector-vector dot product this*u
	double norm2() const; // calculates (this,this)
	double norm() const; // calculates Euclidean norm
	void normalize(); // normalization of the vector
	TTVector<T, d> operator*(const T x) const; //calculates (scalar x)*this
	TTVector<T, d> operator+(const TTVector<T, d> &u) const; // calculates this + u
	TTVector<T, d> operator-(const TTVector<T, d> &u) const; // calculates this - u
	TTVector<T, d> & operator*=(const T x); //multiplies this by x
	TTVector<T, d> & operator+=(const TTVector<T, d> &u); //adds u to this
	TTVector<T, d> & operator-=(const TTVector<T, d> &u); //subtracts u from this
	array<Int, d> size() const;
	array<Int, d - 1> ranks() const;
	void print() const;
	void write(const string &filename) const;
	T & operator[](const Int i);
	const T & operator[](const Int i) const;
	~TTVector();
};

template <class T, Int d>
TTVector<T, d>::TTVector() {

}

template <class T, Int d>
TTVector<T, d>::TTVector(const array<Int, d> &ns, const array<Int, d - 1> &rs) {

}

template <class T, Int d>
TTVector<T, d>::TTVector(const TTVector<T, d> &rhs) {

}

template <class T, Int d>
TTVector<T, d>::TTVector(TTVector<T, d> &&rhs) {

}

template <class T, Int d>
TTVector<T, d> & TTVector<T, d>::operator=(const TTVector<T, d> &rhs) {
	return *this;
}

template <class T, Int d>
TTVector<T, d> & TTVector<T, d>::operator=(TTVector<T, d> &&rhs) {
	return *this;
}

template <class T, Int d>
T TTVector<T, d>::scal(const TTVector<T, d> &u) const {
	//returns the scalar product (this,u)
	return T();
}

template <class T, Int d>
T TTVector<T, d>::dotp(const TTVector<T, d> &u) const {
	//returns the vector-vector dot product this*u
	return T();
}

template <class T, Int d>
double TTVector<T, d>::norm2() const {
	// calculates (this,this)
	return 0.0;
}

template <class T, Int d>
double TTVector<T, d>::norm() const {
	// calculates Euclidean norm
	return 0.0;
}

template <class T, Int d>
void TTVector<T, d>::normalize() {
	// normalization of the vector
}

template <class T, Int d>
TTVector<T, d> TTVector<T, d>::operator*(const T x) const {
	//calculates (scalar x)*this
	return TTVector<T, d>();
}

template <class T, Int d>
TTVector<T, d> TTVector<T, d>::operator+(const TTVector<T, d> &u) const {
	// calculates this + u
	return TTVector<T, d>();
}

template <class T, Int d>
TTVector<T, d> TTVector<T, d>::operator-(const TTVector<T, d> &u) const {
	// calculates this - u
	return TTVector<T, d>();
}

template <class T, Int d>
TTVector<T, d> & TTVector<T, d>::operator*=(const T x) {
	//multiplies this by x
	return *this;
}

template <class T, Int d>
TTVector<T, d> & TTVector<T, d>::operator+=(const TTVector<T, d> &u) {
	//adds u to this
	return *this;
}

template <class T, Int d>
TTVector<T, d> & TTVector<T, d>::operator-=(const TTVector<T, d> &u) {
	//subtracts u from this
	return *this;
}

template <class T, Int d>
array<Int, d> TTVector<T, d>::size() const {
	return array<Int, d>();
}

template <class T, Int d>
array<Int, d - 1> TTVector<T, d>::ranks() const {
	return array<Int, d - 1>();
}

template <class T, Int d>
void TTVector<T, d>::print() const {

}

template <class T, Int d>
void TTVector<T, d>::write(const string &filename) const {

}

template <class T, Int d>
T & TTVector<T, d>::operator[](const Int i) {
	return T();
}

template <class T, Int d>
const T & TTVector<T, d>::operator[](const Int i) const {
	return T();
}

template <class T, Int d>
TTVector<T, d>::~TTVector() {

}