#pragma once

#include "cph.h"

template <class T>
class Vector {
public:
	T *v;
	Vector();
	Vector(const Int n);
	Vector(const Int n, const T val);
	Vector(const Vector &rhs);
	Vector(Vector &&rhs);
	Vector<T> & operator=(const Vector<T> &rhs); //assignment
	Vector<T> & operator=(Vector<T> &&rhs); //move assignment
	T scal(const Vector<T> &u) const; //returns the scalar product (this,u)
	T dotp(const Vector<T> &u) const; //returns the vector-vector dot product this*u
	double norm2() const; // calculates (this,this)
	double norm() const; // calculates Euclidean norm
	void normalize(); // normalization of the vector
	Vector<T> operator*(const T x) const; //calculates (scalar x)*this
	Vector<T> operator+(const Vector<T> &u) const; // calculates this + u
	Vector<T> operator-(const Vector<T> &u) const; // calculates this - u
	Vector<T> & operator*=(const T x); //multiplies this by x
	Vector<T> & operator+=(const Vector<T> &u); //adds u to this
	Vector<T> & operator-=(const Vector<T> &u); //subtracts u from this
	void resize(const Int nnew);
	void setZero(); //sets all the elements eq to zero
	void fill(const T val);
	Int size() const;
	void print() const;
	void write(const string &filename) const;
	std::vector<Vector<T>> makePartition(const std::vector<Int> &inds) const;
	void conj();
	inline T & operator[](const Int i);
	inline const T & operator[](const Int i) const;
	~Vector(void);
	//void setRandom() final;
	//Vector(const string &filename);
	//void dump(const string &filename) const;
protected:
	Int n;
	bool owner = true;
	Vector(const Int n, T *vals);
};

template <class T>
inline T & Vector<T>::operator[](const Int i) {
	return v[i];
}

template <class T>
inline const T & Vector<T>::operator[](const Int i) const {
	return v[i];
}