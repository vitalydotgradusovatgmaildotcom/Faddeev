#pragma once

#include "AMatrix.h"

//Non square matrix

template <class T>
class TrivialM :
	public AMatrix<T> {
public:
	TrivialM(const Int n); //square
	TrivialM(const Int n, const Int m);
	TrivialM(const TrivialM &rhs) = default;
	TrivialM(TrivialM &&rhs) = default;
	TrivialM<T> & operator=(const TrivialM<T> &rhs);
	TrivialM<T> & operator=(const AMatrix<T> &rhs) override;
	TrivialM<T> & operator=(TrivialM<T> &&rhs);
	TrivialM<T> & operator=(AMatrix<T> &&rhs) override;
	void fill(const T val) override;
	void resize(const Int nnew) override;
	void resize(const Int nnew, const Int mnew);
	void solve(Vector<T> & rhssol) override;
	void print() const override;
	void write(const string &filename) const override;
	TrivialM<T> & operator*=(const T c) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	double sizeGb() const override;
	~TrivialM() = default;
protected:
	using AMatrix<T>::n;
	using AMatrix<T>::m;
};


template <class T>
TrivialM<T>::TrivialM(const Int n) {
	this->n = n;
	this->m = n;
}

template <class T>
TrivialM<T>::TrivialM(const Int n, const Int m) {
	this->n = n;
	this->m = m;
}

template <class T>
TrivialM<T> & TrivialM<T>::operator=(const TrivialM<T> &rhs) {
	if (this != &rhs) {
		n = rhs.n;
		m = rhs.m;
	}
	return *this;
}

template <class T>
TrivialM<T> & TrivialM<T>::operator=(const AMatrix<T> &rhs) {
	const TrivialM<T> & rhs_ = dynamic_cast<const TrivialM<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
TrivialM<T> & TrivialM<T>::operator=(TrivialM<T> &&rhs) {
	n = rhs.n; rhs.n = 0;
	m = rhs.m; rhs.m = 0;
	return *this;
}

template <class T>
TrivialM<T> & TrivialM<T>::operator=(AMatrix<T> &&rhs) {
	TrivialM<T> &&rhs_ = dynamic_cast<TrivialM<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void TrivialM<T>::fill(const T val) {
	assert(false);
}

template <class T>
void TrivialM<T>::resize(const Int nnew) {
	n = nnew;
}

template <class T>
void TrivialM<T>::resize(const Int nnew, const Int mnew) {
	n = nnew; m = mnew;
}

template <class T>
void TrivialM<T>::solve(Vector<T> & rhssol) {
	assert(false);
}

template <class T>
void TrivialM<T>::print() const {
	cout << "Trivial " << n << " by " << m << " matrix" << endl;
}

template <class T>
void TrivialM<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	f << "Trivial " << n << " by " << m << " matrix" << endl;
	f.close();
}

template <class T>
TrivialM<T> & TrivialM<T>::operator*=(const T c) {
	return *this;
}

template <class T>
Vector<T> TrivialM<T>::operator*(const Vector<T> &vec) const {
	assert(m == vec.size());
	Vector<T> res(n);
	res.fill(T());
	return res;
}

template <class T>
double TrivialM<T>::sizeGb() const {
	return 0.0;
}