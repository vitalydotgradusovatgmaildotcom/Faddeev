#pragma once

#include "AMatrix.h"

//Non square matrix
template <class T>
class GenMatrix :
	public AMatrix<T> {
public:
	T **matr;
	GenMatrix();
	GenMatrix(const Int n, const Int m);
	GenMatrix(const Int n); //square matrix
	GenMatrix(const GenMatrix<T> &rhs);
	GenMatrix(GenMatrix<T> &&rhs);
	GenMatrix(const GenMatrix<T> &rhs, const Int i1, \
		const Int i2, const Int j1, const Int j2);
	GenMatrix<T> & operator=(const GenMatrix<T> &rhs);
	GenMatrix<T> & operator=(const AMatrix<T> &rhs) override;
	GenMatrix<T> & operator=(GenMatrix<T> &&rhs);
	GenMatrix<T> & operator=(AMatrix<T> &&rhs) override;
	void resize(const Int nnew, const Int mnew);
	void resize(const Int nnew) override;
	void fill(const T val) override;
	void solve(Vector<T> & rhssol) override;
	inline T * operator[](const Int i); //i'th row
	inline const T * operator[](const Int i) const;
	void print() const override;
	void write(const string &filename) const override;
	GenMatrix<T> & operator*=(const T c) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	double sizeGb() const override;
	void getEEV(Vector<Complex> &ev, GenMatrix<Complex> &wl, \
				GenMatrix<Complex> &wr, const bool computeLeft, \
				const bool computeRight);
	//NB! returns right eigenvectors in ROWS of wr
	// left eigenvectors in ROWS of wl
	//eigenvectors normalized to unit 2-norm and largest component real
	void cconj();
	void trans();
	void inv();
	~GenMatrix(void);
	template <class T> friend void genEEV(GenMatrix<T> &, GenMatrix<T> &, \
			GenMatrix<Complex> &, Vector<Complex> &, \
				GenMatrix<Complex> &, const bool computeLeft, \
					const bool computeRight);
protected:
	using AMatrix<T>::n;
	using AMatrix<T>::garbage;
	using AMatrix<T>::m;
	bool LUFactorized = false;
	std::vector<Int> ipiv;
};

template <class T>
inline T * GenMatrix<T>::operator[](const Int i) {
	return matr[i];
}

template <class T>
inline const T * GenMatrix<T>::operator[](const Int i) const {
	return matr[i];
}

