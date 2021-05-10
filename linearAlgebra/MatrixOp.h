#pragma once

#include "GenMatrix.h"
#include "DiagMatrix.h"
#include "TensorProd.h"
#include "SparseMatr.h"

template <class T>
using TensorProd3d = TensorProd<T, 3>;

template <class T>
GenMatrix<T> operator*(const GenMatrix<T> &a, const GenMatrix<T> &b);

template <class T>
GenMatrix<T> operator*(const GenMatrix<T> &a, const DiagMatrix<T> &b);

template <class T>
GenMatrix<T> operator*(const DiagMatrix<T> &a, const GenMatrix<T> &b);

template <class T>
void genEEV(GenMatrix<T> &a, GenMatrix<T> &b, \
		GenMatrix<Complex> &wl, Vector<Complex> &ev, \
			GenMatrix<Complex> &wr, const bool computeLeft, \
				const bool computeRight);

template <class T>
void simultDiagonalize(GenMatrix<T> &a, GenMatrix<T> &b, \
	GenMatrix<Complex> &wbar, DiagMatrix<Complex> &lambda, \
		GenMatrix<Complex> &w);

//template <class T>
//SparseMatr_old<T> ILU0(const SparseMatr_old<T> &a);

template <class T>
SparseMatr<T> ILU0(const SparseMatr<T> &a);
