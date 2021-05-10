#pragma once

#include "AMatrix.h"
#include "Algorithms.h"

//Non square sparse matrix in CSR format
template <class T>
class SparseMatr :
	public AMatrix<T> {
public:
	T *aa; Int *ja; Int *ia;
	SparseMatr(const Int n, const Int m, const long long int nz);
	SparseMatr(const Int n, const long long int nz); //square matrix
	SparseMatr(const SparseMatr<T> &rhs);
	SparseMatr(SparseMatr<T> &&rhs);
	SparseMatr<T> & operator=(const SparseMatr<T> &rhs);
	SparseMatr<T> & operator=(const AMatrix<T> &rhs) override;
	SparseMatr<T> & operator=(SparseMatr<T> &&rhs);
	SparseMatr<T> & operator=(AMatrix<T> &&rhs) override;
	void fill(const T val) override;
	void resize(const Int nnew) override;
	void solve(Vector<T> & rhssol) override;
	void print() const override;
	void write(const string &filename) const override;
	SparseMatr<T> & operator*=(const T c) override;
	Vector<T> operator*(const Vector<T> &vec) const override;
	void order(); //reorders ja (and aa) so that column numbers in each row are sorted
	inline const T & get(Int i, Int j) const;
	inline void set(Int i, Int j, T val);
	inline Int nnonz() const;
	void fillRandom();
	double sizeGb() const override;
	~SparseMatr();
protected:
	using AMatrix<T>::n;
	using AMatrix<T>::m;
	using AMatrix<T>::garbage;
	Int nz;
	//handle and description
	sparse_matrix_t matr;
	struct matrix_descr descr;
	bool owner = true;
	bool LUFactorized = false;
	//pardiso
	bool pardisoFact = false;
	_MKL_DSS_HANDLE_t pt[64];
	Int mtype, maxfct, mnum;
	Int msglvl;
	Int iparm[64];
	SparseMatr();
	void createHandleAndDesc();
	//pardiso
	void pardisoSetParams();
	inline void solvePardiso(Vector<T> &rhssol);
	void pardisoRelease();
	inline void solveLU(Vector<T> &rhssol);
	template <class T> friend SparseMatr<T> ILU0(const SparseMatr<T> &);
};


template <class T>
void orderAndSqueeze(Int &n, T *aa, Int *ja) {
	//sort
	Int *ja_begin, *ja_end;
	T *aa_begin, *aa_end;
	ja_begin = ja;
	ja_end = ja + n;
	aa_begin = aa;
	aa_end = aa + n;
	sortSimult(ja_begin, ja_end, aa_begin, aa_end, \
		[&](Int i, Int j) {return ja_begin[i] < ja_begin[j]; });
	//squeeze
	Int ind1 = 0, ind2 = 0;
	Int j;
	while (ind2 < n) {
		j = ja[ind2];
		aa[ind1] = aa[ind2++];
		ja[ind1] = j;
		while ((ind2 < n) && (ja[ind2] == j)) {
			aa[ind1] += aa[ind2];
			ind2++;
		}
		ind1++;
	}
	n = ind1;
}

template <class T>
inline const T & SparseMatr<T>::get(Int i, Int j) const {
	//assert(false) if element (i,j) is not stored
	Int ind, ind1, ind2;
	ind1 = ia[i]; ind2 = ia[i + 1];

	if (j >= ja[ind1]) {
		while (ind1 != ind2 - 1) {
			ind = (ind1 + ind2) / 2;
			if (ja[ind] > j)
				ind2 = ind;
			else
				ind1 = ind;
			// ja[ind1] <= j < ja[ind2] now
		}
		if (j == ja[ind1])
			return aa[ind1];
	}

	assert(false); return aa[0];
}

template <>
inline void SparseMatr<double>::set(Int i, Int j, double val) {
	//assert(false) if element (i,j) is not stored
	sparse_status_t info = \
		mkl_sparse_d_set_value(matr, i, j, val);
	assert(info == SPARSE_STATUS_SUCCESS);
}

template <>
inline void SparseMatr<Complex>::set(Int i, Int j, Complex val) {
	//assert(false) if element (i,j) is not stored
	sparse_status_t info = \
		mkl_sparse_z_set_value(matr, i, j, val);
	assert(info == SPARSE_STATUS_SUCCESS);
}

template <class T>
inline Int SparseMatr<T>::nnonz() const {
	return nz;
}