
#include "MatrixOp.h"

template GenMatrix<double> operator*(const GenMatrix<double> &, const DiagMatrix<double> &);
template GenMatrix<Complex> operator*(const GenMatrix<Complex> &, const DiagMatrix<Complex> &);
template GenMatrix<double> operator*(const DiagMatrix<double> &, const GenMatrix<double> &);
template GenMatrix<Complex> operator*(const DiagMatrix<Complex> &, const GenMatrix<Complex> &);
template void genEEV(GenMatrix<double> &, GenMatrix<double> &, GenMatrix<Complex> &, \
				Vector<Complex> &, GenMatrix<Complex> &, const bool, const bool);
template void genEEV(GenMatrix<Complex> &, GenMatrix<Complex> &, GenMatrix<Complex> &, \
				Vector<Complex> &, GenMatrix<Complex> &, const bool, const bool);
template void simultDiagonalize(GenMatrix<double> &a, GenMatrix<double> &b, \
	GenMatrix<Complex> &wbar, DiagMatrix<Complex> &lambda, \
		GenMatrix<Complex> &w);
template void simultDiagonalize(GenMatrix<Complex> &a, GenMatrix<Complex> &b, \
	GenMatrix<Complex> &wbar, DiagMatrix<Complex> &lambda, \
		GenMatrix<Complex> &w);
//template SparseMatr_old<double> ILU0(const SparseMatr_old<double> &a);
//template SparseMatr_old<Complex> ILU0(const SparseMatr_old<Complex> &a);
template SparseMatr<double> ILU0(const SparseMatr<double> &a);
template SparseMatr<Complex> ILU0(const SparseMatr<Complex> &a);

template <>
GenMatrix<double> operator*(const GenMatrix<double> &a, const GenMatrix<double> &b) {
	//assert(!a.garbage);
	Int an = a.nrows(); Int am = a.ncols();
	Int bn = b.nrows(); Int bm = b.ncols();
	assert(am == bn);
	GenMatrix<double> res(an, bm);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, an, bm, \
		am, 1.0, a.matr[0], am, b.matr[0], bm, 0.0, res.matr[0], bm);
	return res;
}

template <>
GenMatrix<Complex> operator*(const GenMatrix<Complex> &a, const GenMatrix<Complex> &b) {
	//assert(!a.garbage);
	Int an = a.nrows(); Int am = a.ncols();
	Int bn = b.nrows(); Int bm = b.ncols();
	assert(am == bn);
	GenMatrix<Complex> res(an, bm);
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, an, bm, \
		am, &zone, a.matr[0], am, b.matr[0], bm, &zzero, res.matr[0], bm);
	return res;
}

template <class T>
GenMatrix<T> operator*(const GenMatrix<T> &a, const DiagMatrix<T> &b) {
	Int n = a.nrows(); Int m = b.nrows();
	assert(a.ncols() == m);
	GenMatrix<T> res(n, m);
	for (Int i = 0; i < n; i++)
		for (Int j = 0; j < m; j++)
			res[i][j] = a[i][j]*b[j];
	return res;
}

template <class T>
GenMatrix<T> operator*(const DiagMatrix<T> &a, const GenMatrix<T> &b) {
	Int n = a.nrows(); Int m = b.ncols();
	assert(n == b.nrows());
	GenMatrix<T> res(n, m);
	for (Int i = 0; i < n; i++)
		for (Int j = 0; j < m; j++)
			res[i][j] = a[i]*b[i][j];
	return res;
}

template<>
void genEEV<double>(GenMatrix<double> &a, GenMatrix<double> &b, \
		GenMatrix<Complex> &wl, Vector<Complex> &ev, \
			GenMatrix<Complex> &wr, const bool computeLeft, \
				const bool computeRight) {
//NB! returns right eigenvectors in COLUMNS of wr
// CONJUGATED left eigenvectors in COLUMNS of wl
	assert(!a.garbage); assert(!b.garbage);
	assert(a.n == a.m); assert(b.n == b.m);
	assert(a.n == b.n);

	ev.resize(0); wl.resize(0); wr.resize(0);
	Vector<double> alphai(a.n), alphar(a.n), beta(a.n);
	GenMatrix<double> wl_(0), wr_(0);
	Int ldwl = computeLeft ? a.n : 1;
	wl_.resize(a.n, ldwl);
	Int ldwr = computeRight ? a.n : 1;
	wr_.resize(a.n, ldwr);

	Int info = LAPACKE_dggev(LAPACK_ROW_MAJOR, computeLeft ? 'V' : 'N', \
				computeRight ? 'V' : 'N', a.n, a.matr[0], a.n,
				b.matr[0], b.n, alphar.v, alphai.v, beta.v, \
				wl_.matr[0], ldwl, wr_.matr[0], ldwr);
	assert(info == 0);

	ev.resize(a.n);
	for (Int i = 0; i < a.n; i++) {
		//Checks must be done in fact
		//beta[i] = 0, underflow or overflow can happen
		ev[i] = Complex(alphar[i], alphai[i])/beta[i];
	}
	alphar.resize(0); alphai.resize(0); beta.resize(0);

	if (computeLeft) {
		wl.resize(a.n);
		for (Int j = 0; j < a.n; j++) {
			if ( IS_EPS(ev[j].imag() / ev[j].real()) ) {
				for (Int i = 0; i < a.n; i++)
					wl[i][j] = wl_[i][j];
			} else { //complex pair
				for (Int i = 0; i < a.n; i++) {
					wl[i][j] = Complex(wl_[i][j], wl_[i][j+1]);
					wl[i][j+1] = Complex(wl_[i][j], -wl_[i][j+1]);
				}
				j++;
			}
		}
		wl_.resize(0);
	}

	if (computeRight) {
		wr.resize(a.n);
		for (Int j = 0; j < a.n; j++) {
			if ( IS_EPS(ev[j].imag() / ev[j].real()) ) {
				for (Int i = 0; i < a.n; i++)
					wr[i][j] = wr_[i][j];
			} else { //complex pair
				for (Int i = 0; i < a.n; i++) {
					wr[i][j] = Complex(wr_[i][j], wr_[i][j+1]);
					wr[i][j+1] = Complex(wr_[i][j], -wr_[i][j+1]);
				}
				j++;
			}
		}
		wr_.resize(0);
	}

	a.garbage = true; b.garbage = true;
}

template<>
void genEEV<Complex>(GenMatrix<Complex> &a, GenMatrix<Complex> &b, \
		GenMatrix<Complex> &wl, Vector<Complex> &ev, \
			GenMatrix<Complex> &wr, const bool computeLeft, \
				const bool computeRight) {
//NB! returns right eigenvectors in COLUMNS of wr
// CONJUGATED left eigenvectors in COLUMNS of wl

	assert(!a.garbage); assert(!b.garbage);
	assert(a.n == a.m); assert(b.n == b.m);
	assert(a.n == b.n);

	ev.resize(0);
	Int ldwl = computeLeft ? a.n : 1;
	wl.resize(a.n, ldwl);
	Int ldwr = computeRight ? a.n : 1;
	wr.resize(a.n, ldwr);

	Vector<Complex> alpha(a.n), beta(a.n);
	Int info = LAPACKE_zggev(LAPACK_ROW_MAJOR, computeLeft ? 'V' : 'N', \
				computeRight ? 'V' : 'N', a.n, a.matr[0], a.n, \
				b.matr[0], a.n, alpha.v, beta.v, \
				wl.matr[0], ldwl, wr.matr[0], ldwr);
	assert(info == 0);

	ev.resize(a.n);
	for (Int i = 0; i < a.n; i++) {
		//Checks must be done in fact
		//beta[i] = 0, underflow or overflow can happen
		ev[i] = alpha[i]/beta[i];
	}
	if (!computeLeft) wl.resize(0);
	if (!computeRight) wr.resize(0);

	a.garbage = true; b.garbage = true;
}

template <class T>
void simultDiagonalize(GenMatrix<T> &a, GenMatrix<T> &b, \
	GenMatrix<Complex> &wbar, DiagMatrix<Complex> &lambda, \
		GenMatrix<Complex> &w) {
//return matrices Wbar and W such that
//(Wbar)*A*W = Lambda, (Wbar)*B*W = I
	Vector<Complex> ev;
	Int n = b.nrows();
	assert(n == b.ncols());
	GenMatrix<Complex> b_copy(n);
	for (Int i = 0; i < n; i++)
		for (Int j = 0; j < n; j++)
			b_copy[i][j] = b[i][j];
	genEEV(a, b, wbar, ev, w, true, true);
	wbar.cconj();
	//now (Wbar)*A*W = Lambda_1, (Wbar)*B*W = Lambda_2
	//calculate (Lambda_2)^(-1)
	lambda.resize(n);
	Vector<Complex> tmp(n), tmp2(n);
	for (Int i = 0; i < n; i++) {
		for (Int j = 0; j < n; j++) {
			tmp[j] = w[j][i];
			tmp2[j] = wbar[i][j];
		}
		lambda[i] = zone/tmp2.dotp(b_copy * tmp);
	}
	b_copy.resize(0); tmp.resize(0); tmp2.resize(0);
	wbar = lambda * wbar;
	for (Int i = 0; i < n; i++)
		lambda[i] = ev[i];
}

/*
template <class T>
SparseMatr_old<T> ILU0(const SparseMatr_old<T> &a) {
	static const double PIV_EPS = 1.0e-16; //as in MKL realization
	static const double PIV_MIN = 1.0e-10; //as in MKL realization
	SparseMatr_old<T> lu;
	Int n, i, k, j, kj;
	Int nz, kInd, kjInd;
	std::vector<Int> diagInds;
	T lukj;

	n = a.n; nz = a.nz;

	lu.n = n; lu.nz = nz;
	lu.garbage = true; lu.LUFactorized = true;
	lu.ja = a.ja; lu.ia = a.ia;
	lu.owner = false;
	lu.aa = nz > 0 ? (T *)MALLOC(nz * sizeof(T)) : nullptr;
	std::copy(a.aa, a.aa+nz, lu.aa);

	//find indices of diagonal elements in b.aa, b.ja
	diagInds.reserve(n);
	for (Int i = 0; i < n; i++)
		for (kInd = lu.ia[i]; kInd < lu.ia[i + 1]; kInd++)
			if (lu.ja[kInd] == i)
				diagInds.push_back(kInd);

	//make ilu0 decomposition
	for (Int i = 1; i < n; i++)
		for (kInd = lu.ia[i]; kInd < diagInds[i]; kInd++) { //k = 0, ..., i-1
			k = lu.ja[kInd];
			if (abs(lu.aa[diagInds[k]]) < PIV_EPS) {
				cout << "NB! Zero pivot in ILU0" << endl;
				lu.aa[diagInds[k]] = PIV_MIN;
			} //MKL?
			//assert(abs(lu.aa[diagInds[k]]) == 0.0); //Saad
			lu.aa[kInd] /= lu.aa[diagInds[k]]; //LU(i, k) = LU(i, k) / LU(k, k)
			kjInd = lu.ia[k]; kj = lu.ja[kjInd];
#pragma omp parallel for shared(lu, kInd) private(j, lukj) firstprivate(kj, kjInd) schedule(guided, 2)
			for (Int jInd = kInd + 1; jInd < lu.ia[i + 1]; jInd++) { //j = k+1, ..., n
				j = lu.ja[jInd];
				//find column index kj of first nonzero element
				//in row k with kj >= j (if exists, otherwise last element)
				while (kj < j && kjInd < lu.ia[k+1]-1) {
					kjInd++;
					kj = lu.ja[kjInd];
				}
				if (kj == j)
					lukj = lu.aa[kjInd];
				else
					lukj = T();
				lu.aa[jInd] -= lu.aa[kInd] * lukj;
				//LU(i, j) = LU(i, j) - LU(i, k)*LU(k, j)
			}
		}

	return lu;
}
*/

template <class T>
SparseMatr<T> ILU0(const SparseMatr<T> &a) {
	static const double PIV_EPS = 1.0e-16; //as in MKL realization
	static const double PIV_MIN = 1.0e-10; //as in MKL realization
	SparseMatr<T> lu;
	Int n, i, k, j, kj;
	Int nz, kInd, kjInd;
	std::vector<Int> diagInds;
	T lukj;

	n = a.n; nz = a.nz;

	lu.n = n; lu.nz = nz; lu.m = n;
	lu.garbage = true; lu.LUFactorized = true;
	lu.ja = a.ja; lu.ia = a.ia;
	lu.owner = false;
	lu.aa = nz > 0 ? (T *)MALLOC(nz * sizeof(T)) : nullptr;
	std::copy(a.aa, a.aa + nz, lu.aa);
	lu.createHandleAndDesc();
	lu.descr.type = SPARSE_MATRIX_TYPE_TRIANGULAR;

	//find indices of diagonal elements in b.aa, b.ja
	diagInds.reserve(n);
	for (Int i = 0; i < n; i++)
		for (kInd = lu.ia[i]; kInd < lu.ia[i + 1]; kInd++)
			if (lu.ja[kInd] == i)
				diagInds.push_back(kInd);

	//make ilu0 decomposition
	for (Int i = 1; i < n; i++)
		for (kInd = lu.ia[i]; kInd < diagInds[i]; kInd++) { //k = 0, ..., i-1
			k = lu.ja[kInd];
			if (abs(lu.aa[diagInds[k]]) < PIV_EPS) {
				cout << "NB! Zero pivot in ILU0" << endl;
				lu.aa[diagInds[k]] = PIV_MIN;
			} //MKL?
			//assert(abs(lu.aa[diagInds[k]]) == 0.0); //Saad
			lu.aa[kInd] /= lu.aa[diagInds[k]]; //LU(i, k) = LU(i, k) / LU(k, k)
			kjInd = lu.ia[k]; kj = lu.ja[kjInd];
#pragma omp parallel for shared(lu, kInd) private(j, lukj) firstprivate(kj, kjInd) schedule(guided, 2)
			for (Int jInd = kInd + 1; jInd < lu.ia[i + 1]; jInd++) { //j = k+1, ..., n
				j = lu.ja[jInd];
				//find column index kj of first nonzero element
				//in row k with kj >= j (if exists, otherwise last element)
				while (kj < j && kjInd < lu.ia[k + 1] - 1) {
					kjInd++;
					kj = lu.ja[kjInd];
				}
				if (kj == j)
					lukj = lu.aa[kjInd];
				else
					lukj = T();
				lu.aa[jInd] -= lu.aa[kInd] * lukj;
				//LU(i, j) = LU(i, j) - LU(i, k)*LU(k, j)
			}
		}

	return lu;
}