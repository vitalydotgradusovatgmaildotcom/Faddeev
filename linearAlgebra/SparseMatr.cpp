#include "SparseMatr.h"
#include "RndGenRealUni.h"
#include "RndGenIntUni.h"

template class SparseMatr<double>;
template class SparseMatr<Complex>;

template <class T>
SparseMatr<T>::SparseMatr() : aa(nullptr), ja(nullptr), \
ia(nullptr), nz(0), matr(nullptr) {
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
}

template <class T>
SparseMatr<T>::SparseMatr( \
	const Int n, const Int m, const long long int nz) : nz(nz) {
	this->n = n; this->m = m;
	assert(nz <= maxInt);
	assert((n >= 0) && (m >= 0) \
		&& (nz >= 0) && ((long long int)(n)*m >= nz));

	aa = nz > 0 ? (T *)MALLOC(nz * sizeof(T)) : nullptr;
	ja = nz > 0 ? (Int *)MALLOC(nz * sizeof(Int)) : nullptr;
	ia = n > 0 ? (Int *)MALLOC((n + 1) * sizeof(Int)) : nullptr;
	
	createHandleAndDesc();
	
	/*sparse_index_base_t indexing;
	Int rows, cols;
	Int *rows_start = (Int *)MALLOC(n * sizeof(Int));
	Int *rows_end = (Int *)MALLOC(n * sizeof(Int));
	Int *col_indx = (Int *)MALLOC(nz * sizeof(Int));
	double *values = (double *)MALLOC(nz * sizeof(double));
	info = mkl_sparse_d_export_csr(matr, &indexing, \
		&rows, &cols, &rows_start, &rows_end, \
			&col_indx, &values);*/
	 
}

template <class T>
SparseMatr<T>::SparseMatr(const Int n, const long long int nz) : \
	SparseMatr<T>(n, n, nz) {}

template <class T>
SparseMatr<T>::SparseMatr(const SparseMatr<T> &rhs) : \
		AMatrix<T>(rhs), nz(rhs.nz) {
	assert(!rhs.garbage);
	aa = nz > 0 ? (T *)MALLOC(nz * sizeof(T)) : nullptr;
	ja = nz > 0 ? (Int *)MALLOC(nz * sizeof(Int)) : nullptr;
	ia = n > 0 ? (Int *)MALLOC((n + 1) * sizeof(Int)) : nullptr;
	for (Int ind = 0; ind < nz; ind++) {
		aa[ind] = rhs.aa[ind];
		ja[ind] = rhs.ja[ind];
	}
	for (Int i = 0; i < n + 1; i++) {
		ia[i] = rhs.ia[i];
	}

	createHandleAndDesc();
}

template <class T>
SparseMatr<T>::SparseMatr(SparseMatr<T> &&rhs) : \
			AMatrix<T>(std::move(rhs)), nz(rhs.nz), \
			LUFactorized(rhs.LUFactorized), owner(rhs.owner) {
	if (rhs.pardisoFact) {
		for (Int k = 0; k < 64; k++) {
			pt[k] = rhs.pt[k];
			iparm[k] = rhs.iparm[k];
		}
		mtype = rhs.mtype; maxfct = rhs.maxfct;
		mnum = rhs.mnum; msglvl = rhs.msglvl;
		pardisoFact = true;
	}
	rhs.pardisoFact = false;
	rhs.LUFactorized = false;
	rhs.owner = true;
	aa = rhs.aa; rhs.aa = nullptr;
	ja = rhs.ja; rhs.ja = nullptr;
	ia = rhs.ia; rhs.ia = nullptr;
	matr = rhs.matr; rhs.matr = nullptr;
	descr = rhs.descr;
}

template <class T>
SparseMatr<T> & SparseMatr<T>::operator=(const SparseMatr<T> &rhs) {
	assert(!rhs.garbage);
	if (this != &rhs) {
		pardisoRelease();
		garbage = rhs.garbage;
		owner = true;
		if (nz != rhs.nz) {
			if (aa != nullptr) FREE(aa);
			if (ja != nullptr && owner) FREE(ja);
			nz = rhs.nz;
			aa = nz > 0 ? (T *)MALLOC(nz * sizeof(T)) : nullptr;
			ja = nz > 0 ? (Int *)MALLOC(nz * sizeof(Int)) : nullptr;
		}
		for (Int ind = 0; ind < nz; ind++) {
			aa[ind] = rhs.aa[ind];
			ja[ind] = rhs.ja[ind];
		}
		if (n != rhs.n) {
			if (ia != nullptr && owner) FREE(ia);
			n = rhs.n; m = rhs.m;
			ia = n > 0 ? (Int *)MALLOC((n + 1) * sizeof(Int)) : nullptr;
		}
		for (Int i = 0; i < n + 1; i++)
			ia[i] = rhs.ia[i];
		//recreate handle in any case
		sparse_status_t info = \
			mkl_sparse_destroy(matr);
		assert(info == SPARSE_STATUS_SUCCESS);
		createHandleAndDesc();
	}
	return *this;
}

template <class T>
SparseMatr<T> & SparseMatr<T>::operator=(const AMatrix<T> &rhs) {
	const SparseMatr<T> &rhs_ = dynamic_cast<const SparseMatr<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
SparseMatr<T> & SparseMatr<T>::operator=(SparseMatr<T> &&rhs) {
	pardisoRelease();
	if (rhs.pardisoFact) {
		for (Int k = 0; k < 64; k++) {
			pt[k] = rhs.pt[k];
			iparm[k] = rhs.iparm[k];
		}
		mtype = rhs.mtype; maxfct = rhs.maxfct;
		mnum = rhs.mnum; msglvl = rhs.msglvl;
		pardisoFact = true;
	}
	rhs.pardisoFact = false;
	LUFactorized = rhs.LUFactorized; rhs.LUFactorized = false;
	garbage = rhs.garbage; rhs.garbage = false;
	owner = rhs.owner; rhs.owner = true;
	n = rhs.n; nz = rhs.nz; m = rhs.m;
	if (aa != nullptr) FREE(aa);
	if (ja != nullptr && owner) FREE(ja);
	if (ia != nullptr && owner) FREE(ia);
	aa = rhs.aa; rhs.aa = nullptr;
	ja = rhs.ja; rhs.ja = nullptr;
	ia = rhs.ia; rhs.ia = nullptr;
	matr = rhs.matr; rhs.matr = nullptr;
	descr = rhs.descr;
	return *this;
}

template <class T>
SparseMatr<T> & SparseMatr<T>::operator=(AMatrix<T> &&rhs) {
	SparseMatr<T> && rhs_ = dynamic_cast<SparseMatr<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void SparseMatr<T>::fill(const T val) {
	if (pardisoFact)
		pardisoRelease();
	garbage = false; LUFactorized = false;
	for (Int ind = 0; ind < nz; ind++)
		aa[ind] = val;
}

template <class T>
void SparseMatr<T>::resize(const Int nnew) {
	assert(false);
}

template <class T>
void SparseMatr<T>::solve(Vector<T> & rhssol) {
	assert(n == m);
	assert(rhssol.size() == n);
	if (LUFactorized) {
		solveLU(rhssol);
		return;
	}
	//not LUFactorized, use pardiso
	solvePardiso(rhssol);
}

template <class T>
void SparseMatr<T>::print() const {
	Int ind;
	for (Int i = 0; i < n; i++) {
		ind = ia[i];
		for (Int j = 0; j < m; j++) {
			if (j == ja[ind] && ind != ia[i + 1]) {
				cout << double135<T> << aa[ind] << "  ";
				ind++;
				//if (ind == ia[i+1]) break;
			}
			else
				cout << double135<T> << T() << "  ";
		}
		cout << endl;
	}
}

template <class T>
void SparseMatr<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	f << n << " " << nz << endl;
	Int ind;
	for (Int i = 0; i < n; i++) {
		ind = ia[i];
		for (Int j = 0; j < m; j++) {
			if (j == ja[ind] && ind != ia[i + 1]) {
				f << double2315<T> << aa[ind] << " ";
				ind++;
				//if (ind == ia[i+1]) break;
			}
			else
				f << double2315<T> << T() << " ";
		}
		f << endl;
	}
	f.close();
}

template <class T>
SparseMatr<T> & SparseMatr<T>::operator*=(const T c) {
	assert(!garbage);
	for (Int ind = 0; ind < nz; ind++)
		aa[ind] *= c;
	return *this;
}

template <>
Vector<double> SparseMatr<double>::operator*(\
	const Vector<double> &vec)  const {
	assert(!garbage);
	Int size = vec.size();
	assert(size <= m);
	
	Vector<double> y(n);
	sparse_status_t info = mkl_sparse_d_mv(\
		SPARSE_OPERATION_NON_TRANSPOSE, \
		1.0, matr, descr, vec.v, 0.0, y.v);
	assert(info == SPARSE_STATUS_SUCCESS);

	return y;
}

template <>
Vector<Complex> SparseMatr<Complex>::operator*(\
	const Vector<Complex> &vec)  const {
	assert(!garbage);
	Int size = vec.size();
	assert(size <= m);

	Vector<Complex> y(n);
	sparse_status_t info = mkl_sparse_z_mv(\
		SPARSE_OPERATION_NON_TRANSPOSE, \
		zone, matr, descr, vec.v, zzero, y.v);
	assert(info == SPARSE_STATUS_SUCCESS);

	return y;
}

template <class T>
void SparseMatr<T>::order() {
	sparse_status_t info = \
		mkl_sparse_order(matr);
}

template <class T>
void SparseMatr<T>::fillRandom() {
	assert(n == m); //TODO check if works with n != m
	assert(owner);
	garbage = false;
	RndGenRealUni<T> rgru(-1.0, 1.0);
	rgru.fillArr(aa, nz);
	ia[0] = 0;
	//cout << ia[0] << " ";
	//generate random steps ia[k]-ia[k-1]
	//that are equal to {1,...,n} and sum up to nz
	for (Int i = 1; i <= n; i++)
		ia[i] = 1;
	Int sum = n; Int ind;
	RndGenIntUni rgiu_i(1, n);
	while (sum < nz) {
		ind = rgiu_i.generate();
		if (ia[ind] < n) {
			ia[ind]++;
			sum++;
		}
	}
	//add steps
	for (Int i = 1; i <= n; i++) {
		ia[i] = ia[i - 1] + ia[i];
		//cout << ia[i] << " ";
	}
	//cout << endl;
	RndGenIntUni rgiu_j(0, n - 1); //TODO change (n, m)
	ind = 0; Int ja_new; Int ind_g, ind1, ind2;
	for (Int i = 0; i < n; i++) {
		ja[ind++] = i;
		while (ind < ia[i + 1]) {
			ja_new = rgiu_j.generate();
			//find index of first element greater than ja_new (binary search)
			ind1 = ia[i]; ind2 = ind;
			if (ja_new < ja[ind1]) {
				ind_g = ind1;
			}
			else {
				while (ind1 != ind2 - 1) {
					ind_g = (ind1 + ind2) / 2;
					if (ja[ind_g] > ja_new)
						ind2 = ind_g;
					else
						ind1 = ind_g;
					// ja[ind1] <= ja_new < ja[ind2] now
				}
				ind_g = ind2;
			}
			//insert ja_new if no such element
			if (ja[ind_g - 1] != ja_new) {
				for (Int ind_ = ind; ind_ > ind_g; ind_--)
					ja[ind_] = ja[ind_ - 1];
				ja[ind_g] = ja_new;
				ind++;
			}
		}
	}
	//for (Int ind = 0; ind < nz; ind++)
	//	cout << ja[ind] << " ";
	//cout << endl;
	//MAKE DIAGONALLY DOMINANT
	Int ind_diag;
	double r;
	for (Int i = 0; i < n; i++) {
		ind1 = ia[i]; ind2 = ia[i + 1];
		r = 0.0;
		for (Int ind = ind1; ind < ind2; ind++) {
			r += abs(aa[ind]);
			if (ja[ind] == i)
				ind_diag = ind;
		}
		r *= 1.05;
		aa[ind_diag] = r;
	}
	//for (Int ind = 0; ind < nz; ind++)
	//	cout << aa[ind] << " ";
	//cout << endl;
}

template <class T>
double SparseMatr<T>::sizeGb() const {
	assert(!pardisoFact);
	double res = nz * sizeof(Complex);
	if (owner)
		res += (nz + n + 1) * sizeof(Int);
	return res * bytes2Gbytes;
}

template <>
void SparseMatr<double>::createHandleAndDesc() {
	sparse_status_t info = \
		mkl_sparse_d_create_csr(&matr, SPARSE_INDEX_BASE_ZERO, \
			n, m, ia, ia + 1, ja, aa);
	assert(info == SPARSE_STATUS_SUCCESS);
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
}

template <>
void SparseMatr<Complex>::createHandleAndDesc() {
	sparse_status_t info = \
		mkl_sparse_z_create_csr(&matr, SPARSE_INDEX_BASE_ZERO, \
			n, m, ia, ia + 1, ja, aa);
	assert(info == SPARSE_STATUS_SUCCESS);
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
}

template <>
void SparseMatr<double>::pardisoSetParams() {
	mtype = 11; maxfct = 1; mnum = 1;
	msglvl = 1; //print info
	pardisoinit(pt, &mtype, iparm);
	iparm[0] = 1; //no default
	iparm[1] = 3; //openMP version of fill in analysis
	iparm[5] = 1; //place solution in rhs
	iparm[10] = 0;
	iparm[12] = 0;
	iparm[23] = 1; //openMP factorization
	iparm[24] = 1; //parallel solve
	iparm[34] = 1; //zero-based indexing in ja, ia and perm
}

template <>
void SparseMatr<Complex>::pardisoSetParams() {
	mtype = 13; maxfct = 1; mnum = 1;
	msglvl = 1; //print info
	pardisoinit(pt, &mtype, iparm);
	iparm[0] = 1; //no default
	iparm[1] = 3; //openMP version of fill in analysis
	iparm[5] = 1; //place solution in rhs
	iparm[10] = 0;
	iparm[12] = 0;
	iparm[23] = 1; //openMP factorization
	iparm[24] = 1; //parallel solve
	iparm[34] = 1; //zero-based indexing in ja, ia and perm
}

template <class T>
inline void SparseMatr<T>::solvePardiso(Vector<T> &rhssol) {
	Int phase, error;
	Int nrhs = 1;
	Int *perm = nullptr; //not used ???
	Vector<T> x_tmp(0);
	if (!pardisoFact) { //initialize and factorize
		assert(owner);
		pardisoSetParams();
		phase = 12;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, \
			aa, ia, ja, perm, &nrhs, iparm, &msglvl, \
			rhssol.v, x_tmp.v, &error);
		if (error != 0)
			cout << "Pardiso initialize and factorize error = " << \
				error << endl;
		assert(error == 0);
		pardisoFact = true;
		if (aa != nullptr) FREE(aa);
		garbage = true;
	}
	//solve
	phase = 33;
	x_tmp.resize(rhssol.size());
	//documentation: x_tmp is always used ???
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, \
		aa, ia, ja, perm, &nrhs, iparm, &msglvl, \
		rhssol.v, x_tmp.v, &error);
	x_tmp.resize(0);
	if (error != 0)
		cout << "Pardiso solve error = " << \
			error << endl;
	assert(error == 0);
	cout << endl;
}

template <class T>
void SparseMatr<T>::pardisoRelease() {
	if (!pardisoFact)
		return;
	pardisoFact = false;
	garbage = false;
	Int phase = -1;
	Int error;
	Int nrhs = 1;
	Int *perm = nullptr; //not used ???
	Vector<T> x_tmp(0);
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, \
		aa, ia, ja, perm, &nrhs, iparm, &msglvl, \
		x_tmp.v, x_tmp.v, &error);
	assert(error == 0);
}

template <>
inline void SparseMatr<double>::solveLU(Vector<double> &rhssol) {
	Vector<double> tmp(rhssol.size());
	descr.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
	descr.mode = SPARSE_FILL_MODE_LOWER;
	descr.diag = SPARSE_DIAG_UNIT;
	sparse_status_t info = \
		mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, \
			1.0, matr, descr, rhssol.v, tmp.v);
	assert(info == SPARSE_STATUS_SUCCESS);
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;
	info = \
		mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, \
			1.0, matr, descr, tmp.v, rhssol.v);
	assert(info == SPARSE_STATUS_SUCCESS);
}

template <>
inline void SparseMatr<Complex>::solveLU(Vector<Complex> &rhssol) {
	Vector<Complex> tmp(rhssol.size());
	descr.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
	descr.mode = SPARSE_FILL_MODE_LOWER;
	descr.diag = SPARSE_DIAG_UNIT;
	sparse_status_t info = \
		mkl_sparse_z_trsv(SPARSE_OPERATION_NON_TRANSPOSE, \
			zone, matr, descr, rhssol.v, tmp.v);
	assert(info == SPARSE_STATUS_SUCCESS);
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;
	info = \
		mkl_sparse_z_trsv(SPARSE_OPERATION_NON_TRANSPOSE, \
			zone, matr, descr, tmp.v, rhssol.v);
	assert(info == SPARSE_STATUS_SUCCESS);
}

template <class T>
SparseMatr<T>::~SparseMatr() {
	if (matr != nullptr) {
		sparse_status_t info = \
			mkl_sparse_destroy(matr);
		assert(info == SPARSE_STATUS_SUCCESS);
	}
	pardisoRelease();
	if (aa != nullptr) FREE(aa);
	if (owner) {
		if (ja != nullptr) FREE(ja);
		if (ia != nullptr) FREE(ia);
	}
}
