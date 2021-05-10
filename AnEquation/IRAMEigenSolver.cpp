#include "IRAMEigenSolver.h"
#include "MatrixOp.h"
#include "Algorithms.h"

template class IRAMEigenSolver<double>;
template class IRAMEigenSolver<Complex>;

template <class T>
IRAMEigenSolver<T>::IRAMEigenSolver() \
	: ks(nullptr), ritz(0), resids(0), \
		eps23(pow(eps, 2.0/3)), epsx2(LAPACKE_dlamch('p')) { }

template <>
void IRAMEigenSolver<double>::getHEigAndResid() {
	GenMatrix<double> Hmbar(ks->getHbar());
	double hmp1m = Hmbar.matr[m][m - 1];

	//get Schur decomposition of Hmbar
	GenMatrix<double> z(m); //Schur vectors of Hmbar
	Vector<double> wr(m), wi(m);
	Int info = LAPACKE_dhseqr(LAPACK_ROW_MAJOR, 'S', 'I', \
		m, 1, m, Hmbar.matr[0], m, wr.v, wi.v, z.matr[0], m);
	if (info != 0)
		if (info < 0)
			cout << -info << " th parameter in ?hseqr has illegal value" << endl;
		else //info > 0
			cout << "?hseqr failed to compute all of the eigenvalues" << endl;
	assert(info == 0);
	//Hmbar now contains the upper triangular matrix T: Hmbar=Z*T*Z^H
	//Complex conjugate pairs of ev appear consequtively in wr, wi
	for (Int i = 0; i < m; i++)
		ritz[i] = Complex(wr[i], wi[i]);

	//find eigenvectors of Hmbar
	Int *select = nullptr; double *vl = nullptr; Int m_;
	info = LAPACKE_dtrevc(LAPACK_ROW_MAJOR, 'R', 'B', \
		select, m, Hmbar.matr[0], m, vl, m, z.matr[0], m, m, &m_);
	if (info != 0)
		cout << -info << " th parameter in ?trevc has illegal value" << endl;
	assert(info == 0);
	//if complex pair of ev, consequtive z columns contain
	//real and imaginary parts of eigenvector

	//normalize eigenvectors and obtain residues
	Vector<Complex> evec(m);
	for (Int j = 0; j < m; j++) {
		if (ritz[j].imag() == 0.0) {
			for (Int i = 0; i < m; i++)
				evec[i] = z.matr[i][j];
			resids[j] = abs(evec[m - 1]) / evec.norm();
		} else {
			for (Int i = 0; i < m; i++)
				evec[i] = Complex(z.matr[i][j], z.matr[i][j+1]);
			resids[j] = ABS(evec[m - 1]) / evec.norm();
			resids[++j] = resids[j];
		}
	}
	resids *= hmp1m;
}

template <>
void IRAMEigenSolver<Complex>::getHEigAndResid() {
	GenMatrix<Complex> Hmbar(ks->getHbar());
	double hmp1m = (Hmbar.matr[m][m - 1]).real();

	//get Schur decomposition of Hmbar
	GenMatrix<Complex> z(m); //Schur vectors of Hmbar
	Int info = LAPACKE_zhseqr(LAPACK_ROW_MAJOR, 'S', 'I', \
		m, 1, m, Hmbar.matr[0], m, ritz.v, z.matr[0], m);
	if (info != 0)
		if (info < 0)
			cout << -info << " th parameter in ?hseqr has illegal value" << endl;
		else //info > 0
			cout << "?hseqr failed to compute all of the eigenvalues" << endl;
	assert(info == 0);
	//Hmbar now contains the upper triangular matrix T: Hmbar=Z*T*Z^H

	//find eigenvectors of Hmbar
	Int *select = nullptr; Complex *vl = nullptr; Int m_;
	info = LAPACKE_ztrevc(LAPACK_ROW_MAJOR, 'R', 'B', \
		select, m, Hmbar.matr[0], m, vl, m, z.matr[0], m, m, &m_);
	if (info != 0)
		cout << -info << " th parameter in ?trevc has illegal value" << endl;
	assert(info == 0);

	//normalize eigenvectors and obtain residues
	Vector<Complex> evec(m);
	for (Int j = 0; j < m; j++) {
		for (Int i = 0; i < m; i++)
			evec[i] = z.matr[i][j];
		resids[j] = ABS(evec[m - 1]) / evec.norm();
	}
	resids *= hmp1m;
}

template <>
void IRAMEigenSolver<double>::applyShifts() {
	//get current Hmbar
	GenMatrix<double> & Hmbar = ks->Hbar;
	//make Q identity matrix
	GenMatrix<double> Q(m); Q.fill(0.0);
	for (Int i = 0; i < m; i++)
		Q[i][i] = 1.0;

	//apply q shifts
	bool cmplxTheta = false;
	Int istart, iend, nele, info;
	double h11, h21, c, f, g, s, r, tmp;
	double h12, h22, h32, ab, xu[3], sumth, tau;
	Complex theta;
	double *work = nullptr;
	for (Int jshift = 0; jshift < q; jshift++) {
		theta = ritz[jshift];
		if (theta.imag() != 0.0)
			cmplxTheta = true;
		if ((jshift == q - 1) && cmplxTheta) {
			//last shift is complex without pair => k++, q--, dont apply it
			k++; q--;
			continue;
		}
		istart = 0;
		do { //apply shift jshift to current block

			 //search for next block
			iend = m - 1; //default
			for (Int i = istart; i < m - 1; i++) {
				double sumdiag = abs(Hmbar[i][i]) + \
												abs(Hmbar[i + 1][i + 1]);
				//|Hm_{i+1,i+1}| + |Hm_{i+2,i+2}|
				if (sumdiag == 0.0)
					sumdiag = LAPACKE_dlange(LAPACK_ROW_MAJOR, \
						'1', m - jshift, m - jshift, Hmbar.matr[0], m); //MAGIC! Why m-jshift?
				if (abs(Hmbar[i + 1][i]) <= max(epsx2*sumdiag, small)) {
					//Hm_{i+2,i+1} is zero in fact, the matrix is splitted
					iend = i; //change iend
					Hmbar[i + 1][i] = 0.0;
				}
			} //now next block is defined by istart, iend

			//no sense to apply shift to block of size 1
			if (istart == iend) {
				istart = iend + 1; //start index of next block
				continue;
			}

			//no sense to apply double shift to 2x2 block
			if ((istart + 1 == iend) && (theta.imag() != 0.0)) {
				istart = iend + 1; //start index of next block
				continue;
			}

			//Start implicitly shifted QR
			h11 = Hmbar[istart][istart];
			h21 = Hmbar[istart + 1][istart];
			if (cmplxTheta) { //double shift
				h12 = Hmbar[istart][istart+1];
				h22 = Hmbar[istart+1][istart+1];
				h32 = Hmbar[istart+2][istart+1];
				sumth = 2.0 * theta.real();
				ab = (h11*(h11 - sumth) + norm(theta)) / h21 + h12;
				xu[1] = h11 + h22 - sumth;
				xu[2] = h32;
				for (Int i = istart; i < iend; i++) {
					//only 1 element to zero out in the last but 2 column
					nele = min((Int)3, iend - i + 1);
					//find Householder reflection matrix
					info = LAPACKE_dlarfg(nele, &ab, xu+1, 1, &tau);
					assert(info == 0);
					xu[0] = 1.0;
					if (i > istart) { //return Hm to Hessenberg form
						Hmbar[i][i - 1] = ab;
						Hmbar[i + 1][i - 1] = 0.0;
						if (i < iend - 1) //not last but 2 column
							Hmbar[i + 2][i - 1] = 0.0;
					}
					//sum rows i, i+1 and i+2 of Hm
					//work is not referenced since Householder matrix order is 3 or 2
					//info = LAPACKE_dlarfx(LAPACK_ROW_MAJOR, \
						'L', nele, m-i, xu, tau, Hmbar.matr[i]+i, m, work);
					info = LAPACKE_dlarfx_work(LAPACK_ROW_MAJOR, \
						'L', nele, m - i, xu, tau, Hmbar.matr[i] + i, m, work);
					assert(info == 0);
					//sum columns i, i+1 and i+2 of Hm
					//info = LAPACKE_dlarfx(LAPACK_ROW_MAJOR, \
						'R', min(i+3, iend)+1, nele, xu, tau, Hmbar.matr[0]+i, m, work);
					info = LAPACKE_dlarfx_work(LAPACK_ROW_MAJOR, \
						'R', min(i + 3, iend) + 1, nele, xu, tau, Hmbar.matr[0] + i, m, work);
					assert(info == 0);
					//sum columns i, i+1 and i+2 of Q
					//info = LAPACKE_dlarfx(LAPACK_ROW_MAJOR, \
						'R', m, nele, xu, tau, Q.matr[0]+i, m, work);
					info = LAPACKE_dlarfx_work(LAPACK_ROW_MAJOR, \
						'R', m, nele, xu, tau, Q.matr[0] + i, m, work);
					assert(info == 0);
					//next vector to apply Householder reflection
					if (i < iend - 1) {
						ab = Hmbar[i + 1][i];
						xu[1] = Hmbar[i+2][i];
						if (i < iend - 2)
							xu[2] = Hmbar[i+3][i];
					}
				} //end of double shift
			} else { //single shift
				f = h11 - theta.real(); g = h21;
				for (Int i = istart; i < iend; i++) {
					//find Givens rotation matrix
					cblas_drotg(&f, &g, &c, &s);
					r = f;
					if (i > istart) { //return Hm to Hessenberg form
						//make subdiagonal elements of Hm nonnegative
						if (r < 0.0) {
							r = -r; c = -c; s = -s;
						}
						Hmbar[i][i - 1] = r;
						Hmbar[i + 1][i - 1] = 0.0;
					}
					//sum rows i and i+1 of Hm
					for (Int j = i; j < m; j++) {
						tmp = c * Hmbar[i][j] + s * Hmbar[i + 1][j];
						Hmbar[i + 1][j] = -s * Hmbar[i][j] + c * Hmbar[i + 1][j];
						Hmbar[i][j] = tmp;
					}
					//sum columns j = i and j+1 of Hm
					Int j = i;
					for (Int i_ = 0; i_ < min(i + 3, m); i_++) {
						tmp = c * Hmbar[i_][j] + s * Hmbar[i_][j + 1];
						Hmbar[i_][j + 1] = -s * Hmbar[i_][j] + c * Hmbar[i_][j + 1];
						Hmbar[i_][j] = tmp;
					}
					//sum columns j = i and j+1 of Q
					for (Int i_ = 0; i_ < min(i + 2 + jshift, m); i_++) { //Q adds diagonal after each shift
						tmp = c * Q[i_][j] + s * Q[i_][j + 1];
						Q[i_][j + 1] = -s * Q[i_][j] + c * Q[i_][j + 1];
						Q[i_][j] = tmp;
					}
					//next vector to apply Givens rotation
					if (i < iend - 1) {
						f = Hmbar[i + 1][i];
						g = Hmbar[i + 2][i];
					}
				} //end of single shift
			} //end of implicitly shifted QR
			istart = iend + 1; //start index of next block
		} while (iend < m - 1); //go to next block
		if (cmplxTheta) {
			jshift++;
			cmplxTheta = false;
		}
	} //end applying current shift, go to next

	//perform similarity transformation to make
	//elements of new Hk below diagonal non negative
	for (Int i = 0; i < k; i++) {
		if ( Hmbar[i + 1][i] < 0.0 ) {
			cblas_dscal(m - i, -1.0, Hmbar[i + 1] + i, 1);
			cblas_dscal(min(i + 3, m), -1.0, Hmbar[0] + i + 1, m);
			cblas_dscal(min(i + 2 + q, m), -1.0, Q[0] + i + 1, m);
		}
	}

	//Final check for splitting
	for (Int i = 0; i < k; i++) {
		double sumdiag = abs(Hmbar[i][i]) + abs(Hmbar[i + 1][i + 1]);
		//|Hm_{i+1,i+1}| + |Hm_{i+2,i+2}|
		if (sumdiag == 0.0)
			sumdiag = LAPACKE_dlange(LAPACK_ROW_MAJOR, \
				'1', k, k, Hmbar.matr[0], m);
		if (Hmbar[i + 1][i] <= max(epsx2*sumdiag, small))
			Hmbar[i + 1][i] = 0.0; //Hm_{i+2,i+1} is zero in fact
	}

	//Compute (k+1)-st column of new Vk if nonzero Hm_{k+1, k}
	unique_ptr<AFunction<double>> vkp1;
	if (Hmbar[k][k - 1] > 0.0) {
		vkp1 = (*ks)[0].clone();
		*vkp1 *= Q[0][k];
		for (Int j = 1; j < m; j++) {
			unique_ptr<AFunction<double>> temp = (*ks)[j].clone();
			*temp *= Q[j][k];
			*vkp1 += *temp;
		}
	}

	//Compute first k columns of new Vk, store in last columns of old Vm
	for (Int jc = k - 1; jc >= 0; jc--) {
		*(ks->V[jc + q]) *= Q[jc + q][jc];
		for (Int j = jc + q - 1; j >= 0; j--) {
			unique_ptr<AFunction<double>> temp = (*ks)[j].clone();
			*temp *= Q[j][jc];
			*(ks->V[jc + q]) += *temp;
		}
	}

	//"Copy" new columns in appropriate place
	for (Int jc = 0; jc < k; jc++)
		ks->V[jc] = std::move(ks->V[m - k + jc]);

	//Compute residue of Arnoldi process
	//\hat{v}_{k+1} = Hm_{k+1, k}[Vm*Q]_{:,k+1}+
	//									Q_{m,k}Hmbar_{m+1,m}v_{m+1}
	*(ks->V[m]) *= Hmbar[m][m - 1] * Q[m - 1][k - 1];
	ks->V[k] = std::move(ks->V[m]);
	if (Hmbar[k][k - 1] > 0.0) {
		*vkp1 *= Hmbar[k][k - 1];
		*(ks->V[k]) += *vkp1;
	}

	//Hkbar_{k+1,k} is the \hat{v}_{k+1} norm
	Hmbar[k][k - 1] = (*ks)[k].coef.norm();
	//Compute v_{k+1}
	ks->V[k]->coef.normalize();
}

template <>
void IRAMEigenSolver<Complex>::applyShifts() {
	//get current Hmbar
	GenMatrix<Complex> & Hmbar = ks->Hbar;
	//make Q identity matrix
	GenMatrix<Complex> Q(m); Q.fill(zzero);
	for (Int i = 0; i < m; i++)
		Q[i][i] = zone;

	//apply q shifts
	Int istart, iend;
	Complex theta, h11, h21, f, g, s, r, tmp;
	double c;
	for (Int jshift = 0; jshift < q; jshift++) {
		theta = ritz[jshift];
		istart = 0;
		do { //apply shift jshift to current block

			 //search for next block
			iend = m-1; //default
			for (Int i = istart; i < m - 1; i++) {
				double sumdiag = cblas_dcabs1(&(Hmbar[i][i])) + \
					cblas_dcabs1(&(Hmbar[i+1][i+1]));
				//|Re Hm_{i+1,i+1}|+|Im Hm_{i+1,i+1}|+|Re Hm_{i+2,i+2}|+|Im Hm_{i+2,i+2}|
				if (sumdiag == 0.0)
					sumdiag = LAPACKE_zlange(LAPACK_ROW_MAJOR, \
						'1', m-jshift, m-jshift, Hmbar.matr[0], m); //MAGIC! Why m-jshift?
				if (abs(Hmbar[i + 1][i].real()) <= max(epsx2*sumdiag, small)) {
					//Hm_{i+2,i+1} is zero in fact, the matrix is splitted
					iend = i; //change iend
					Hmbar[i + 1][i] = zzero;
				}
			} //now next block is defined by istart, iend

			//no sense to apply shift to block of size 1
			if (istart == iend) {
				istart = iend + 1; //start index of next block
				continue;
			}

			//Start implicitly shifted QR
			h11 = Hmbar[istart][istart];
			h21 = Hmbar[istart+1][istart];
			f = h11 - theta; g = h21;
			for (Int i = istart; i < iend; i++) {
				//find Givens rotation matrix
				cblas_zrotg(&f, &g, &c, &s);
				r = f;
				if (i > istart) { //return Hm to Hessenberg form
					Hmbar[i][i - 1] = r;
					Hmbar[i + 1][i - 1] = zzero;
				}
				//sum rows i and i+1 of Hm
				for (Int j = i; j < m; j++) {
					tmp = c * Hmbar[i][j] + s * Hmbar[i + 1][j];
					Hmbar[i+1][j] = -CONJ(s) * Hmbar[i][j] + c * Hmbar[i + 1][j];
					Hmbar[i][j] = tmp;
				}
				//sum columns j = i and j+1 of Hm
				Int j = i;
				for (Int i_ = 0; i_ < min(i+3, m); i_++) {
					tmp = c * Hmbar[i_][j] + CONJ(s) * Hmbar[i_][j + 1];
					Hmbar[i_][j+1] = -s * Hmbar[i_][j] + c * Hmbar[i_][j + 1];
					Hmbar[i_][j] = tmp;
				}
				//sum columns j = i and j+1 of Q
				for (Int i_ = 0; i_ < min(i+2+jshift, m); i_++) { //Q adds diagonal after each shift
					tmp = c * Q[i_][j] + CONJ(s) * Q[i_][j + 1];
					Q[i_][j + 1] = -s * Q[i_][j] + c * Q[i_][j + 1];
					Q[i_][j] = tmp;
				}
				//next vector to apply Givens rotation
				if (i < iend - 1) {
					f = Hmbar[i+1][i];
					g = Hmbar[i+2][i];
				}
			} //end of implicitly shifted QR

			istart = iend + 1; //start index of next block
		} while (iend < m-1); //go to next block
	} //end applying current shift, go to next

	//perform similarity transformation to make
	//elements of new Hk below diagonal real and non negative
	for (Int i = 0; i < k; i++) {
		if ((Hmbar[i + 1][i].real() < 0.0) || (Hmbar[i + 1][i].imag() != 0.0)) {
			tmp = Hmbar[i + 1][i] / ABS(Hmbar[i+1][i]);
			cblas_zscal(min(i+3, m), &tmp, Hmbar[0]+i+1, m);
			cblas_zscal(min(i+2+q, m), &tmp, Q[0]+i+1, m);
			tmp = CONJ(tmp);
			cblas_zscal(m-i, &tmp, Hmbar[i+1]+i, 1);
			Hmbar[i + 1][i] = Complex(Hmbar[i+1][i].real(), 0.0);
		}
	}

	//Final check for splitting
	for (Int i = 0; i < k; i++) {
		double sumdiag = cblas_dcabs1(&(Hmbar[i][i])) + \
			cblas_dcabs1(&(Hmbar[i + 1][i + 1]));
		//|Re Hm_{i+1,i+1}|+|Im Hm_{i+1,i+1}|+|Re Hm_{i+2,i+2}|+|Im Hm_{i+2,i+2}|
		if (sumdiag == 0.0)
			sumdiag = LAPACKE_zlange(LAPACK_ROW_MAJOR, \
				'1', k, k, Hmbar.matr[0], m);
		if (Hmbar[i + 1][i].real() <= max(epsx2*sumdiag, small))
			Hmbar[i + 1][i] = zzero; //Hm_{i+2,i+1} is zero in fact
	}

	//Compute (k+1)-st column of new Vk if nonzero Hm_{k+1, k}
	unique_ptr<AFunction<Complex>> vkp1;
	if (Hmbar[k][k-1].real() > 0.0) {
		vkp1 = (*ks)[0].clone();
		*vkp1 *= Q[0][k];
		for (Int j = 1; j < m; j++) {
			unique_ptr<AFunction<Complex>> temp = (*ks)[j].clone();
			*temp *= Q[j][k];
			*vkp1 += *temp;
		}
	}

	//Compute first k columns of new Vk, store in last columns of old Vm
	for (Int jc = k - 1; jc >= 0; jc--) {
		*(ks->V[jc+q]) *= Q[jc + q][jc];
		for (Int j = jc +q - 1; j >= 0; j--) {
			unique_ptr<AFunction<Complex>> temp = (*ks)[j].clone();
			*temp *= Q[j][jc];
			*(ks->V[jc + q]) += *temp;
		}
	}

	//"Copy" new columns in appropriate place
	for (Int jc = 0; jc < k; jc++)
		ks->V[jc] = std::move( ks->V[m-k+jc] );

	//Compute residue of Arnoldi process
	//\hat{v}_{k+1} = Hm_{k+1, k}[Vm*Q]_{:,k+1}+
	//									Q_{m,k}Hmbar_{m+1,m}v_{m+1}
	*(ks->V[m]) *= Hmbar[m][m - 1] * Q[m - 1][k - 1];
	ks->V[k] = std::move(ks->V[m]);
	if (Hmbar[k][k - 1].real() > 0.0) {
		*vkp1 *= Hmbar[k][k - 1];
		*(ks->V[k]) += *vkp1;
	}

	//Hkbar_{k+1,k} is the \hat{v}_{k+1} norm
	Hmbar[k][k - 1] = (*ks)[k].coef.norm();
	//Compute v_{k+1}
	ks->V[k]->coef.normalize();
}

template <>
void IRAMEigenSolver<double>::getEEV(AnEigenProblem<double> &ep, \
	const Int nev, const probType &task, const bool computeVec) {
	header("Start IR Arnoldi process ...");
	AnOperator<double> &op = *(ep.op);
	Vector<Complex> &eval = ep.eval;
	std::vector<shared_ptr<AFunction<double>>> &evec = ep.evec;
	const Int n = op.getRank();
	assert(n > 0); assert((0 < nev) && (nev <= n));
	if (computeVec)
		assert(false);
	assert(!evec.empty());
	AFunction<double> &initv = *(evec[0]);
	assert(initv.size() == n);

	this->task = task;
	small = LAPACKE_dlamch('s') * (n / epsx2);
	m = min(n, (Int)KRYLOV_DIM_MAX_IRAM - 1);
	//m = min(m, 2 * nev+1);
	assert(m >= nev + 2);
	k = nev; //initial numbers of "wanted" eigenvalues ...
	q = m - k; //... and shifts per iteration, can change dynamically
	Int iter = 0; //iteration number

	//generate initial vector
	Int info;
	if (initv.coef.norm2() == 0.0) {
		Int iseed[4] = { 1, 3, 5, 7 };
		info = LAPACKE_dlarnv(2, iseed, n, initv.coef.v);
		assert(info == 0);
	}
	//force initial vector to be in the range of op
	unique_ptr<AFunction<double>> v0 = initv.clone();
	cout << "Step 0 of " << m << " IR Arnoldi" << endl;
	op.times(initv, *v0);
	//create Krylov subspace entity
	ks = new KrylovSubspace<double>(op, *v0, m + 1);
	v0.reset();

	//perform k initial steps of Arnoldi algorithm, so that dim = k + 1
	for (Int j = 0; j < k; j++) {
		cout << "Step " << j + 1 << \
			" of " << m << " IR Arnoldi" << endl;
		ks->makeGramSchmidtStep();
	}

	//main iteration loop
	ritz.resize(m);
	resids.resize(m);
	while (true) {
		iter++;
		//compute q additional steps of Arnoldi algorithm, so that dim = m + 1
		for (Int j = 0; j < q; j++) {
			cout << "Step " << k + j + 1 << \
				" of " << m << " IR Arnoldi" << endl;
			ks->makeGramSchmidtStep();
		}

		//restore initial values of q and k
		k = nev;
		q = m - k;

		//compute eigenvalues of the current matrix H and compute residues
		getHEigAndResid();

		//sort: "unwanted" eigenvalues first, "wanted" last
		//reorder residues appropriately
		reorderEig();
		//possibly k++, q-- to exclude complex shift without pair
		if ((ritz[q - 1].real() == ritz[q].real()) && \
				(ritz[q - 1].imag() == -ritz[q].imag())) {
			k++; q--;
		}
		//then sort "unwanted" part to make eigenvalues with largest residue first
		sortSimult(resids.v, resids.v+q, ritz.v, ritz.v+q, [this](Int i, Int j) {return resids[i] > resids[j]; });
		k0 = k; //nev or nev+1

		//get number of converged eigenvalues among k0 "wanted"
		nconv = 0;
		for (Int j = q; j < m; j++) {
			if (resids[j] <= reltol * max(eps23, ABS(ritz[j])))
				nconv++;
		}

		//if zero residue corresponds to "unwanted" Ritz value make q--, k++
		Int q_old = q;
		for (Int j = 0; j < q_old; j++)
			if (resids[j] == 0.0) {
				q--; k++;
			}

		cout << " Iter = " << iter \
			<< " num converged = " << nconv << " Residues: ";
		for (Int j = q_old; j < m; j++)
			cout << double135<double> << resids[j] << "  ";
		cout << endl;

		//check if next iteration is to be done, exit if not
		bool condRestart = !(nconv == k0) && (iter <= ARNOLDI_ITER_MAX_IRAM) && (q > 0);
		if (!condRestart)
			break;

		//if not converged, possibly increase k. MAGIC!
		Int k_old = k;
		k += min(nconv, q / 2);
		if (k == 1) {
			if (m >= 6)
				k = m / 2;
			else if (m > 3)
				k = 2;
		}
		if (k == m - 1)
			k = m - 2; //keep space for possible k++ in applyShifts() below
		q = m - k;
		if (k > k_old) {//k increased, resort eigenvalues and residues
			reorderEig();
			sortSimult(resids.v, resids.v+q, ritz.v, ritz.v+q, [this](Int i, Int j) {return resids[i] > resids[j]; });
		}

		//apply shifts, possibly k++, q--
		applyShifts();

		ks->reduce(k + 1);
	}

	//get converged "wanted" eigenvalues
	getConvWantedEig();
	eval.resize(nconv);
	for (Int i = 0; i < nconv; i++)
		eval[i] = ritz[i];

	resids.resize(0);  ritz.resize(0);
	delete ks;
}

template <>
void IRAMEigenSolver<Complex>::getEEV(AnEigenProblem<Complex> &ep, \
	const Int nev, const probType &task, const bool computeVec) {
	header("Start IR Arnoldi process ...");
	AnOperator<Complex> & op = *(ep.op);
	Vector<Complex> &eval = ep.eval;
	std::vector<shared_ptr<AFunction<Complex>>> &evec = ep.evec;
	const Int n = op.getRank();
	assert(n > 0); assert( (0 < nev) && (nev <= n) );
	if (computeVec)
		assert(false);
	assert(!evec.empty());
	AFunction<Complex> &initv = *(evec[0]);
	assert(initv.size() == n);

	this->task = task;
	small = LAPACKE_dlamch('s') * (n / epsx2);
	m = min(n, (Int)KRYLOV_DIM_MAX_IRAM-1);
	//m = min(m, 2 * nev + 1);
	assert(m >= nev+1);
	k = nev; //initial numbers of "wanted" eigenvalues ...
	q = m - k; //... and shifts per iteration, can change dynamically
	Int iter = 0; //iteration number

	//generate initial vector
	Int info;
	if (initv.coef.norm2() == 0.0) {
		Int iseed[4] = {1, 3, 5, 7};
		info = LAPACKE_zlarnv(2, iseed, n, initv.coef.v);
		assert(info == 0);
	}
	//force initial vector to be in the range of op
	unique_ptr<AFunction<Complex>> v0 = initv.clone();
	cout << "Step 0 of " << m << " IR Arnoldi" << endl;
	op.times(initv, *v0);
	//create Krylov subspace entity
	ks = new KrylovSubspace<Complex>(op, *v0, m+1);
	v0.reset();
	
	//perform k initial steps of Arnoldi algorithm, so that dim = k + 1
	for (Int j = 0; j < k; j++) {
		cout << "Step " << j + 1 << \
			" of " << m << " IR Arnoldi" << endl;
		ks->makeGramSchmidtStep();
	}

	//main iteration loop
	ritz.resize(m);
	resids.resize(m);
	while (true) {
		iter++;
		//compute q additional steps of Arnoldi algorithm, so that dim = m + 1
		for (Int j = 0; j < q; j++) {
			cout << "Step " << k+j+1 << \
				" of " << m << " IR Arnoldi" << endl;
			ks->makeGramSchmidtStep();
		}

		//restore initial values of q and k
		k = nev;
		q = m - k;

		//compute eigenvalues of the current matrix H and compute residues
		getHEigAndResid();

		//sort: "unwanted" eigenvalues first, "wanted" last
		reorderEig();
		//then sort "unwanted" part to make eigenvalues with largest residue first
		sortSimult(resids.v, resids.v+q, ritz.v, ritz.v+q, [this](Int i, Int j) {return resids[i] > resids[j]; });
		//reorder residues appropriately
		k0 = k;

		//get number of converged eigenvalues among "wanted"
		nconv = 0;
		for (Int j = q; j < m; j++) {
			if (resids[j] <= reltol * max(eps23, ABS(ritz[j])))
				nconv++;
		}

		//if zero residue corresponds to "unwanted" Ritz value make q--, k++
		Int q_old = q;
		for (Int j = 0; j < q_old; j++)
			if (resids[j] == 0.0) {
				q--; k++;
			}

		cout << " Iter = " << iter \
			<< " num converged = " << nconv << " Residues: ";
		for (Int j = q_old; j < m; j++)
			cout << double135<double> << resids[j] << "  ";
		cout << endl;

		//check if next iteration is to be done, exit if not
		bool condRestart =  !(nconv == k0) && (iter <= ARNOLDI_ITER_MAX_IRAM) && (q > 0);
		if (!condRestart)
			break;

		//if not converged, possibly increase k. MAGIC!
		Int k_old = k;
		k += min(nconv, q/2);
		if (k == 1) {
			if (m >= 6)
				k = m / 2;
			else if (m > 3)
				k = 2;
		}
		q = m - k;
		if (k > k_old) {//k increased, resort eigenvalues and residues
			reorderEig();
			sortSimult(resids.v, resids.v+q, ritz.v, ritz.v+q, [this](Int i, Int j) {return resids[i] > resids[j]; });
		}

		//apply shifts
		applyShifts();

		ks->reduce(k+1);
	}

	//get converged "wanted" eigenvalues
	getConvWantedEig();
	eval.resize(nconv);
	for (Int i = 0; i < nconv; i++)
		eval[i] = ritz[i];

	resids.resize(0);  ritz.resize(0);
	delete ks;
}

template <class T>
void IRAMEigenSolver<T>::reorderEig() {
	//sort eigenvalues according to "task"
	switch (task) {
		case LarR:
			sortSimult(ritz.v, ritz.v+m, resids.v, resids.v+m, [this](Int i, Int j) {return ritz[i].real() < ritz[j].real(); });
			break;
		case LowR:
			sortSimult(ritz.v, ritz.v + m, resids.v, resids.v + m, [this](Int i, Int j) {return ritz[i].real() > ritz[j].real(); });
			break;
		case LarM:
			sortSimult(ritz.v, ritz.v + m, resids.v, resids.v + m, [this](Int i, Int j) {return ABS(ritz[i]) < ABS(ritz[j]); });
			break;
		case LowM:
			sortSimult(ritz.v, ritz.v + m, resids.v, resids.v + m, [this](Int i, Int j) {return ABS(ritz[i]) > ABS(ritz[j]); });
			break;
		default:
			assert(false);
	}
}

template <class T>
void IRAMEigenSolver<T>::getConvWantedEig() {
	//reverse eigenvalues and residues arrays
	std::reverse(ritz.v, ritz.v+m);
	std::reverse(resids.v, resids.v+m);

	//At first sort "wanted" eigenvalues wrt relative residues
	//to make those with smallest relative residues first
	for (Int i = 0; i < k0; i++)
		resids[i] = resids[i] / max(eps23, ABS(ritz[i]));
	sortSimult(resids.v, resids.v+k0, ritz.v, ritz.v+k0, [this](Int i, Int j) {return resids[i] < resids[j]; });
	for (Int i = 0; i < k0; i++)
		resids[i] = resids[i] * max(eps23, ABS(ritz[i]));

	//Now converged "wanted" eigenvalues are first
	//Sort them wrt "task" to make most "wanted" last
	switch (task) {
		case LarR:
			sortSimult(ritz.v, ritz.v+nconv, resids.v, resids.v+nconv, [this](Int i, Int j) {return ritz[i].real() < ritz[j].real(); });
			break;
		case LowR:
			sortSimult(ritz.v, ritz.v + nconv, resids.v, resids.v + nconv, [this](Int i, Int j) {return ritz[i].real() > ritz[j].real(); });
			break;
		case LarM:
			sortSimult(ritz.v, ritz.v + nconv, resids.v, resids.v + nconv, [this](Int i, Int j) {return ABS(ritz[i]) < ABS(ritz[j]); });
			break;
		case LowM:
			sortSimult(ritz.v, ritz.v + nconv, resids.v, resids.v + nconv, [this](Int i, Int j) {return ABS(ritz[i]) > ABS(ritz[j]); });
			break;
		default:
			assert(false);
	}
}

template <class T>
IRAMEigenSolver<T>::~IRAMEigenSolver() { }
