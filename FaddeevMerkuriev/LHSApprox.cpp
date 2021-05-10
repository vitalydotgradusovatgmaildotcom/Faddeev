#include "LHSApprox.h"
#include "IdentityM.h"
#include "BlockDiagMatr.h"
#include "TriDiagM.h"
#include "PermutationM.h"
#include "MatrixProd.h"
#include "SpecialFunctions.h"

LHSApprox::LHSApprox( \
		const shared_ptr<const CompsDiscr> &cdiscr, \
		const array<shared_ptr<const APotential<Complex>>, 3> &vc,
		const array<shared_ptr<const APotential<Complex>>, 3> &vs,
		const array<const APotential<Complex> *, 3> &vcy, \
		const System3Body &sys, double E) \
		: VectOperator<Complex, 3>(cdiscr), sys(sys), J(sys.J), Mproj(sys.Mproj), tau(sys.tau), \
		Mmin((1-sys.tau)/2), Mmax(sys.J), nM(Mmax-Mmin+1), cdiscr(cdiscr), \
		vcy(vcy), E(E), matr(nullptr), solComps(cdiscr->getSolComps()), nComp(cdiscr->getNComp()) {
	for (Int alpha = 0; alpha < 3; alpha++) {
		this->vs[alpha] = vs[alpha];
		this->vc[alpha] = vc[alpha];
	}
	getDiscrParams();
	makeMatrix();
}

void LHSApprox::times(const AFunction<Complex> &u, \
				AFunction<Complex> &res) {
	res.coef = *matr * u.coef;
}

void LHSApprox::solve(AFunction<Complex> &rhssol) {
	matr->solve(rhssol.coef);
}

void LHSApprox::makeMatrix() {
	std::vector<shared_ptr<AMatrix<Complex>>> blocks;
	blocks.reserve(solComps.size());
	for (Int alpha : solComps)
		blocks.push_back(makeMatrix(alpha));
	matr = new BlockDiagMatr<Complex>(blocks);
	cout << "Matrix size: " << fixed << \
		sizeGb << " Gbytes" << endl << endl;
//===================================================
}

shared_ptr<AMatrix<Complex>> LHSApprox::makeMatrix(const Int alpha) {
	std::vector<shared_ptr<GenMatrix<Complex>>> wzbar, wz;
	std::vector<shared_ptr<GenMatrix<Complex>>> wxbar, wx;
	std::vector<shared_ptr<DiagMatrix<Complex>>> lambdax;
	std::vector<shared_ptr<GenMatrix<Complex>>> wybar, wy;
	std::vector<shared_ptr<DiagMatrix<Complex>>> lambday;
	
	makeZoperators(alpha, wzbar, wz);
	makeYoperator(alpha, wybar, wy, lambday);
	makeXoperator(alpha, lambday, wxbar, wx, lambdax);

	for (Int k = 0; k < wz.size(); k++)
		sizeGb += wz[k]->sizeGb() + wzbar[k]->sizeGb();
	for (Int k = 0; k < wy.size(); k++)
		sizeGb += wy[k]->sizeGb() + wybar[k]->sizeGb();
	for (Int k = 0; k < wx.size(); k++)
		sizeGb += wx[k]->sizeGb() + wxbar[k]->sizeGb();

	Int n = nM*nx[alpha]*ny[alpha]*nz[alpha];

	std::array<std::shared_ptr<AMatrix<Complex>>, 3> arr3;
	arr3[1] = \
		make_shared<IdentityM<Complex>>(ny[alpha]);
	arr3[2] = \
		make_shared<IdentityM<Complex>>(nx[alpha]);
	std::vector<shared_ptr<AMatrix<Complex>>> wz_, wzbar_;
	wz_.reserve(wz.size()); wzbar_.reserve(wzbar.size());
	for (Int k = 0; k < wz.size(); k++) {
		arr3[0] = wz[k];
		wz_.push_back(make_shared<TensorProd<Complex, 3>>(arr3));
		arr3[0] = wzbar[k];
		wzbar_.push_back(make_shared<TensorProd<Complex, 3>>(arr3));
	}
	shared_ptr<BlockDiagMatr<Complex>> Wz = \
		make_shared<BlockDiagMatr<Complex>>(wz_);
	shared_ptr<BlockDiagMatr<Complex>> Wzbar = \
		make_shared<BlockDiagMatr<Complex>>(wzbar_);

	std::array<std::shared_ptr<AMatrix<Complex>>, 2> arr2;
	arr2[1] = make_shared<IdentityM<Complex>>(nx[alpha]);
	std::vector<shared_ptr<AMatrix<Complex>>> wy_, wybar_;
	wy_.reserve(nM*nz[alpha]); wybar_.reserve(nM*nz[alpha]);
	Int lmin = Mmin;
	for (Int im = Mmin; im <= Mmax; im++)
		for (Int il = im; il < im+nz[alpha]; il++) {
			arr2[0] = wy[il-lmin];
			wy_.push_back(make_shared<TensorProd<Complex, 2>>(arr2));
			arr2[0] = wybar[il-lmin];
			wybar_.push_back(make_shared<TensorProd<Complex, 2>>(arr2));
		}
	shared_ptr<BlockDiagMatr<Complex>> Wy = \
		make_shared<BlockDiagMatr<Complex>>(wy_);
	shared_ptr<BlockDiagMatr<Complex>> Wybar = \
		make_shared<BlockDiagMatr<Complex>>(wybar_);

	std::vector<shared_ptr<AMatrix<Complex>>> wx_, wxbar_;
	wx_.reserve(nM*nz[alpha]*ny[alpha]);
	wxbar_.reserve(nM*nz[alpha]*ny[alpha]);
	Int ind;
	for (Int im = Mmin; im <= Mmax; im++)
		for (Int il = im; il < im+nz[alpha]; il++)
			for (Int iy = 0; iy < ny[alpha]; iy++) {
				ind = (il-lmin)*ny[alpha]+iy;
				wx_.push_back(wx[ind]);
				wxbar_.push_back(wxbar[ind]);
			}
	shared_ptr<BlockDiagMatr<Complex>> Wx = \
		make_shared<BlockDiagMatr<Complex>>(wx_);
	shared_ptr<BlockDiagMatr<Complex>> Wxbar = \
		make_shared<BlockDiagMatr<Complex>>(wxbar_);

	
	//cout << Wz->nrows() << endl;
	//cout << Wy->nrows() << endl;
	//cout << Wx->nrows() << endl;

	//GenMatrix<Complex> Lambda_(n);
	//Lambda_.fill(zzero);
	Int nBlock;
	if (Mmin == Mmax)
		nBlock = 1;
	else
		nBlock = nx[alpha]*ny[alpha]*nz[alpha] - nx[alpha]*ny[alpha];
	// = half band width of non permuted matrix
	std::vector<std::shared_ptr<TriDiagM<Complex>>> blocks;
	blocks.reserve(nBlock); Int ii = 0;
	//cout << "Start bsize check ..." << endl;
	for (Int ib = 0; ib < nBlock; ib++) {
		Int bsize = (n - 1 - ib)/nBlock + 1;
		//cout << bsize << "  ";
		std::shared_ptr<TriDiagM<Complex>> block = \
			make_shared<TriDiagM<Complex>>(bsize);
		block->fill(zzero);
		blocks.push_back(block);
	}
	//cout << endl << "End bsize check ..." << endl;
	//diagonal elements
	Int cntr = 0;
	double jjmm2;
	for (Int im = Mmin; im <= Mmax; im++) {
		jjmm2 = J*(J+1.0)-2.0*im*im;
		for (Int il = im; il < im+nz[alpha]; il++)
			for (Int iy = 0; iy < ny[alpha]; iy++) {
				ind = (il-lmin)*ny[alpha]+iy;
				for (Int ix = 0; ix < nx[alpha]; ix++) {
					//Lambda_[cntr][cntr] = (*lambdax[ind])[ix] + jjmm2;
					blocks[cntr%nBlock]->set(cntr/nBlock, cntr/nBlock, (*lambdax[ind])[ix] + jjmm2);
					cntr++;
				}
			}
	}
	//superdiagonal elements
	cntr = 0;
	Int sh = nx[alpha]*ny[alpha]*nz[alpha] - nx[alpha]*ny[alpha];
	for (Int im = Mmin; im <= Mmax-1; im++) {
		cntr += nx[alpha]*ny[alpha];
		for (Int il = im+1; il < im+nz[alpha]; il++)
			for (Int iy = 0; iy < ny[alpha]; iy++) {
				ind = (il-lmin)*ny[alpha]+iy;
				for (Int ix = 0; ix < nx[alpha]; ix++) {
					//Lambda_[cntr][cntr+sh] = -lambdaJM(il, im)*(lambdaJM(J, im)*sqrt(1.0+DELTA(im, 0))*(1.0-DELTA(im, 0)*DELTA(tau, -1)));
					blocks[cntr%nBlock]->set(cntr/nBlock, cntr/nBlock+1, -lambdaJM(il, im)*(lambdaJM(J, im)*sqrt(1.0+DELTA(im, 0))*(1.0-DELTA(im, 0)*DELTA(tau, -1))));
					cntr++;
				}
			}
	}
	//subdiagonal elements
	cntr = nx[alpha]*ny[alpha]*nz[alpha];
	for (Int im = Mmin+1; im <= Mmax; im++) {
		for (Int il = im; il < im+nz[alpha]-1; il++) {
			for (Int iy = 0; iy < ny[alpha]; iy++) {
				ind = (il-lmin)*ny[alpha]+iy;
				for (Int ix = 0; ix < nx[alpha]; ix++) {
					//Lambda_[cntr][cntr-sh] = lambdaJM(il, -im)*(-lambdaJM(J, -im)*sqrt(1.0+DELTA(im, 1))*(1.0-DELTA(im, 1)*DELTA(tau, -1)));
					blocks[cntr%nBlock]->set(cntr/nBlock, cntr/nBlock-1, lambdaJM(il, -im)*(-lambdaJM(J, -im)*sqrt(1.0+DELTA(im, 1))*(1.0-DELTA(im, 1)*DELTA(tau, -1))));
					cntr++;
				}
			}
		}
		cntr += nx[alpha]*ny[alpha];
	}

	std::vector<std::shared_ptr<AMatrix<Complex>>> blocks_;
	blocks_.reserve(nBlock);
	for (Int ib = 0; ib < nBlock; ib++)
		blocks_.push_back(blocks[ib]);
	std::shared_ptr<BlockDiagMatr<Complex>> bdm = \
		make_shared<BlockDiagMatr<Complex>>(blocks_);

	//bdm->print();
	//Lambda_.print();
	
	std::vector<Int> perm, perm_inv;
	perm.resize(n); perm_inv.resize(n);
	for (Int ib = 0; ib < nBlock; ib++) {
		Int i = ib;
		while (i < n) {
			perm[i] = ii;
			perm_inv[ii++] = i;
			i += nBlock;
		}
	}

	std::vector<shared_ptr<AMatrix<Complex>>> mults;
	mults.resize(3);
	mults[0] = make_shared<PermutationM<Complex>>(std::move(perm_inv));
	mults[1] = bdm;
	mults[2] = make_shared<PermutationM<Complex>>(std::move(perm));
	//mults[2]->write("LAMBDA_PERM.dat");
	//mults[0]->write("LAMBDA_PERM_INV.dat");
	std::shared_ptr<MatrixProd<Complex>> Lambda = \
		make_shared<MatrixProd<Complex>>(mults);

	for (Int ib = 0; ib < blocks.size(); ib++)
		sizeGb += blocks[ib]->sizeGb();
	sizeGb += mults[0]->sizeGb() + mults[2]->sizeGb();
//===========================================
//Проверка

	/*
	GenMatrix<Complex> Lambda_gm(n);
	Vector<Complex> vec(n), vvv(n);

	for (Int i = 0; i < n; i++) {
		if (i % 10 == 0)
			cout << i << endl;
		for (Int j = 0; j < n; j++) {
			vec.fill(zzero);
			vec[j] = zone;
			//vvv = *mults[1] * vec;
			vvv = *Lambda * vec;
			vec.fill(zzero);
			vec[i] = zone;
			Lambda_gm[i][j] = vec.scal(vvv);
		}
	}

	Lambda_gm.write("LAMBDA_3diag.dat");
	*/
	/*for (Int i = 0; i < n; i++)
		for (Int j = 0; j < n; j++) {
			vec.fill(zzero);
			vec[j] = zone;
			vvv = *Lambda*vec;
			vec.fill(zzero);
			vec[i] = zone;
			Lambda_gm[i][j] = vec.scal(vvv);
		}*/

	/*ofstream f;
	f.open("LAMBDA_.dat");
	for (Int i = 0; i < Lambda_.nrows(); i++) {
		for (Int j = 0; j < Lambda_.ncols(); j++) {
			f << Lambda_[i][j] << " ";
			//if (Lambda_[i][j] != zzero)
			//	f << 1;
			//else
			//	f << 0;
		}
		f << endl;
	}
	f.close();

	f.open("LAMBDA.dat");
	for (Int i = 0; i < Lambda_gm.nrows(); i++) {
		for (Int j = 0; j < Lambda_gm.ncols(); j++) {
			f << Lambda_gm[i][j] << " ";
			//if (Lambda_gm[i][j] != zzero)
			//	f << 1;
			//else
			//	f << 0;
		}
		f << endl;
	}
	f.close();*/

	/*double norma = 0.0;
	for (Int i = 0; i < n; i++)
		for (Int j = 0; j < n; j++)
			norma += abs(Lambda_gm[i][j]-Lambda_[i][j]);
	cout << "!!! NORMA = " << norma << endl;*/

//===========================================
	//The final matrix is
	// Wzbar^(-1)*Wybar^(-1)*Wxbar^(-1)*Lambda*Wx^(-1)*Wy^(-1)*Wz^(-1);
	std::vector<std::shared_ptr<AMatrix<Complex>>> res;
	std::vector<bool> inverse;
	res.resize(7); inverse.resize(7);
	res[0] = Wzbar; inverse[0] = true;
	res[1] = Wybar; inverse[1] = true;
	res[2] = Wxbar; inverse[2] = true;
	res[3] = Lambda; inverse[3] = false;
	res[4] = Wx; inverse[4] = true;
	res[5] = Wy; inverse[5] = true;
	res[6] = Wz; inverse[6] = true;

	return make_shared<MatrixProd<Complex>>(res, inverse);
}

void LHSApprox::makeZoperators(const Int alpha, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &wbar, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &w) {
//discretizes z operators and decomposes matrices

	wbar.resize(nM); w.resize(nM);

	//NB! dm1 and dp1 are needed only for tests!
	GenMatrix<Complex> s(0), d(0), dm1(0), dp1(0), dp1_prev(0);
	double zi; Int j;
	Complex fi, dfi;
	NumbersList limits;
	NumbersList::const_iterator c;
	std::vector<Int> inds; inds.resize(nz[alpha]);
	for (Int im = Mmin; im <= Mmax; im++) {
		s.resize(nz[alpha]); d.resize(nz[alpha]);
		dm1.resize(nz[alpha]); dp1.resize(nz[alpha]);
		s.fill(zzero); d.fill(zzero);
		dm1.fill(zzero); dp1.fill(zzero);
		for (Int i = 0; i < nz[alpha]; i++) {
			zi = zcGrid[alpha][i];
			//cout << double135<double> << zi << " ";
			zbasis[alpha]->getNonzero(zi, limits);
			for ( c = limits.begin(); c != limits.end(); c++ ) {
				j = *c;
				fi = zbasis[alpha]->f(zi, j);
				dfi = zbasis[alpha]->d(zi, j);
				s[i][j] = fi;
				/*d[i][j] = -(1.0-zi*zi)*zbasis[alpha]->dd(zi, j) \
					+ 2.0*zi*dfi \
					+ fi*im*im/(1-zi*zi);
				dm1[i][j] = sqrt(1.0-zi*zi)*dfi \
					+ (im-1.0)*zi*fi/sqrt(1.0-zi*zi);
				dp1[i][j] = -sqrt(1.0-zi*zi)*dfi \
					+ (im+1.0)*zi*fi/sqrt(1.0-zi*zi);*/
				d[i][j] = -(1.0 - zi * zi)*zbasis[alpha]->dd(zi, j) \
					+ 2.0*(im+1)*zi*dfi + im*(im+1.0)*fi;
				dm1[i][j] = dfi;
				dp1[i][j] = -(1.0 - zi * zi)*dfi \
					+ 2.0*(im + 1.0)*zi*fi;
			}
		}
		dm1 *= lambdaJM(J, -im)*sqrt(1.0+DELTA(im, 1))*(1.0-DELTA(im, 1)*DELTA(tau, -1));
		dp1 *= lambdaJM(J, im)*sqrt(1.0+DELTA(im, 0))*(1.0-DELTA(im, 0)*DELTA(tau, -1));
		GenMatrix<Complex> d_copy(d);
		GenMatrix<Complex> s_copy(s);

		wbar[im-Mmin] = make_shared<GenMatrix<Complex>>(0);
		w[im-Mmin] = make_shared<GenMatrix<Complex>>(0);
		shared_ptr<DiagMatrix<Complex>> lambda \
			= make_shared<DiagMatrix<Complex>>(0);

		simultDiagonalize(d, s, *wbar[im-Mmin], *lambda, *w[im-Mmin]);

		//sort eigenvalues by absolute value
		std::iota(inds.begin(), inds.end(), 0);
		std::sort(inds.begin(), inds.end(), [&](Int i, Int j) \
			{return abs((*lambda)[i]) < abs((*lambda)[j]);});
		Int icurr, inew;
		for (Int i = 0; i < nz[alpha]; i++) {
			if (inds[i] == -1) continue;
			icurr = i;
			while (inds[icurr] != i) {
				SWAP((*lambda)[icurr], (*lambda)[inds[icurr]]);
				for (Int j = 0; j < nz[alpha]; j++) {
					SWAP((*wbar[im-Mmin])[icurr][j], (*wbar[im-Mmin])[inds[icurr]][j]);
					SWAP((*w[im-Mmin])[j][icurr], (*w[im-Mmin])[j][inds[icurr]]);
				}
				inew = inds[icurr];
				inds[icurr] = -1;
				icurr = inew;
			}
			inds[icurr] = -1;
		}

		//normalize eigenvectors
		//so that \Int_{-1}^{1} dz (1-z^2)^m ( g_m^l(z) )^2 = 1
		//and fix phase using
		//sign( g_m^l(0) ) = (-1)^k if m+l = 2k
		//sign( g'_m^l(0) ) = (-1)^k if m+l = 2k+1
		//GenMatrix<double> ovlpM(zbasis[alpha]->getOverlap());
		GenMatrix<Complex> ovlpM(\
			zbasis[alpha]->getOverlapW(\
				[&im](double z) {return pow(1.0 - z * z, im); }, 2 * im));
		Vector<Complex> coefs(nz[alpha]);
		Complex mult, phase;
		for (Int j = 0; j < nz[alpha]; j++) {
			for (Int i = 0; i < nz[alpha]; i++)
				coefs[i] = (*w[im-Mmin])[i][j];
			mult = zzero;
			for (Int ifu = 0; ifu < nz[alpha]; ifu++)
				for (Int jfu = 0; jfu < nz[alpha]; jfu++)
					mult += conj(coefs[ifu])*ovlpM[ifu][jfu]*coefs[jfu];
			mult = zone/sqrt( mult );
			
			phase = zzero;
			zbasis[alpha]->getNonzero(0.0, limits);
			Int k = (2*im+j)/2; //m+l = 2k or 2k+1
			if ( (2*im+j)%2 == 0 ) //m+l = 2k
				for ( c = limits.begin(); c != limits.end(); c++ )
					phase += coefs[*c]*zbasis[alpha]->f(0.0, *c);
			else //m+l = 2k+1
				for ( c = limits.begin(); c != limits.end(); c++ )
					phase += coefs[*c]*zbasis[alpha]->d(0.0, *c);
			assert(abs(phase) != 0.0);
			phase = pow(-1.0, k)*phase/abs(phase);

			for (Int i = 0; i < nz[alpha]; i++) {
				(*w[im-Mmin])[i][j] *= mult/phase;
				(*wbar[im-Mmin])[j][i] /= mult/phase;
			}
		}
//=============================================================
#ifdef	FILTER_EIGS_LHSAPPROX
		shared_ptr<GenMatrix<Complex>> lam_tilde;
		lam_tilde = make_shared<GenMatrix<Complex>>( \
			(*wbar[im - Mmin])*d_copy*(*w[im - Mmin]));
		double llp1, relerr;
		bool cond;
		for (Int l = max((Int)1, im); l < im + nz[alpha]; l++) {
			llp1 =  l * (l + 1.0);
			relerr = abs( ((*lam_tilde)[l - im][l - im] - llp1) / llp1 );
			cond = relerr > FILTER_EIGS_THRESH;
			if (cond) {
				for (Int l_ = l; l_ < im+nz[alpha]; l_++)
					for (Int i = 0; i < nz[alpha]; i++) {
						(*w[im - Mmin])[i][l_ - im] = zzero;
						(*wbar[im - Mmin])[l_ - im][i] = zzero;
					}
				break;
			}
		}
#endif
//=============================================================
//Проверка
		/*
		shared_ptr<GenMatrix<Complex>> lambdam1, lambdap1, lambda_;
		
		lambda_ = make_shared<GenMatrix<Complex>>( \
					(*wbar[im-Mmin])*d_copy*(*w[im-Mmin]) );
		if (im != Mmin) {
			lambdam1 = make_shared<GenMatrix<Complex>>( \
						(*wbar[im-Mmin])*dm1*(*w[im-Mmin-1]) );
			lambdap1 = make_shared<GenMatrix<Complex>> ( \
						(*wbar[im-Mmin-1])*dp1_prev*(*w[im-Mmin]) );
		}

		cout << "alpha = " << alpha << ", M' = " << im << endl << endl;
		double avrg;
		Complex val;
		if (im != Mmin) {

			//lambdam1->print();
			cout << "Approx upper diag:" << endl;
			double sum1 = 0.0;
			for (Int i = 1; i < nz[alpha]; i++) {
				sum1 += abs((*lambdam1)[i-1][i]);
				cout << (*lambdam1)[i-1][i] << " ";
			}
			cout << endl;
			avrg = 0.0;
			cout << "Exact upper diag:" << endl;
			for (Int m = 0; m < nz[alpha] - 1; m++) {
				val= lambdaJM(im + m, -im)*(-lambdaJM(J, -im)*sqrt(1.0 + DELTA(im, 1))*(1.0 - DELTA(im, 1)*DELTA(tau, -1)));
				cout << val << " ";
				if (abs(val)!= 0.0)
					avrg += abs( ((*lambdam1)[m][m + 1] - val) / val );
			}
			cout << endl;
			avrg /= nz[alpha] - 1;

			double sum = 0.0;
			for (Int i = 0; i < nz[alpha]; i++)
				for (Int j = 0; j < nz[alpha]; j++)
					if (i != j-1)
						sum += abs((*lambdam1)[i][j]);

			//cout << "LambdaUp: " << sum/sum1 << endl;
			cout << "LambdaUp: " << Complex( (sum/(nz[alpha]*nz[alpha] - (nz[alpha]-1)))/(sum1/(nz[alpha]-1)), 0.0) << endl;
			cout << "Average relerr = " << avrg << endl;
			cout << endl;

			//lambdap1->print();
			cout << "Approx lower diag:" << endl;
			sum1 = 0.0;
			for (Int i = 1; i < nz[alpha]; i++) {
				sum1 += abs((*lambdap1)[i][i-1]);
				cout << (*lambdap1)[i][i-1];
			}
			cout << endl;
			avrg = 0.0;
			cout << "Exact lower diag:" << endl;
			for (Int m = 1; m < nz[alpha]; m++) {
				val = Complex(-lambdaJM(im - 1 + m, im - 1)*(lambdaJM(J, im - 1)*sqrt(1.0 + DELTA(im - 1, 0))*(1.0 - DELTA(im - 1, 0)*DELTA(tau, -1))), 0.0);
				cout << val << " ";
				if (abs(val) != 0.0)
					avrg += abs(((*lambdap1)[m][m - 1] - val) / val);
			}
			cout << endl;
			avrg /= nz[alpha] - 1;
			sum = 0.0;
			for (Int i = 0; i < nz[alpha]; i++)
				for (Int j = 0; j < nz[alpha]; j++)
					if (i != j+1)
						sum += abs((*lambdap1)[i][j]);

			//cout << "LambdaDown: " << sum/sum1 << endl;
			cout << "LambdaDown: " << (sum/(nz[alpha]*nz[alpha] - (nz[alpha]-1)))/(sum1/(nz[alpha]-1)) << endl;
			cout << "Average relerr = " << avrg << endl;
		}
		cout << endl;

		cout << "------------------------" << endl << endl;
		//lambda_->print();
		cout << "Approx eigenvalues:" << endl;
		for (Int m = 0; m < nz[alpha]; m++)
			cout << (*lambda_)[m][m] << " ";
		cout << endl;
		avrg = 0.0;
		cout << "Exact eigenvalues:" << endl;
		for (Int l = im; l < im + nz[alpha]; l++) {
			val = l * (l + 1.0);
			cout << val << " ";
			if (abs(val) != 0.0)
				avrg += abs((val - (*lambda_)[l - im][l - im]) / val);
		}
		cout << endl;
		avrg /= nz[alpha];
		cout << "Average relerr = " << avrg << endl;
		cout << "=================================" << endl;
		cout << endl;
		dp1_prev = dp1;
		*/
//==============================================================
//Check preconditioner
	/*
		GenMatrix<Complex> wInv(*w[im - Mmin]), wbarInv(*wbar[im - Mmin]);
		DiagMatrix<Complex> lambda_diag(nz[alpha]);
		wInv.inv();  wbarInv.inv();

		for (Int im_ = im; im_ < im + nz[alpha]; im_++)
			lambda_diag[im_ - im] = im_*(im_ + 1.0);

		(wbarInv * lambda_diag * wInv).write("MATR_DZ_" + std::to_string(alpha) + \
			"_" + std::to_string(im)+ ".dat");

		if (im == Mmin)
			(wbarInv * wInv).write("MATR_SZ_" + std::to_string(alpha) + ".dat");
		
		if (im != Mmin) {
			GenMatrix<Complex> wm1Inv(*w[im - Mmin-1]);
			wm1Inv.inv();
			Int n_ = wbarInv.nrows();
			Int m_ = wm1Inv.nrows();
			GenMatrix<Complex> lambdam1(n_, m_);
			lambdam1.fill(zzero);
			for (Int im_ = im; im_ < im + min(n_, m_ - 1); im_++)
				lambdam1[im_ - im][im_ - im + 1] = lambdaJM(im_, -im)*(-lambdaJM(J, -im)*sqrt(1.0 + DELTA(im, 1))*(1.0 - DELTA(im, 1)*DELTA(tau, -1)));
			//lambdam1 *= zzero;
			//lambdam1.print();
			(wbarInv * lambdam1 * wm1Inv).write("MATR_DZM1_" + std::to_string(alpha) + \
				"_" + std::to_string(im) + ".dat");

			GenMatrix<Complex> wbarm1Inv(*wbar[im - Mmin - 1]);
			wbarm1Inv.inv();
			n_ = wbarm1Inv.nrows();
			m_ = wInv.nrows();
			GenMatrix<Complex> lambdap1(n_, m_);
			lambdap1.fill(zzero);
			for (Int im_ = im + 1; im_ < im + min(n_ - 1, m_) + 1; im_++)
				lambdap1[im_ - im][im_ - im - 1] = -lambdaJM(im_ - 1, im - 1)*(lambdaJM(J, im - 1)*sqrt(1.0 + DELTA(im - 1, 0))*(1.0 - DELTA(im - 1, 0)*DELTA(tau, -1)));
			//lambdap1 *= zzero;
			(wbarm1Inv * lambdap1 * wInv).write("MATR_DZP1_" + std::to_string(alpha) + \
				"_" + std::to_string(im-1) + ".dat");
			//lambdap1.print();
			//cout << "=======" << endl;
			//((*wbar[im - Mmin - 1])*dp1_prev*(*w[im-Mmin])).print();
		}
		//dp1_prev = dp1;
		*/
//==============================================================
	} 
}

void LHSApprox::makeYoperator(const Int alpha, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &wbar, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &w, \
			std::vector<shared_ptr<DiagMatrix<Complex>>> &lambda) {
//discretizes y operator and decomposes matrices

	Int lmin = Mmin; Int lmax = Mmax+nz[alpha]-1;
	wbar.resize(lmax-lmin+1); w.resize(lmax-lmin+1);
	lambda.resize(lmax-lmin+1);

	GenMatrix<Complex> s(0), d(0);
	double yi; Int j;
	Complex fi;
	NumbersList limits;
	NumbersList::const_iterator c;
	for (Int il = lmin; il <= lmax; il++) {
		s.resize(ny[alpha]); d.resize(ny[alpha]);
		s.fill(zzero); d.fill(zzero);
		for (Int i = 0; i < ny[alpha]; i++) {
			yi = ycGrid[alpha][i];
			ybasis[alpha]->getNonzero(yi, limits);
			for (c = limits.begin(); c != limits.end(); c++) {
				j = *c;
				fi = ybasis[alpha]->f(yi, j);
				s[i][j] = fi;
				d[i][j] = -ybasis[alpha]->dd(yi, j) \
					+ fi*( il*(il+1.0)/(yi*yi) + (*vcy[alpha])(yi) );
//====================================
				//Check preconditioner, V(y)
				//d[i][j] = -ybasis[alpha]->dd(yi, j) \
					+ fi * (il*(il + 1.0) / (yi*yi));
//====================================
			}
		}
		wbar[il-lmin] = make_shared<GenMatrix<Complex>>(0);
		w[il-lmin] = make_shared<GenMatrix<Complex>>(0);
		lambda[il-lmin] = make_shared<DiagMatrix<Complex>>(0);

		
		simultDiagonalize(d, s, *wbar[il-lmin], *lambda[il-lmin], *w[il-lmin]);

//=============================================================
//Проверка
		/*
		cout << "alpha = " << alpha << ", l = " << il << endl;
		//if (il == 2) {
			lambda[il-lmin]->print();

			//( (*wbar[il-lmin])*d_copy*(*w[il-lmin]) ).print();
			//cout << "------------------------" << endl;
			//( (*wbar[il-lmin])*s_copy*(*w[il-lmin]) ).print();
		//}
*/
		
//=============================================================
	}
}

void LHSApprox::makeXoperator(const Int alpha, \
			const std::vector<shared_ptr<DiagMatrix<Complex>>> &lambday, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &wbar, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &w, \
			std::vector<shared_ptr<DiagMatrix<Complex>>> &lambda) {
//discretizes x operator and decomposes matrices
	Int lmin = Mmin; Int lmax = Mmax+nz[alpha]-1;
	wbar.resize((lmax-lmin+1)*ny[alpha]);
	w.resize((lmax-lmin+1)*ny[alpha]);
	lambda.resize((lmax-lmin+1)*ny[alpha]);
	
	GenMatrix<Complex> sr2(0), d(0);
	double xi; Int j;
	Complex fi;
	NumbersList limits;
	NumbersList::const_iterator c;
	Int ind;
	for (Int il = lmin; il <= lmax; il++)
		for (Int iy = 0; iy < ny[alpha]; iy++) {
			sr2.resize(nx[alpha]); d.resize(nx[alpha]);
			sr2.fill(zzero); d.fill(zzero);
			for (Int i = 0; i < nx[alpha]; i++) {
				xi = xcGrid[alpha][i];
				xbasis[alpha]->getNonzero(xi, limits);
				for (c = limits.begin(); c != limits.end(); c++) {
					j = *c;
					fi = xbasis[alpha]->f(xi, j);
					sr2[i][j] = fi/(xi*xi);
					d[i][j] = -xbasis[alpha]->dd(xi, j) \
						+ fi*( (*lambday[il-lmin])[iy] \
							+ (*vc[alpha])(xi) + (*vs[alpha])(xi) - E );
				}
			}

			ind = (il-lmin)*ny[alpha]+iy;
			wbar[ind] = make_shared<GenMatrix<Complex>>(0);
			w[ind] = make_shared<GenMatrix<Complex>>(0);
			lambda[ind] = make_shared<DiagMatrix<Complex>>(0);

			simultDiagonalize(d, sr2, *wbar[ind], \
								*lambda[ind], *w[ind]);
			
//=============================================================
//Проверка
			/*
			cout << "alpha = " << alpha << ", l = " << il;
			cout << ", iy = " << iy << ", lambday = " << (*lambday[il-lmin])[iy] << endl;
			if (il == 2 && iy == ny/2) {
				lambda[ind]->print();

				//( (*wbar[ind])*d_copy*(*w[ind]) ).print();
				//cout << "------------------------" << endl;
				//lambda[ind]->print();
				//cout << "------------------------" << endl;
				//( (*wbar[ind])*s_copy*(*w[ind]) ).print();
			}*/
//=============================================================
			for (Int i = 0; i < nx[alpha]; i++)
				(*lambda[ind])[i] +=  il*(il+1.0);
		}
}

double LHSApprox::lambdaJM(const Int j, const Int m) const {
	return sqrt(j*(j+1.0)-m*(m+1.0));
}

void LHSApprox::getDiscrParams() {
	nx.reserve(nComp); ny.reserve(nComp); nz.reserve(nComp);
	xbasis.reserve(nComp); ybasis.reserve(nComp); zbasis.reserve(nComp);
	xcGrid.reserve(nComp); ycGrid.reserve(nComp); zcGrid.reserve(nComp);
	for (Int alpha = 0; alpha < nComp; alpha++) {
		nz.push_back( cdiscr->get(alpha, Mmin).getNi(0) );
		ny.push_back( cdiscr->get(alpha, Mmin).getNi(1) );
		nx.push_back( cdiscr->get(alpha, Mmin).getNi(2) );
		const Collocation<Complex, 3, Complex> &cd = \
			dynamic_cast<const Collocation<Complex, 3, Complex> &> \
			(cdiscr->get(alpha, Mmin));
		zbasis.push_back(cd.bases[0].get());
		ybasis.push_back(cd.bases[1].get());
		xbasis.push_back(cd.bases[2].get());
		zcGrid.push_back(cd.cGrids[0]);
		ycGrid.push_back(cd.cGrids[1]);
		xcGrid.push_back(cd.cGrids[2]);
	}
}

LHSApprox::~LHSApprox(void) {
	if (matr != nullptr)
		delete matr;
}
