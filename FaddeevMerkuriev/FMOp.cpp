#include "FMOp.h"
#include "BlockDiagMatr.h"
#include "BlockMatr.h"
#include "TensorProd.h"
#include "IdentityM.h"
#include "TrivialM.h"
#include "DiagMatrix.h"
#include "MatrixProd.h"

FMOp::FMOp( \
		const shared_ptr<const CompsDiscr> &cdiscr, \
		const array<shared_ptr<const APotential<Complex>>, 3> &vc,
		const array<shared_ptr<const APotential<Complex>>, 3> &vs,
		const array<const MerkurievCutoff<Complex> *, 3> &mc, \
		const System3Body &sys, double energy) \
			: VectOperator<Complex, 3>(cdiscr), cdiscr(cdiscr), sys(sys), \
			E(energy), J(sys.J), Mproj(sys.Mproj), tau(sys.tau), Mmin((1-sys.tau)/2), \
			Mmax(sys.J), nComp(cdiscr->getNComp()), mc(mc), matr(nullptr), \
				matrNStor(nullptr), solComps(cdiscr->getSolComps()), nM(Mmax-Mmin+1), \
					pmz({1.0, 1.0, 1.0}) {

	for (Int alpha = 0; alpha < 3; alpha++) {
		this->vs[alpha] = vs[alpha];
		this->vc[alpha] = vc[alpha];
	}
	for (Int alpha = 0; alpha < 3; alpha++) {
		pm1M[alpha].reserve(nM);
		for (Int im = Mmin; im <= Mmax; im++)
			pm1M[alpha].push_back(1.0);
	}

	getDiscrParams();
	prepareIdentical();
	makeMatrs();
}

Components FMOp::operator*(const Components &u) const {
	//Components res(u.getDiscr());
	Components res(u);
#ifdef	STORE_MATRIX_FMOP
	res.coef = *matr * u.coef;
#else
	res.coef = *matrNStor * u.coef;
#endif
	return res;
}

void FMOp::times(const AFunction<Complex> &u, AFunction<Complex> &res) {
const Components & u_ = \
		dynamic_cast<const Components &>(u);
	Components & res_ = \
		dynamic_cast<Components &>(res);
	res_ = *this * u_;
}

void FMOp::solve(AFunction<Complex> &rhssol) {
#ifdef	STORE_MATRIX_FMOP
	matr->solve(rhssol.coef);
#else
	matrNStor->solve(rhssol.coef);
#endif
}

void FMOp::makeMatrs() {

	//matrices in M space
	std::vector<Complex> m2;
	std::vector<Complex> lamm, lamp;
	GenMatrix<double> f00pm, f00mm;
	std::vector<Complex> mhat, m2hat, jm, lamM;
	m2.reserve(nM);
	f00pm.resize(nM);
	f00mm.resize(nM);
	mhat.reserve(nM); m2hat.reserve(nM);
	jm.reserve(nM); lamM.reserve(nM);
	for (Int im = Mmin; im <= Mmax; im++) {
		m2.push_back(Complex(im*im, 0.0));
		lamm.push_back(Complex(lambdaJM(J, -im)*sqrt(1.0 + DELTA(im, 1))*(1.0 - DELTA(im, 1)*DELTA(tau, -1)), 0.0));
		lamp.push_back(Complex(lambdaJM(J, im)*sqrt(1.0 + DELTA(im, 0))*(1.0 - DELTA(im, 0)*DELTA(tau, -1)), 0.0));
		for (Int imbar = Mmin; imbar <= Mmax; imbar++) {
			f00pm[im - Mmin][imbar - Mmin] = pow(-1.0, imbar - im)*2.0 / (sqrt(2.0 + 2.0*DELTA(imbar, 0))*sqrt(2.0 + 2.0*DELTA(im, 0)));
			f00mm[im - Mmin][imbar - Mmin] = f00pm[im - Mmin][imbar - Mmin] * tau*pow(-1.0, im);
		}
		mhat.push_back(Complex(2.0*(im+1), 0.0));
		m2hat.push_back(Complex(im*(im+1.0), 0.0));
		jm.push_back(Complex(J*(J+1.0)-2.0*im*im, 0.0));
		lamM.push_back(mhat.back()*lamp.back());
	}

	std::vector<GenMatrix<Complex>> sz;
	std::vector<std::vector<GenMatrix<Complex>>> \
		dsz, dszm1, dszp1;
	sz.resize(nComp);
	dsz.resize(nComp);
	dszm1.resize(nComp); dszp1.resize(nComp);
	std::vector<GenMatrix<Complex>> d1z, zd1z, z2d1z, zsz, dsz0;
	d1z.resize(nComp); zd1z.resize(nComp);
	z2d1z.resize(nComp); zsz.resize(nComp);
	dsz0.resize(nComp);
	std::vector<shared_ptr<SparseMatr<Complex>>> \
		sxsy, xr2sxsy, xr2yr2, dsxdsy;
	std::vector<GenMatrix<Complex>> vsxydotxy;
	sxsy.resize(nComp);
	xr2sxsy.resize(nComp);
	xr2yr2.resize(nComp); dsxdsy.resize(nComp);
	vsxydotxy.resize(nComp);
	Int j;
	NumbersList limits;
	NumbersList::const_iterator c;
	double xi, yi, zi_;
	std::array<Int, 3> nnzx = { 0, 0, 0 };
	std::array<Int, 3> nnzy = { 0, 0, 0 };
	std::array<Int, 3> nnzz = { 0, 0, 0 };
	for (Int alpha : solComps) {
		//1d MATRICES FOR LHS OPERATOR (BLOCK ALPHA)
		GenMatrix<Complex> sx(nx[alpha]), sy(ny[alpha]);
		GenMatrix<Complex> xr2sx(nx[alpha]), yr2sy(ny[alpha]);
		GenMatrix<Complex> dsx(nx[alpha]);
		GenMatrix<Complex> dsy(ny[alpha]);
		sz[alpha].resize(nz[alpha]);
		dsz[alpha].reserve(nM);
		dszm1[alpha].reserve(nM);
		dszp1[alpha].reserve(nM);
		for (Int im = Mmin; im <= Mmax; im++) {
			dsz[alpha].push_back(GenMatrix<Complex>(nz[alpha]));
			dszm1[alpha].push_back(GenMatrix<Complex>(nz[alpha]));
			dszp1[alpha].push_back(GenMatrix<Complex>(nz[alpha]));
		}
		d1z[alpha].resize(nz[alpha]); zd1z[alpha].resize(nz[alpha]);
		z2d1z[alpha].resize(nz[alpha]); zsz[alpha].resize(nz[alpha]);
		dsz0[alpha].resize(nz[alpha]);
		
		Complex fi, dfi;
		double xi2, yi2;
		Complex vxi, vsxi;
		for (Int i = 0; i < nx[alpha]; i++) {
			xi = xcGrid[alpha][i];
			xbasis[alpha]->getNonzero(xi, limits);
			vxi = (*vc[alpha])(xi) + (*vs[alpha])(xi);
			xi2 = xi*xi;
			nnzx[alpha] += limits.size();
			for ( c = limits.begin(); c != limits.end(); c++ ) {
				j = *c;
				sx[i][j] = fi = xbasis[alpha]->f(xi, j);
				dsx[i][j] = -xbasis[alpha]->dd(xi, j) + vxi * fi;
				xr2sx[i][j] = fi/xi2;
			}
		}
		for (Int i = 0; i < ny[alpha]; i++) {
			yi = ycGrid[alpha][i];
			ybasis[alpha]->getNonzero(yi, limits);
			yi2 = yi*yi;
			nnzy[alpha] += limits.size();
			for ( c = limits.begin(); c != limits.end(); c++ ) {
				j = *c;
				sy[i][j] = fi = ybasis[alpha]->f(yi, j);
				yr2sy[i][j] = fi/yi2;
				dsy[i][j] = -ybasis[alpha]->dd(yi, j) - E*fi;
			}
		}
		Complex szdiv, sq1mzi2, zidivsq;
		Complex dzsq, szzidiv;
		for (Int i = 0; i < nz[alpha]; i++) {
			zi_ = zcGrid[alpha][i];
			zbasis[alpha]->getNonzero(zi_, limits);
			nnzz[alpha] += limits.size();
			for (c = limits.begin(); c != limits.end(); c++) {
				j = *c;
				sz[alpha][i][j] = fi = zbasis[alpha]->f(zi_, j);
				dfi = zbasis[alpha]->d(zi_, j);
				dsz0[alpha][i][j] = -(1.0 - zi_ * zi_)*zbasis[alpha]->dd(zi_, j);
				for (Int im = Mmin; im <= Mmax; im++) {
					dsz[alpha][im - Mmin][i][j] = dsz0[alpha][i][j] + \
						2.0*(im+1.0)*zi_*dfi + im*(im+1.0)*fi;
					dszm1[alpha][im - Mmin][i][j] = \
						lamm[im - Mmin] * dfi;
					dszp1[alpha][im - Mmin][i][j] = \
						lamp[im - Mmin] * (-(1.0 - zi_ * zi_)*dfi + 2.0*(im + 1.0)*zi_*fi);
				}
				d1z[alpha][i][j] = dfi;
				zd1z[alpha][i][j] = zi_ * dfi;
				z2d1z[alpha][i][j] = (1.0 - zi_ * zi_) * dfi;
				zsz[alpha][i][j] = zi_ * fi;
			}
		}
//============================================
//Check preconditioner
	/*
		ifstream f;
		string filename;
		filename = "MATR_SZ_" + std::to_string(alpha) + ".dat";
		f.open(filename);
		for (Int i = 0; i < sz[alpha].nrows(); i++)
			for (Int j = 0; j < sz[alpha].ncols(); j++)
				f >> sz[alpha][i][j];
		f.close();
		
		for (Int im = Mmin; im <= Mmax; im++) {
			filename = "MATR_DZ_" + std::to_string(alpha) + \
				"_" + std::to_string(im) + ".dat";
			f.open(filename);
			for (Int i = 0; i < dsz[alpha][im - Mmin].nrows(); i++)
				for (Int j = 0; j < dsz[alpha][im - Mmin].ncols(); j++)
					f >> dsz[alpha][im - Mmin][i][j];
			f.close();
			//dsz[alpha][im - Mmin].print();
			if (im > Mmin) {
				filename = "MATR_DZM1_" + std::to_string(alpha) + \
				"_" + std::to_string(im) + ".dat";
				f.open(filename);
				for (Int i = 0; i < dszm1[alpha][im - Mmin].nrows(); i++)
					for (Int j = 0; j < dszm1[alpha][im - Mmin].ncols(); j++)
						f >> dszm1[alpha][im - Mmin][i][j];
				f.close();
				//dszm1[alpha][im - Mmin].print();

				filename = "MATR_DZP1_" + std::to_string(alpha) + \
					"_" + std::to_string(im - 1) + ".dat";
				f.open(filename);
				for (Int i = 0; i < dszp1[alpha][im-1 - Mmin].nrows(); i++)
					for (Int j = 0; j < dszp1[alpha][im-1 - Mmin].ncols(); j++)
						f >> dszp1[alpha][im-1 - Mmin][i][j];
				f.close();
				//dszp1[alpha][im - 1 - Mmin].print();
			}
		}
		nnzz[alpha] = nz[alpha] * nz[alpha];
		*/
//============================================
		//2d MATRICES FOR LHS OPERATOR (BLOCK ALPHA)
		Int n = nx[alpha]*ny[alpha];
		long long int nnz = (long long int)nnzx[alpha]*nnzy[alpha];
		sxsy[alpha] = make_shared<SparseMatr<Complex>>(n, nnz);
		xr2sxsy[alpha] = make_shared<SparseMatr<Complex>>(n, nnz);
		xr2yr2[alpha] = make_shared<SparseMatr<Complex>>(n, nnz);
		dsxdsy[alpha] = make_shared<SparseMatr<Complex>>(n, nnz);
		vsxydotxy[alpha].resize(nx[alpha], ny[alpha]);
		sxsy[alpha]->ia[0] = 0;
		xr2sxsy[alpha]->ia[0] = 0;
		xr2yr2[alpha]->ia[0] = 0;
		dsxdsy[alpha]->ia[0] = 0;
		NumbersList limitsx, limitsy;
		NumbersList::const_iterator cx, cy;
		Int jx, jy; Int i_ = 0; Int cntr = 0;
		for (Int iy = 0; iy < ny[alpha]; iy++) {
			yi = ycGrid[alpha][iy];
			ybasis[alpha]->getNonzero(yi, limitsy);
			for (Int ix = 0; ix < nx[alpha]; ix++) {
				xi = xcGrid[alpha][ix];
				xbasis[alpha]->getNonzero(xi, limitsx);
				vxi = (*vc[alpha])(xi);
				vsxi = (*vs[alpha])(xi);
//====================================
				//Check preconditioner, rhs
				//vxi = zzero;
				//vsxi = zzero;
//====================================
				for ( cy = limitsy.begin(); cy != limitsy.end(); cy++ ) {
					jy = *cy;
					for ( cx = limitsx.begin(); cx != limitsx.end(); cx++ ) {
						jx = *cx;
						j = jy*nx[alpha]+jx;
						sxsy[alpha]->aa[cntr] = sx[ix][jx]*sy[iy][jy];
						xr2sxsy[alpha]->aa[cntr] = xr2sx[ix][jx]*sy[iy][jy];
						xr2yr2[alpha]->aa[cntr] = xr2sxsy[alpha]->aa[cntr] + \
							sx[ix][jx]* yr2sy[iy][jy];
						dsxdsy[alpha]->aa[cntr] = \
							sx[ix][jx]*dsy[iy][jy] + dsx[ix][jx]*sy[iy][jy];
						xr2yr2[alpha]->ja[cntr] = dsxdsy[alpha]->ja[cntr] = \
							sxsy[alpha]->ja[cntr] = xr2sxsy[alpha]->ja[cntr] = j;
						cntr++;
					}
				}
				i_++;
				xr2yr2[alpha]->ia[i_] = dsxdsy[alpha]->ia[i_] = \
					sxsy[alpha]->ia[i_] = xr2sxsy[alpha]->ia[i_] = cntr;
				vsxydotxy[alpha][ix][iy] = \
					( vxi*(*mc[alpha])(xi, yi)+vsxi )*xi*yi;
			}
		}
	}

	//3D MATRICES FOR RHS OPERATOR
	NumbersList limsxrot, limsyrot, limszrot;
	jacobiCoo xyza, xyzb;
	double xroti, yroti, zroti, zroti2;
	NumbersList::const_iterator cx, cy, cz;
	Complex firot;
	std::array<std::vector<std::vector<Complex>>, 3> szrot, sxdivxrot, sydivyrot;
	std::array<std::vector<std::vector<Complex>>, 3> mz2powmhalf;
	for (Int alpha : solComps) {
		sxdivxrot[alpha].reserve((nComp-1)*nx[alpha]*ny[alpha]*nz[alpha]);
		sydivyrot[alpha].reserve((nComp - 1)*nx[alpha] * ny[alpha] * nz[alpha]);
		szrot[alpha].reserve((nComp - 1)*nx[alpha] * ny[alpha] * nz[alpha]);
		mz2powmhalf[alpha].reserve((nComp - 1)*nx[alpha] * ny[alpha] * nz[alpha]);
		for (Int iz = 0; iz < nz[alpha]; iz++) {
			xyza.z = zcGrid[alpha][iz];
			for (Int iy = 0; iy < ny[alpha]; iy++) {
				xyza.y = ycGrid[alpha][iy];
				for (Int ix = 0; ix < nx[alpha]; ix++) {
					xyza.x = xcGrid[alpha][ix];
					for (Int beta = 0; beta < nComp; beta++) {
						if (beta != alpha) {
							sys.rotReduced(xyza, alpha, xyzb, beta);
							xroti = xyzb.x; yroti = xyzb.y;
							//zroti = identical ? abs(xyzb.z) : xyzb.z;
							zroti = xyzb.z;
							zroti *= pmz[beta];
							xbasis[beta]->getNonzero(xroti, limsxrot);
							ybasis[beta]->getNonzero(yroti, limsyrot);
							zbasis[beta]->getNonzero(zroti, limszrot);

							std::vector<Complex> sxdivxrot_, sydivyrot_, szrot_;
							szrot_.reserve(limszrot.size());
							sxdivxrot_.reserve(limsxrot.size()); sydivyrot_.reserve(limsyrot.size());
							for (cz = limszrot.begin(); cz != limszrot.end(); cz++) {
								firot = zbasis[beta]->f(zroti, *cz);
								//if (xyzb.z < 0.0)
								//	firot *= symmetry;
								szrot_.push_back(firot);
							}
							for (cy = limsyrot.begin(); cy != limsyrot.end(); cy++)
								sydivyrot_.push_back(ybasis[beta]->f(yroti, *cy) / yroti);
							for ( cx = limsxrot.begin(); cx != limsxrot.end(); cx++ )
								sxdivxrot_.push_back(xbasis[beta]->f(xroti, *cx)/xroti);
							sxdivxrot[alpha].push_back(sxdivxrot_);
							sydivyrot[alpha].push_back(sydivyrot_);
							szrot[alpha].push_back(szrot_);
							std::vector<Complex> mz2powmhalf_;
							mz2powmhalf_.reserve(nM);
							zroti2 = zroti * zroti;
							for (Int im = Mmin; im <= Mmax; im++)
								mz2powmhalf_.push_back(pow(1.0-zroti2, 0.5*im));
							mz2powmhalf[alpha].push_back(mz2powmhalf_);
						}
					}
				}
			}
		}
	}

#ifdef STORE_MATRIX_FMOP
	
	//CALCULATE TOTAL NUMBER OF NONZERO ELEMENTS
	NumbersList limitsx, limitsy, limitsz;
	Int *ja = new Int[2 * rx*ry*rz];
	Complex *aa = new Complex[2 * rx*ry*rz];
	Int ind;
	long long int nnz_lhs, nnz_rhs, nnz = 0;
	for (Int alpha : solComps) {
		for (Int iz = 0; iz < nz[alpha]; iz++) {
			xyza.z = zcGrid[alpha][iz];
			for (Int iy = 0; iy < ny[alpha]; iy++) {
				xyza.y = ycGrid[alpha][iy];
				for (Int ix = 0; ix < nx[alpha]; ix++) {
					ind = 0;
					xyza.x = xcGrid[alpha][ix];
					//beta != alpha
					for (Int beta = 0; beta < nComp; beta++) {
						if (beta != alpha) {
							sys.rotReduced(xyza, alpha, xyzb, beta);
							xroti = xyzb.x; yroti = xyzb.y;
							//zroti = identical ? abs(xyzb.z) : xyzb.z;
							zroti = xyzb.z;
							zroti *= pmz[beta];
							xbasis[beta]->getNonzero(xroti, limsxrot);
							ybasis[beta]->getNonzero(yroti, limsyrot);
							zbasis[beta]->getNonzero(zroti, limszrot);
							for (cz = limszrot.begin(); cz != limszrot.end(); cz++)
								for (cy = limsyrot.begin(); cy != limsyrot.end(); cy++)
									for (cx = limsxrot.begin(); cx != limsxrot.end(); cx++)
										ja[ind++] = \
										cdiscr->getRaw(beta, Mmin, { *cz, *cy, *cx });
						}
					}
					orderAndSqueeze(ind, aa, ja);
					nnz_rhs = ind;
					//beta = alpha
					xbasis[alpha]->getNonzero(xyza.x, limitsx);
					ybasis[alpha]->getNonzero(xyza.y, limitsy);
					zbasis[alpha]->getNonzero(xyza.z, limitsz);
//============================================
					//Check preconditioner

					//limitsz.resize(nz[alpha]);
					//std::iota(limitsz.begin(), limitsz.end(), 0);

//============================================
					nnz_lhs = 0; //those that are not in rhs
					for (cz = limitsz.begin(); cz != limitsz.end(); cz++)
						for (cy = limitsy.begin(); cy != limitsy.end(); cy++)
							for (cx = limitsx.begin(); cx != limitsx.end(); cx++) {
								j = cdiscr->getRaw(alpha, Mmin, { *cz, *cy, *cx });
								if (!std::binary_search(ja, ja + ind, j))
									nnz_lhs++;
							}
					nnz += nnz_rhs * nM*nM + \
									nnz_lhs * ( 3 * nM - 2 );
				}
			}
		}
	}
	delete[] ja; delete[] aa;

	//COMPUTE MATRIX
	Int jx, jy, jz;
	double inv1mz2pow;
	matr = \
		new SparseMatr<Complex>(cdiscr->getN(), nnz);
	matr->ia[0] = 0;
	long long int cntr_ = 0;
	Int i = 0;
	std::array<Int, 3> js;
	Int cntr;
	double wba;
	Complex dsxyv, f00mv, vlxy, vsxydotxyz;
	ja = new Int[rx*ry*rz*(2*nM+3)];
	aa = new Complex[rx*ry*rz*(2 * nM + 3)];
	for (Int alpha : solComps) {
		for (Int im = Mmin; im <= Mmax; im++) {
			SparseMatr<Complex> dsxdsy_(*xr2sxsy[alpha]);
			dsxdsy_ *= J*(J+1.0)-2.0*im*im;
			for (Int ind = 0; ind < dsxdsy_.nnonz(); ind++)
				dsxdsy_.aa[ind] += dsxdsy[alpha]->aa[ind];
			Int cntr_xyzb = 0;
			for (Int iz = 0; iz < nz[alpha]; iz++) {
				xyza.z = zi_ = zcGrid[alpha][iz];
				zbasis[alpha]->getNonzero(zi_, limitsz);
				inv1mz2pow = 1.0 / pow(1.0 - zi_ * zi_, 0.5*im);
//============================================
//Check preconditioner
				
				//limitsz.resize(nz[alpha]);
				//std::iota(limitsz.begin(), limitsz.end(), 0);
				
//============================================
				Int cntr_xy = 0; Int cntr_xy_old;
				for (Int iy = 0; iy < ny[alpha]; iy++) {
					xyza.y = yi = ycGrid[alpha][iy];
					ybasis[alpha]->getNonzero(yi, limitsy);
					for (Int ix = 0; ix < nx[alpha]; ix++) {
						xyza.x = xi = xcGrid[alpha][ix];
						xbasis[alpha]->getNonzero(xi, limitsx);
						vsxydotxyz = inv1mz2pow * vsxydotxy[alpha][ix][iy];
						cntr = 0;
						for (Int beta = 0; beta < nComp; beta++) {
							if (beta == alpha) { //beta == alpha
								vlxy = zzero;
								for (Int gamma = 0; gamma < nComp; gamma++) {
									sys.rotReduced(xyza, alpha, xyzb, gamma);
									xroti = xyzb.x; yroti = xyzb.y;
									if (gamma != alpha) {
										vlxy += (*vc[gamma])(xroti)*mc[gamma]->tail(xroti, yroti);
									}
								}
								for (Int gamma = nComp; gamma < 3; gamma++) {
									sys.rotReduced(xyza, alpha, xyzb, gamma);
									xroti = xyzb.x;
									vlxy += (*vc[gamma])(xroti) + (*vs[gamma])(xroti);
								}
//====================================
//Check preconditioner, tails
								//vlxy = zzero;
//====================================
								cntr_xy_old = cntr_xy;
								if (im != Mmin) {
									for ( cy = limitsy.begin(); cy != limitsy.end(); cy++ ) {
										js[1] = jy = *cy;
										for ( cx = limitsx.begin(); cx != limitsx.end(); cx++ ) {
											js[2] = jx = *cx;
											for ( cz = limitsz.begin(); cz != limitsz.end(); cz++ ) {
												js[0] = jz = *cz;
												j = cdiscr->getRaw(alpha, im-1, js);
												aa[cntr] = dszm1[alpha][im-Mmin][iz][jz]* xr2sxsy[alpha]->aa[cntr_xy];
												ja[cntr] = j; //NB not sorted!
												cntr++;
											}
											cntr_xy++;
										}
									}
								}
								cntr_xy = cntr_xy_old;
								for ( cy = limitsy.begin(); cy != limitsy.end(); cy++ ) {
									js[1] = jy = *cy;
									for ( cx = limitsx.begin(); cx != limitsx.end(); cx++ ) {
										js[2] = jx = *cx;
										dsxyv = dsxdsy_.aa[cntr_xy] + vlxy*sxsy[alpha]->aa[cntr_xy];
										for ( cz = limitsz.begin(); cz != limitsz.end(); cz++ ) {
											js[0] = jz = *cz;
											j = cdiscr->getRaw(alpha, im, js);
											aa[cntr] = sz[alpha][iz][jz]*dsxyv + \
												dsz[alpha][im-Mmin][iz][jz]*xr2yr2[alpha]->aa[cntr_xy];
											ja[cntr] = j; //NB not sorted!
											cntr++;
										}
										cntr_xy++;
									}
								}
								if (im != Mmax) {
									cntr_xy = cntr_xy_old;
									for ( cy = limitsy.begin(); cy != limitsy.end(); cy++ ) {
										js[1] = jy = *cy;
										for ( cx = limitsx.begin(); cx != limitsx.end(); cx++ ) {
											js[2] = jx = *cx;
											for ( cz = limitsz.begin(); cz != limitsz.end(); cz++ ) {
												js[0] = jz = *cz;
												j = cdiscr->getRaw(alpha, im+1, js);
												aa[cntr] = dszp1[alpha][im-Mmin][iz][jz]* xr2sxsy[alpha]->aa[cntr_xy];
												ja[cntr] = j; //NB not sorted!
												cntr++;
											}
											cntr_xy++;
										}
									}
								}
							} else { //beta != alpha
								sys.rotReduced(xyza, alpha, xyzb, beta);
								xroti = xyzb.x; yroti = xyzb.y;
								//zroti = identical ? abs(xyzb.z) : xyzb.z;
								zroti = xyzb.z;
								zroti *= pmz[beta];
								xbasis[beta]->getNonzero(xroti, limsxrot);
								ybasis[beta]->getNonzero(yroti, limsyrot);
								zbasis[beta]->getNonzero(zroti, limszrot);
								
								std::vector<Complex> sxyzdivxyrot;
								sxyzdivxyrot.reserve(limsxrot.size()*limsyrot.size()*limszrot.size());
								std::vector<Complex>::const_iterator cxyz;
								Complex szydivyrot;
								for ( cz = limszrot.begin(); cz != limszrot.end(); cz++ ) {
									for ( cy = limsyrot.begin(); cy != limsyrot.end(); cy++ ) {
										szydivyrot = szrot[alpha][cntr_xyzb][cz-limszrot.begin()]\
													*sydivyrot[alpha][cntr_xyzb][cy-limsyrot.begin()];
										for ( cx = limsxrot.begin(); cx != limsxrot.end(); cx++ ) {
											sxyzdivxyrot.push_back(szydivyrot*sxdivxrot[alpha][cntr_xyzb][cx-limsxrot.begin()]);
										}
									}
								}
								wba = acos( (sys.s[beta][alpha] * yi*zi_ + sys.c[beta][alpha] * xi) / xroti );
								if ((beta - alpha + 1) % 3 != 0)
									wba = 2 * PI - wba;
								for (Int imbar = Mmin; imbar <= Mmax; imbar++) {
									f00mv = f00pm[im-Mmin][imbar-Mmin]*wignerDSmall(J, imbar, im, wba) + \
												f00mm[im-Mmin][imbar-Mmin]*wignerDSmall(J, imbar, -im, wba);
									f00mv *= pm1M[beta][imbar-Mmin] * \
										vsxydotxyz*mz2powmhalf[alpha][cntr_xyzb][imbar-Mmin];
									cxyz = sxyzdivxyrot.begin();
									for ( cz = limszrot.begin(); cz != limszrot.end(); cz++ ) {
										js[0] = *cz;
										for ( cy = limsyrot.begin(); cy != limsyrot.end(); cy++ ) {
											js[1] = *cy;
											for ( cx = limsxrot.begin(); cx != limsxrot.end(); cx++ ) {
												js[2] = *cx;
												j = cdiscr->getRaw(beta, imbar, js);
												aa[cntr] = f00mv*(*cxyz);
												ja[cntr] = j;
												cntr++;
												cxyz++;
											}
										}
									}
								}
								cntr_xyzb++;
							}
						}
						orderAndSqueeze(cntr, aa, ja);
						std::copy(ja, ja+cntr, matr->ja+cntr_);
						std::copy(aa, aa + cntr, matr->aa + cntr_);
						cntr_ += cntr;
						matr->ia[++i] = cntr_;
					}
				}
			}
		}
	}
	delete[] aa; delete[] ja;
	nnzTot = matr->nnonz();
	sizeGb = matr->sizeGb();
#else
	//assert(Mmin == Mmax);
	std::vector<shared_ptr<AMatrix<Complex>>>	 lhsrhs;
	lhsrhs.reserve(2);

	//LHS TERM
	std::vector<shared_ptr<AMatrix<Complex>>> blocks;
	blocks.reserve(solComps.size());
	shared_ptr<SparseMatr<Complex>> msz, mdz;
	for (Int alpha : solComps) {
		//more 2d MATRICES FOR LHS OPERATOR (BLOCK ALPHA)
		msz = make_shared<SparseMatr<Complex>>( \
								nM*nz[alpha], (long long int)(3*nM-2)*nnzz[alpha]);
		mdz = make_shared<SparseMatr<Complex>>( \
								nM*nz[alpha], (long long int)nM*nnzz[alpha]);
		Int cntr1 = 0, cntr2 = 0;
		msz->ia[0] = 0; mdz->ia[0] = 0;
		Int i_ = 0;
		for (Int im = Mmin; im <= Mmax; im++) {
			for (Int i = 0; i < nz[alpha]; i++) {
				zi_ = zcGrid[alpha][i];
				zbasis[alpha]->getNonzero(zi_, limits);
				for (c = limits.begin(); c != limits.end(); c++) {
					j = *c;
					if (im != Mmin) {
						msz->aa[cntr1] = lamm[im-Mmin] * d1z[alpha][i][j];
						msz->ja[cntr1++] = (im - Mmin - 1)*nz[alpha] + j;
					}
					msz->aa[cntr1] = jm[im-Mmin]*sz[alpha][i][j];
					msz->ja[cntr1++] = (im-Mmin)*nz[alpha]+j;
					mdz->aa[cntr2] = mhat[im - Mmin] * zd1z[alpha][i][j] + \
												m2hat[im - Mmin] * sz[alpha][i][j];
					mdz->ja[cntr2++] = (im - Mmin)*nz[alpha] + j;
					if (im != Mmax) {
						msz->aa[cntr1] = -lamp[im - Mmin] * z2d1z[alpha][i][j] + \
							lamM[im - Mmin] * zsz[alpha][i][j];
						msz->ja[cntr1++] = (im - Mmin + 1)*nz[alpha] + j;
					}
				}
				msz->ia[++i_] = cntr1;
				mdz->ia[i_] = cntr2;
			}
		}

		//more 3d MATRICES FOR LHS OPERATOR (BLOCK ALPHA)
		shared_ptr<SparseMatr<Complex>> dxyz = \
			make_shared<SparseMatr<Complex>>( \
					nx[alpha]*ny[alpha]*nz[alpha], \
						(long long int)nnzx[alpha]*nnzy[alpha]*nnzz[alpha]);
		NumbersList limitsx, limitsy, limitsz;
		std::array<Int, 3> js;
		Complex vlxy;
		Int cntr = 0;
		dxyz->ia[0] = 0; i_ = 0;
		Int nxyz = nx[alpha] * ny[alpha] * nz[alpha];
		Int jx, jy, jz;
		Complex dsxyv;
		for (Int iz = 0; iz < nz[alpha]; iz++) {
			xyza.z = zi_ = zcGrid[alpha][iz];
			zbasis[alpha]->getNonzero(zi_, limitsz);
			Int cntr_xy = 0;
			for (Int iy = 0; iy < ny[alpha]; iy++) {
				xyza.y = yi = ycGrid[alpha][iy];
				ybasis[alpha]->getNonzero(yi, limitsy);
				for (Int ix = 0; ix < nx[alpha]; ix++) {
					xyza.x = xi = xcGrid[alpha][ix];
					xbasis[alpha]->getNonzero(xi, limitsx);
					vlxy = zzero;
					for (Int gamma = 0; gamma < nComp; gamma++) {
						sys.rotReduced(xyza, alpha, xyzb, gamma);
						xroti = xyzb.x; yroti = xyzb.y;
						if (gamma != alpha) {
							vlxy += (*vc[gamma])(xroti)*mc[gamma]->tail(xroti, yroti);
						}
					}
					for (Int gamma = nComp; gamma < 3; gamma++) {
						sys.rotReduced(xyza, alpha, xyzb, gamma);
						xroti = xyzb.x;
						vlxy += (*vc[gamma])(xroti) + (*vs[gamma])(xroti);
					}
					for (cy = limitsy.begin(); cy != limitsy.end(); cy++) {
						js[1] = jy = *cy;
						for (cx = limitsx.begin(); cx != limitsx.end(); cx++) {
							js[2] = jx = *cx;
							dsxyv = dsxdsy[alpha]->aa[cntr_xy] + \
											vlxy * sxsy[alpha]->aa[cntr_xy];
							for (cz = limitsz.begin(); cz != limitsz.end(); cz++) {
								js[0] = jz = *cz;
								j = jz*nx[alpha]*ny[alpha]+jy*nx[alpha]+jx;
								dxyz->aa[cntr] = sz[alpha][iz][jz] * dsxyv + \
									dsz0[alpha][iz][jz] * xr2yr2[alpha]->aa[cntr_xy];
								dxyz->ja[cntr] = j;
								cntr++;
							}
							cntr_xy++;
						}
					}
					dxyz->ia[++i_] = cntr;
				}
			}
		}

		std::vector<shared_ptr<AMatrix<Complex>>> lterms;
		lterms.reserve(3);
		std::array<shared_ptr<AMatrix<Complex>>, 2> ID;
		ID[0] = make_shared<IdentityM<Complex>>(nM);
		ID[1] = dxyz;
		std::array<shared_ptr<AMatrix<Complex>>, 2> MSSYS;
		MSSYS[0] = msz;
		MSSYS[1] = xr2sxsy[alpha];
		std::array<shared_ptr<AMatrix<Complex>>, 2> MDXSYS;
		MDXSYS[0] = mdz;
		MDXSYS[1] = xr2yr2[alpha];

		lterms.push_back(make_shared<TensorProd<Complex, 2>>(ID));
		lterms.push_back(make_shared<TensorProd<Complex, 2>>(MSSYS));
		lterms.push_back(make_shared<TensorProd<Complex, 2>>(MDXSYS));

		blocks.push_back(make_shared<MatrixSum<Complex>>(lterms));
		nnzTot += msz->nnonz() + mdz->nnonz() + \
			xr2sxsy[alpha]->nnonz() + xr2yr2[alpha]->nnonz() + \
				dxyz->nnonz();
		sizeGb += msz->sizeGb() + mdz->sizeGb() + \
			xr2sxsy[alpha]->sizeGb() + xr2yr2[alpha]->sizeGb() + \
			dxyz->sizeGb();
	}
	lhsrhs.push_back( \
		make_shared<BlockDiagMatr<Complex>>(blocks));

	//RHS TERM
	//CALCULATE TOTAL NUMBER OF NONZERO ELEMENTS
	{
		NumbersList limitsx, limitsy, limitsz;
		array<array<long long int, 3>, 3> nnz;
		for (Int alpha : solComps) {
			nnz[alpha] = { 0, 0, 0 };
			for (Int iz = 0; iz < nz[alpha]; iz++) {
				xyza.z = zcGrid[alpha][iz];
				for (Int iy = 0; iy < ny[alpha]; iy++) {
					xyza.y = ycGrid[alpha][iy];
					for (Int ix = 0; ix < nx[alpha]; ix++) {
						xyza.x = xcGrid[alpha][ix];
						//beta != alpha
						for (Int beta = 0; beta < nComp; beta++) {
							if (beta != alpha) {
								sys.rotReduced(xyza, alpha, xyzb, beta);
								xroti = xyzb.x; yroti = xyzb.y;
								//zroti = identical ? abs(xyzb.z) : xyzb.z;
								zroti = xyzb.z;
								zroti *= pmz[beta];
								xbasis[beta]->getNonzero(xroti, limsxrot);
								ybasis[beta]->getNonzero(yroti, limsyrot);
								zbasis[beta]->getNonzero(zroti, limszrot);
								nnz[alpha][beta] += \
									limsxrot.size()*limsyrot.size()*limszrot.size();
							}
						}
					}
				}
			}
		}

		//COMPUTE MATRIX
		Int nxyz;
		vector<vector<shared_ptr<AMatrix<Complex>>>> blocks;
		blocks.reserve(solComps.size());
		std::array<Int, 3> js;
		array<long long int, 3> cntr;
		double wba;
		Complex vsxydotxyz, f00mv;
		for (Int alpha : solComps) {
			blocks.push_back( \
				vector<shared_ptr<AMatrix<Complex>>>());
			Int i = 0;
			cntr = { 0, 0, 0 };
			nxyz = nx[alpha] * ny[alpha] * nz[alpha];
			vector<shared_ptr<SparseMatr<Complex>>> sab;
			vector<shared_ptr<TensorProd<Complex, 2>>> imsab;
			vector<shared_ptr<BlockMatr<Complex>>> fab;
			vector<vector<vector<shared_ptr<DiagMatrix<Complex>>>>> fabBl;
			vector<double> fabBlSizeGb;
			vector<Int> fabBlNNZ;
			sab.resize(nComp); fab.resize(nComp);
			fabBl.resize(nComp); imsab.resize(nComp);
			fabBlSizeGb.resize(nComp); fabBlNNZ.resize(nComp);
			for (Int beta = 0; beta < nComp; beta++) {
				if (beta != alpha) {
					sab[beta] = \
						make_shared<SparseMatr<Complex>>(nxyz, \
							nx[beta]*ny[beta]*nz[beta], nnz[alpha][beta]);
					sab[beta]->ia[0] = 0;
					fabBl[beta].resize(nM);
					for (Int im = Mmin; im <= Mmax; im++) {
						fabBl[beta][im - Mmin].reserve(nM);
						for (Int imbar = Mmin; imbar <= Mmax; imbar++)
							fabBl[beta][im - Mmin].push_back(\
								make_shared<DiagMatrix<Complex>>(nxyz));
					}
				}
			}
			Int cntr_xyzb = 0;
			for (Int iz = 0; iz < nz[alpha]; iz++) {
				xyza.z = zi_ = zcGrid[alpha][iz];
				zbasis[alpha]->getNonzero(zi_, limitsz);
				for (Int iy = 0; iy < ny[alpha]; iy++) {
					xyza.y = yi = ycGrid[alpha][iy];
					ybasis[alpha]->getNonzero(yi, limitsy);
					for (Int ix = 0; ix < nx[alpha]; ix++) {
						xyza.x = xi = xcGrid[alpha][ix];
						xbasis[alpha]->getNonzero(xi, limitsx);
						for (Int beta = 0; beta < nComp; beta++) {
							if (beta != alpha) { //beta != alpha
								sys.rotReduced(xyza, alpha, xyzb, beta);
								xroti = xyzb.x; yroti = xyzb.y;
								//zroti = identical ? abs(xyzb.z) : xyzb.z;
								zroti = xyzb.z;
								zroti *= pmz[beta];
								xbasis[beta]->getNonzero(xroti, limsxrot);
								ybasis[beta]->getNonzero(yroti, limsyrot);
								zbasis[beta]->getNonzero(zroti, limszrot);

								std::vector<Complex> sxyzdivxyrot;
								sxyzdivxyrot.reserve(limsxrot.size()*limsyrot.size()*limszrot.size());
								std::vector<Complex>::const_iterator cxyz;
								Complex szydivyrot;
								for (cz = limszrot.begin(); cz != limszrot.end(); cz++) {
									for (cy = limsyrot.begin(); cy != limsyrot.end(); cy++) {
										szydivyrot = szrot[alpha][cntr_xyzb][cz - limszrot.begin()]\
											*sydivyrot[alpha][cntr_xyzb][cy - limsyrot.begin()];
										for (cx = limsxrot.begin(); cx != limsxrot.end(); cx++) {
											sxyzdivxyrot.push_back(szydivyrot*sxdivxrot[alpha][cntr_xyzb][cx - limsxrot.begin()]);
										}
									}
								}

								cxyz = sxyzdivxyrot.begin();
								for (cz = limszrot.begin(); cz != limszrot.end(); cz++) {
									js[0] = *cz;
									for (cy = limsyrot.begin(); cy != limsyrot.end(); cy++) {
										js[1] = *cy;
										for (cx = limsxrot.begin(); cx != limsxrot.end(); cx++) {
											js[2] = *cx;
											j = js[0] * nx[beta] * ny[beta] + js[1] * nx[beta] + js[2];
											sab[beta]->aa[cntr[beta]] = *cxyz; //S(zrot)*( S(xrot)/xrot )*( S(yrot)/yrot )
											sab[beta]->ja[cntr[beta]] = j;
											cntr[beta]++;
											cxyz++;
										}
									}
								}

								wba = acos((sys.s[beta][alpha] * yi*zi_ + sys.c[beta][alpha] * xi) / xroti);
								if ((beta - alpha + 1) % 3 != 0)
									wba = 2 * PI - wba;
								for (Int im = Mmin; im <= Mmax; im++) {
									vsxydotxyz = \
										vsxydotxy[alpha][ix][iy] / pow(1.0 - zi_ * zi_, 0.5*im);
									for (Int imbar = Mmin; imbar <= Mmax; imbar++) {
										f00mv = f00pm[im - Mmin][imbar - Mmin] * wignerDSmall(J, imbar, im, wba) + \
											f00mm[im - Mmin][imbar - Mmin] * wignerDSmall(J, imbar, -im, wba);
										f00mv *= pm1M[beta][imbar - Mmin] * \
											vsxydotxyz*mz2powmhalf[alpha][cntr_xyzb][imbar - Mmin];
										//NB!  1/(xrot*yrot) in S_{alpha beta}
										(*fabBl[beta][im - Mmin][imbar - Mmin])[i] = f00mv;
									}
								}
								cntr_xyzb++;
								sab[beta]->ia[i + 1] = cntr[beta];
							}
						}
						i++;
					}
				}
			}
			
			array<shared_ptr<AMatrix<Complex>>, 2> imsab_;
			vector<vector<shared_ptr<AMatrix<Complex>>>> fabBl_;
			for (Int beta = 0; beta < nComp; beta++)
				if (beta != alpha) {//beta != alpha
					imsab_[0] = \
						make_shared<IdentityM<Complex>>(nM);
					imsab_[1] = sab[beta];
					fabBl_.resize(nM);
					for (Int im = Mmin; im <= Mmax; im++) {
						fabBl_[im-Mmin].reserve(nM);
						for (Int imbar = Mmin; imbar <= Mmax; imbar++)
							fabBl_[im-Mmin].push_back(fabBl[beta][im-Mmin][imbar-Mmin]);
					}
					imsab[beta] = \
						make_shared<TensorProd<Complex, 2>>(imsab_);
					fab[beta] = \
						make_shared<BlockMatr<Complex>>(fabBl_);
					fabBlSizeGb[beta] = fabBl[beta][0][0]->sizeGb();
					fabBlNNZ[beta] = fabBl[beta][0][0]->nrows();
					fabBl_.resize(0);
				}

			vector<vector<Int>> betas;
			betas.resize(solComps.size());
			Int beta_ = -1;
			for (Int beta = 0; beta < nComp; beta++)
				if (beta != alpha)
					if (binary_search(solComps.begin(), solComps.end(), beta))
						betas[++beta_].push_back(beta);
					else //beta not in solComps
						betas[beta_].push_back(beta);
				else //beta = alpha
					beta_++;
			blocks.back().reserve(solComps.size());
			shared_ptr<AMatrix<Complex>> tmp;
			for (beta_ = 0; beta_ < solComps.size(); beta_++) {
				if (betas[beta_].size() == 0)
					tmp = make_shared<TrivialM<Complex>>(nM*nxyz);
				else { //1 or 2 components
					if (betas[beta_].size() == 1) {//1 component
						vector<shared_ptr<AMatrix<Complex>>> mults;
						mults.reserve(2);
						mults.push_back(fab[betas[beta_][0]]);
						mults.push_back(imsab[betas[beta_][0]]);
						tmp = \
							make_shared<MatrixProd<Complex>>(mults);
						nnzTot += \
							fab[betas[beta_][0]]->nBlock * fabBlNNZ[betas[beta_][0]] + \
								sab[betas[beta_][0]]->nnonz();
						sizeGb += \
							fab[betas[beta_][0]]->nBlock * fabBlSizeGb[betas[beta_][0]] + \
								sab[betas[beta_][0]]->sizeGb();
					} else { // 2 components
						vector<shared_ptr<AMatrix<Complex>>> terms;
						terms.reserve(betas[beta_].size());
						for (auto gamma : betas[beta_]) {
							vector<shared_ptr<AMatrix<Complex>>> mults;
							mults.reserve(2);
							mults.push_back(fab[gamma]);
							mults.push_back(imsab[gamma]);
							terms.push_back(\
								make_shared<MatrixProd<Complex>>(mults));
							nnzTot += \
								fab[gamma]->nBlock * fabBlNNZ[gamma] + \
									sab[gamma]->nnonz();
							sizeGb += \
								fab[gamma]->nBlock * fabBlSizeGb[gamma] + \
									sab[gamma]->sizeGb();
						}
						tmp = make_shared<MatrixSum<Complex>>(terms);
					}
				}
				blocks.back().push_back(tmp);
			}
		}
		lhsrhs.push_back( \
			make_shared<BlockMatr<Complex>>(blocks));
	}
	matrNStor = new MatrixSum<Complex>(lhsrhs);
#endif
	cout << "Matrix size: " << \
		fixed << sizeGb << " Gbytes, " << \
			nnzTot << " nonzero elements" << endl << endl;
}

void FMOp::prepareIdentical() {
	if (sys.id == two_sym || sys.id == two_asym) {
		pmz[1] = -1.0;
		double p = (sys.id == two_sym) ? 1.0 : -1.0;
		for (Int im = Mmin; im <= Mmax; im++)
			pm1M[1][im - Mmin] = p * tau * pow(-1.0, J + im);
	}
}

void FMOp::getDiscrParams() {
	nx.reserve(nComp); ny.reserve(nComp); nz.reserve(nComp);
	xbasis.reserve(nComp); ybasis.reserve(nComp); zbasis.reserve(nComp);
	xcGrid.reserve(nComp); ycGrid.reserve(nComp); zcGrid.reserve(nComp);
	rx = 0; ry = 0; rz = 0;
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
		rz = max(rz, zbasis[alpha]->nFuncOnInterval());
		ry = max(ry, ybasis[alpha]->nFuncOnInterval());
		rx = max(rx, xbasis[alpha]->nFuncOnInterval());
	}
}

double FMOp::lambdaJM(const Int j, const Int m) const {
	return sqrt(j*(j+1.0)-m*(m+1.0));
}

FMOp::~FMOp(void) {
	if (matr != nullptr)
		delete matr;
	if (matrNStor != nullptr)
		delete matrNStor;
}
