#include "FMIdOp.h"
#include "TensorProd.h"
#include "BlockDiagMatr.h"

FMIdOp::FMIdOp( \
	const shared_ptr<const CompsDiscr> &cdiscr, \
		const System3Body &sys) \
		: VectOperator<Complex, 3>(cdiscr), cdiscr(cdiscr), matr(nullptr), \
				Mmin((1 - sys.tau) / 2), Mmax(sys.J), nComp(cdiscr->getNComp()), \
					solComps(cdiscr->getSolComps()), nM(Mmax-Mmin+1) {
	getDiscrParams();
	makeMatrix();
}

void FMIdOp::times(const AFunction<Complex> &u, AFunction<Complex> &res) {
	res.coef = *matr * u.coef;
}

void FMIdOp::solve(AFunction<Complex> &rhssol) {
	matr->solve(rhssol.coef);
}

void FMIdOp::makeMatrix() {
	std::vector<shared_ptr<AMatrix<Complex>>> blocks;
	blocks.reserve(solComps.size()*nM);
	
	for (Int alpha : solComps) {
		shared_ptr<AMatrix<Complex>> block = \
			makeMatrix(alpha);
		//NB! All blocks of component alpha are the same matrix in memory
		//TODO makeMatrix(const Int alpha, const Int im) ???
		for (Int im = Mmin; im <= Mmax; im++)
			blocks.push_back(block);
	}

	matr = new BlockDiagMatr<Complex>(blocks);
}

std::shared_ptr<AMatrix<Complex>> FMIdOp::makeMatrix(const Int alpha) {
	std::shared_ptr < GenMatrix<Complex>> sx, sy, sz;
	double zi, xi, yi; Int j;
	NumbersList limits;
	NumbersList::const_iterator c;

	sz = make_shared<GenMatrix<Complex>>(nz[alpha]);
	sz->fill(zzero);
	for (Int i = 0; i < nz[alpha]; i++) {
		zi = zcGrid[alpha][i];
		zbasis[alpha]->getNonzero(zi, limits);
		for (c = limits.begin(); c != limits.end(); c++) {
			j = *c;
			(*sz)[i][j] = zbasis[alpha]->f(zi, j);
		}
	}

	sx = make_shared<GenMatrix<Complex>>(nx[alpha]);
	sx->fill(zzero);
	for (Int i = 0; i < nx[alpha]; i++) {
		xi = xcGrid[alpha][i];
		xbasis[alpha]->getNonzero(xi, limits);
		for (c = limits.begin(); c != limits.end(); c++) {
			j = *c;
			(*sx)[i][j] = xbasis[alpha]->f(xi, j);
		}
	}

	sy = make_shared<GenMatrix<Complex>>(ny[alpha]);
	sy->fill(zzero);
	for (Int i = 0; i < ny[alpha]; i++) {
		yi = ycGrid[alpha][i];
		ybasis[alpha]->getNonzero(yi, limits);
		for (c = limits.begin(); c != limits.end(); c++) {
			j = *c;
			(*sy)[i][j] = ybasis[alpha]->f(yi, j);
		}
	}
	
	std::array<std::shared_ptr<AMatrix<Complex>>, 3> sxyz = {sz, sy, sx};

	return make_shared<TensorProd<Complex, 3>>(sxyz);
}

void FMIdOp::getDiscrParams() {
	nx.reserve(nComp); ny.reserve(nComp); nz.reserve(nComp);
	xbasis.reserve(nComp); ybasis.reserve(nComp); zbasis.reserve(nComp);
	xcGrid.reserve(nComp); ycGrid.reserve(nComp); zcGrid.reserve(nComp);
	for (Int alpha = 0; alpha < nComp; alpha++) {
		nz.push_back(cdiscr->get(alpha, Mmin).getNi(0));
		ny.push_back(cdiscr->get(alpha, Mmin).getNi(1));
		nx.push_back(cdiscr->get(alpha, Mmin).getNi(2));
		const Collocation<Complex, 3, Complex> &cd = \
			dynamic_cast<const Collocation<Complex, 3, Complex> &> \
				(cdiscr->get(alpha, Mmin));
		zbasis.push_back( cd.bases[0].get() );
		ybasis.push_back( cd.bases[1].get() );
		xbasis.push_back( cd.bases[2].get() );
		zcGrid.push_back( cd.cGrids[0] );
		ycGrid.push_back( cd.cGrids[1] );
		xcGrid.push_back( cd.cGrids[2] );
	}
}

FMIdOp::~FMIdOp() {
	if (matr != nullptr)
		delete matr;
}
