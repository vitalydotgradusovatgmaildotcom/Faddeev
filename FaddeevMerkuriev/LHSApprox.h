#pragma once

#include "VectOperator.h"
#include "CompsDiscr.h"
#include "Collocation.h"
#include "APotential.h"
#include "MatrixOp.h"
#include "System3Body.h"

//#define FILTER_EIGS_LHSAPPROX
#define FILTER_EIGS_THRESH 0.05

class LHSApprox :
	public VectOperator<Complex, 3> {
public:
	LHSApprox( \
		const shared_ptr<const CompsDiscr> &cdiscr, \
		const array<shared_ptr<const APotential<Complex>>, 3> &vc,
		const array<shared_ptr<const APotential<Complex>>, 3> &vs,
		const array<const APotential<Complex> *, 3> &vcy, \
		const System3Body &sys, double energy);
	void times(const AFunction<Complex> &u, AFunction<Complex> &res) override;
	void solve(AFunction<Complex> &rhssol) override;
	~LHSApprox(void);
protected:
	AMatrix<Complex> *matr;
	shared_ptr<const CompsDiscr> cdiscr;
	const System3Body &sys;
	Int J, Mproj, tau;
	Int Mmin, Mmax, nM;
	Int nComp;
	double E;
	std::vector<Int> nx, ny, nz;
	std::vector<const ABasis<Complex> *> xbasis, ybasis, zbasis;
	std::vector<Grid> xcGrid, ycGrid, zcGrid;
	array<shared_ptr<const APotential<Complex>>, 3> vs;
	array<shared_ptr<const APotential<Complex>>, 3> vc;
	array<const APotential<Complex> *, 3> vcy;
	vector<Int> solComps;
	double sizeGb = 0.0;
	void makeMatrix();
	shared_ptr<AMatrix<Complex>> makeMatrix(const Int alpha);
	void makeZoperators(const Int alpha, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &wbar, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &w);
	void makeYoperator(const Int alpha, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &wbar, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &w, \
			std::vector<shared_ptr<DiagMatrix<Complex>>> &lambda);
	void makeXoperator(const Int alpha, \
			const std::vector<shared_ptr<DiagMatrix<Complex>>> &lambday, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &wbar, \
			std::vector<shared_ptr<GenMatrix<Complex>>> &w, \
			std::vector<shared_ptr<DiagMatrix<Complex>>> &lambda);
	double lambdaJM(const Int j, const Int m) const;
	void getDiscrParams();
};

