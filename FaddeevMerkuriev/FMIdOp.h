#pragma once

#include "VectOperator.h"
#include "CompsDiscr.h"
#include "System3Body.h"
#include "Collocation.h"

class FMIdOp :
	public VectOperator<Complex, 3> {
public:
	FMIdOp( \
		const shared_ptr<const CompsDiscr> &cdiscr, \
			const System3Body &sys);
	void times(const AFunction<Complex> &u, AFunction<Complex> &res) override;
	void solve(AFunction<Complex> &rhssol) override;
	~FMIdOp();
protected:
	shared_ptr<const CompsDiscr> cdiscr;
	std::vector<Int> nx, ny, nz;
	std::vector<const ABasis<Complex> *> xbasis, ybasis, zbasis;
	std::vector<Grid> xcGrid, ycGrid, zcGrid;
	Int Mmin, Mmax, nM;
	Int nComp;
	AMatrix<Complex> *matr;
	vector<Int> solComps;
	void makeMatrix();
	std::shared_ptr<AMatrix<Complex>> makeMatrix(const Int alpha);
	void getDiscrParams();
};

