#pragma once

#include "VectOperator.h"
#include "MerkurievCutoff.h"
#include "Components.h"
#include "APotential.h"
#include "Collocation.h"
#include "SparseMatr.h"
#include "MatrixSum.h"
#include "System3Body.h"
#include "SpecialFunctions.h"
#include "ABasis.h"

//#define STORE_MATRIX_FMOP

class FMOpPrec;

class FMOp :
	public VectOperator<Complex, 3> {
public:
	FMOp(const shared_ptr<const CompsDiscr> &cdiscr, \
		const array<shared_ptr<const APotential<Complex>>, 3> &vc, \
		const array<shared_ptr<const APotential<Complex>>, 3> &vs, \
		const array<const MerkurievCutoff<Complex> *, 3> &mc, \
		const System3Body &sys, double energy);
	Components operator*(const Components &u) const;
	void times(const AFunction<Complex> &u, AFunction<Complex> &res) override;
	void solve(AFunction<Complex> &rhssol) override;
	~FMOp(void);
protected:
	shared_ptr<const CompsDiscr> cdiscr;
	const System3Body &sys;
	Int nComp;
	Int J, Mproj, tau;
	Int Mmin, Mmax, nM;
	double E;
	array<shared_ptr<const APotential<Complex>>, 3> vs;
	array<shared_ptr<const APotential<Complex>>, 3> vc;
	array<const MerkurievCutoff<Complex> *, 3> mc;
	std::vector<Int> nx, ny, nz;
	Int rx, ry, rz; //highest possible number of nonzero basis functions
	//at any point of Interval of definition for all basis sets
	std::vector<const ABasis<Complex> *> xbasis, ybasis, zbasis;
	std::vector<Grid> xcGrid, ycGrid, zcGrid;
	SparseMatr<Complex> *matr;
	MatrixSum<Complex> *matrNStor;
	
	std::vector<Int> solComps;
	std::array<double, 3> pmz;
	std::array<std::vector<double>, 3> pm1M;
	long long int nnzTot = 0;
	double sizeGb = 0.0;
	void makeMatrs();
	void getDiscrParams();
	void prepareIdentical();
	double lambdaJM(const Int j, const Int m) const;
	friend class FMOpPrec;
};