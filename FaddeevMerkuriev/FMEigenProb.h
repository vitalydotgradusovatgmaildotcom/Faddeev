#pragma once

#include "AnEigenProblem.h"
#include "FaddeevConfigurator.h"
#include "System3Body.h"
#include "CompsDiscr.h"
#include "ScaledPot.h"
#include "CoulombTail.h"
#include "MerkurievCutoff.h"
#include "PrecFMOp.h"
#include "FMIdOp.h"
#include "LHSApprox.h"

//#define DIRECT_SOLV_FMEIGEN
//choose direct or iterative with preconditioning solution method

//#define ILU0_PREC_FMEIGEN

//#define POLY_BASIS_Z_FMEIGEN

//Eigenvalue problem for Faddeev operator
//discretized via spline collocation method
class FMEigenProb :
	public AnEigenProblem<Complex> {
public:
	FMEigenProb(const System3Body &sys, const double energy);
	Vector<Complex> getEval() const override;
	shared_ptr<AFunction<Complex>> getEvec(const Int i) override;
	~FMEigenProb();
protected:
	const FaddeevConfigurator &config;
	const System3Body &sys;
	double ene; //energy
	double eConv;
	Int nComp; //number of components in FM equations
	shared_ptr<const CompsDiscr> cdiscr;
	std::vector<Int> solComps; //solution components
	//number of equations and solutions solComps.size() <= nComp
	array<shared_ptr<const APotential<Complex>>, 3> vs; //short-range potential
	array<shared_ptr<const APotential<Complex>>, 3> vc;
	array<const APotential<Complex> *, 3> vcy;
	array<const MerkurievCutoff<Complex> *, 3> mc;
	PrecFMOp *precFM;
	FMIdOp *idop;
	void prepareIdentical();
	void makeDiscretization();
	Grid makeXGrid(const Int alpha, const double xmax);
	Grid makeYGrid(const Int alpha, const double ymax);
	Grid makeZGrid(const Int alpha);
	void makePotentials();
	void makeOperators();
	//scale from physical to computational coordinates
	void scaleCoo(double &x, double &y, const Int alpha) const;
	//void polyCollocGridZ(const Int alpha, \
		Collocation<Complex, 3, Complex> &colloc);
};

//S * (H-E)^(-1)
class IterativeOp :
	public VectOperator<Complex, 3> {
public:
	IterativeOp(const shared_ptr<const CompsDiscr> &cdiscr, \
		PrecFMOp *precFM, FMIdOp *idop);
	void times(const AFunction<Complex> &u, \
		AFunction<Complex> &res) override;
	void solve(AFunction<Complex> &rhssol) override;
	~IterativeOp(void) = default;
protected:
	PrecFMOp *precFM;
	FMIdOp *idop;
};