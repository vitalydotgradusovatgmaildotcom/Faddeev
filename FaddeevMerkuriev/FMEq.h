#pragma once

#include "ALinearEq.h"
#include "FaddeevConfigurator.h"
#include "System3Body.h"
#include "CompsDiscr.h"
#include "ScaledPot.h"
#include "CoulombTail.h"
#include "MerkurievCutoff.h"
#include "LHSApprox.h"
#include "PrecFMOp.h"
#include "FMOp.h"
#include "BinChannel.h"

#define WRITE_PROJ_NPOI_FMEQ 2001
#define WRITE_SOL_NPOI_FMEQ 150

//#define POLY_BASIS_Z_FMEQ

//#define CORRECTION_FMEQ //hack

//system of Faddeev equations for scattering calculations
//discretized via spline collocation method
class FMEq :
	public ALinearEq<Complex> {
public:
	FMEq(const System3Body &sys, const double energy);
	//void renewEne(const double energy);
	bool needsSolution() override;
	void postprocess();
	void calcCrossSect();
	~FMEq(void);
protected:
	const FaddeevConfigurator &config;
	const System3Body &sys;
	double ene; //energy
	double eConv;
	Int nComp; //number of components in FM equations
	shared_ptr<const CompsDiscr> cdiscr;
	std::vector<Int> solComps; //solution components
	//number of equations and solutions solComps.size() <= nComp
	Int Mmin, Mmax, nM;
	array<shared_ptr<const APotential<Complex>>, 3> vs; //short-range potential
	array<shared_ptr<const APotential<Complex>>, 3> vc;
	array<const APotential<Complex> *, 3> vcy;
	array<const APotential<Complex> *, 3> vcy_sing;
	array<const MerkurievCutoff<Complex> *, 3> mc;
	std::vector<Int> nx, ny, nz;
	std::vector<const ABasis<Complex> *> xbasis, ybasis, zbasis;
	std::vector<Grid> xcGrid, ycGrid, zcGrid;
	BinChannel chan0;
	array<std::vector<BinChannel>, 3> chans;
	array<std::vector<BinChannel>, 3> chans_nond; //nondegenerate channels
	vector<vector<shared_ptr<Collocation<Complex, 1, Complex>>>> discry;
	vector<vector<shared_ptr<Collocation<Complex, 1, Complex>>>> discrz;
	std::vector<Int> lambdat;
	std::vector<Int>::const_iterator lamtit;
	std::array<std::vector<std::vector<Int>>, 3> lambda;
	std::array<std::vector<GenMatrix<Complex>>, 3> blamlam;
	void prepareIdentical();
	void makeChannels();
	void preparePartial();
	void makeDiscretization();
	Grid makeXGrid(const Int alpha, const double xmax);
	Grid makeYGrid(const Int alpha, const double ymax);
	Grid makeZGrid(const Int alpha);
	void makePotentials();
	void makeOperators();
	void makeRHS();
	void getDiscrParams();
	void projectBound();
	void writeProjections(const Int alpha, const Int mbar, const Int iCh, \
				const vector<Function<Complex, 1>> &projs);
	//void writeSolution(const Int alpha, const double z);
	void getAmpls(const Int alpha, const Int mbar, const Int iCh, \
		const vector<Function<Complex, 1>> &projs);
	//scale from physical to computational coordinates
	void scaleCoo(double &x, double &y, const Int alpha);
	double lambdaJM(const Int j, const Int m) const;
	//void polyCollocGridZ(const Int alpha, \
		Collocation<Complex, 3, Complex> &colloc);
};