#pragma once

#include "cph.h"
#include "ALinearEq.h"
#include "Pair.h"
#include "FaddeevConfigurator.h"
#include "Collocation.h"
#include "ScaledPot.h"
#include "CWFCalculator.h"

#define WRITE_SOL_NPOI_RADEQ 2001

//Radial Schroedinger equation
//scattering calculation
class RadEq :
	public ALinearEq<Complex> {
public:
	RadEq(const Pair &pair, const Int l, const double ene);
	void writeSol() const;
	//void test() const;
	bool needsSolution() override;
	~RadEq();
	using ALinearEq<Complex>::solved;
protected:
	const FaddeevConfigurator &config;
	const Pair &pair;
	const Int l; //angular momentum
	const double ene; //energy
	double eConv;
	double k, eta;
	shared_ptr<const Collocation<Complex, 1, Complex>> discr;
	Int nx;
	const ABasis<Complex>* basis;
	//unique_ptr<ABasis<Complex>> vvv;
	Grid cGrid;
	shared_ptr<const APotential<Complex>> vs; //short-range potential
	shared_ptr<const APotential<Complex>> vc;
	CWFCalculator *coul;
	void makeDiscretization();
	void makePotentials();
	//scale from physical to computational coordinates
	void scaleCoo(double &x) const;
	void solve();
};

//The exact solution (the scattered wave) with sharply cut Coulomb tail:
	//CWFCalculator coul0(l, zzero);
	//Complex asmall = exp(zi*sigmaL(l, eta)) * \
		(coul->fl(k*2.5)*coul->dulp(k*2.5) - coul->dfl(k*2.5)*coul->ulp(k*2.5)) / \
		(coul0.fl(k*2.5)*coul->dulp(k*2.5) - coul0.dfl(k*2.5)*coul->ulp(k*2.5));
	//Complex abig = exp(zi*sigmaL(l, eta)) * \
		(coul->fl(k*2.5)*coul0.dfl(k*2.5) - coul->dfl(k*2.5)*coul0.fl(k*2.5)) / \
		(coul0.fl(k*2.5)*coul->dulp(k*2.5) - coul0.dfl(k*2.5)*coul->ulp(k*2.5));
	//if (xs[0] >= 2.5)
	//	val = abig*coul->ulp(k*xs[0]);
	//else
	//	val = asmall*coul0.fl(k*xs[0]) - coul->flExp(k*xs[0]);