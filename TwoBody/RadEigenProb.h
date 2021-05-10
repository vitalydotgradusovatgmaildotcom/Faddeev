#pragma once

#include "AnEigenProblem.h"
#include "Pair.h"
#include "FaddeevConfigurator.h"
#include "Collocation.h"
#include "ScaledPot.h"

#define WRITE_SOL_NPOI_RADEP 2001
#define IM_EV_THRESHOLD 1.0e-12

//Radial Schroedinger equation
//returns negative real eigenvalues first
//with eigenvectors normalized
class RadEigenProb :
	public AnEigenProblem<Complex> {
public:
	RadEigenProb(const Pair &pair, const Int l);
	Vector<Complex> getEval() const override;
	shared_ptr<AFunction<Complex>> getEvec(const Int i) override;
	Int numberOfBound() const;
	void writeSolution() const;
	double rms(const Int i) const; //in physical coordinates
	~RadEigenProb() = default;
protected:
	const FaddeevConfigurator &config;
	const Pair &pair;
	const Int l; //angular momentum
	double eConv;
	shared_ptr<const Collocation<Complex, 1, double>> discr;
	Int nx;
	const ABasis<double>* basis;
	Grid cGrid;
	shared_ptr<const APotential<Complex>> vs; //short-range potential
	shared_ptr<const APotential<Complex>> vc;
	Int nBound = 0; //number of bound states
	void makeDiscretization();
	void makePotentials();
	//scale from physical to computational coordinates
	void scaleCoo(double &x) const;
	void solve();
	void normalizeBound();
};

