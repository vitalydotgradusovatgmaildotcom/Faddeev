#pragma once

#include "FaddeevConfigurator.h"
#include "Pair.h"
#include "RadEigenProb.h"

#define CHI_NPOI_GRID_CALC 10001
//NB! must be odd for Simpson

class BoundGrCalculator {
public:
	BoundGrCalculator(const Pair &pair);
	BoundGrCalculator(const BoundGrCalculator &rhs) = delete;
	BoundGrCalculator(BoundGrCalculator &&rhs) = delete;
	BoundGrCalculator & operator=(const BoundGrCalculator &rhs) = delete;
	BoundGrCalculator & operator=(BoundGrCalculator &&rhs) = delete;
	~BoundGrCalculator() = default;
protected:
	Grid gr;
	shared_ptr<const APotential<Complex>> vs; //short-range potential
	shared_ptr<const APotential<Complex>> vc;
	const FaddeevConfigurator &config;
	const Pair &pair;
	double eConv;
	void makeDiscretization();
	void makePotentials();
	void makeChi();
	//scale from physical to computational coordinates
	void scaleCoo(double &x) const;
};

