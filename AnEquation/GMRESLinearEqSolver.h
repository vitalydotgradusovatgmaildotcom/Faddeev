
#pragma once

#include "ALinearEqSolver.h"
#include "KrylovSubspace.h"
#include "UTMatrix.h"

#define KRYLOV_DIM_MAX_GMRES 700
#define RESTART_NUM_MAX_GMRES 40
//NB after restart, stagnation is possible!

template <class T>
class GMRESLinearEqSolver :
	public ALinearEqSolver<T> {
public:
	GMRESLinearEqSolver(void);
	void calculateSolution(ALinearEq<T> &eq) override;
	//using initial guess x0 = eq.sol
	~GMRESLinearEqSolver(void);
protected:
	using ALinearEqSolver<T>::reltol;
	void makeRot();
	Int dimmax;
	std::vector<T> si; //double in fact!
	std::vector<T> ci;
	UTMatrix<T> Hmrot;
	Vector<T> gmrot;
	KrylovSubspace<T> * ks;
};

