#pragma once

#include "AnEigenSolver.h"
#include "KrylovSubspace.h"

//Implicitly Restarted Arnoldi's Method
//Finds ev A x = lambda x with largest/lowest moduli/real parts
//No Shift-and-Invert atm

#define KRYLOV_DIM_MAX_IRAM 20 //30
#define ARNOLDI_ITER_MAX_IRAM 200

template <class T>
class IRAMEigenSolver :
	public AnEigenSolver<T> {
public:
	IRAMEigenSolver();
	IRAMEigenSolver(const IRAMEigenSolver &rhs) = delete;
	IRAMEigenSolver(IRAMEigenSolver &&rhs) = delete;
	IRAMEigenSolver & operator=(const IRAMEigenSolver &rhs) = delete;
	IRAMEigenSolver & operator=(IRAMEigenSolver &&rhs) = delete;
	void getEEV(AnEigenProblem<T> &ep, \
		const Int nev, const  probType &task, const bool computeVec) override;
	//Look comments to AnEigenSolver::getEEV
	//ep.evec[0] contains initial vector for Arnoldi algorithm
	//	it MUST have size op.getRank()
	//	(if it has zero norm then it is randomly generated)
	//NB! Computation of eigenvectors is not realized atm
	//NB! The number of converged eigenvalues can be both less and greater
	//	(by one in double case) than nev
	~IRAMEigenSolver();
protected:
	Int m = 0;
	Int k0 = 0; //initial number of "wanted" eigenvalues
	Int k = 0; //current number of "wanted" eigenvalues
	Int q = 0; //current number of shifts per iteration
	//(= number of "unwanted" eigenvalues)
	Int nconv = 0; //number of converged "wanted" ev
	KrylovSubspace<T> * ks;
	const double eps23;
	const double epsx2; //base*eps
	double small; //measure of small number for given task
	probType task;
	Vector<Complex> ritz;
	Vector<double> resids;
	void getHEigAndResid();
	void reorderEig();
	void getConvWantedEig();
	void applyShifts();
};

