#pragma once

#include "AnEigenProblem.h"

#define DEFAULT_RELTOL_EIGEN 1.0e-14

enum probType { LarM, LowM, LarR, LowR };

template <class T>
class AnEigenSolver {
public:
	AnEigenSolver();
	AnEigenSolver(const AnEigenSolver &rhs) = delete;
	AnEigenSolver(AnEigenSolver &&rhs) = delete;
	AnEigenSolver & operator=(const AnEigenSolver &rhs) = delete;
	AnEigenSolver & operator=(AnEigenSolver &&rhs) = delete;
	void setRelTol(const double reltol);
	virtual void getEEV(AnEigenProblem<T> &ep, \
		const Int nev, const  probType &task, const bool computeVec) = 0;
	//nev - number of wanted eigenvalues
	//task - defines wanted eigenvalues
	//task = LarM/LowM - largest/lowest moduli, LarR/LowR - largest/lowest real parts
	//If computeVec = true, eigenvectors are to be found
	//On exit, ep.eval and ep.evec contain eigenvalues and corresponding eigenvectors
	virtual ~AnEigenSolver() = default;
protected:
	double reltol;
};

template <class T>
AnEigenSolver<T>::AnEigenSolver() : reltol(DEFAULT_RELTOL_EIGEN) {}

template <class T>
void AnEigenSolver<T>::setRelTol(const double reltol) {
	this->reltol = reltol;
}