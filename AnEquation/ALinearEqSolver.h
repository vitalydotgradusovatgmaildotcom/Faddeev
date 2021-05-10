#pragma once

#include "cph.h"
#include "ALinearEq.h"

#define DEFAULT_RELTOL_SOLV 1.0e-12

template <class T>
class ALinearEqSolver {
public:
	ALinearEqSolver(void);
	ALinearEqSolver(const ALinearEqSolver &rhs) = delete;
	ALinearEqSolver(ALinearEqSolver &&rhs) = delete;
	ALinearEqSolver & operator=(const ALinearEqSolver &rhs) = delete;
	ALinearEqSolver & operator=(ALinearEqSolver &&rhs) = delete;
	virtual void calculateSolution(ALinearEq<T> &eq) = 0;
	void setRelTol(const double reltol);
	virtual ~ALinearEqSolver(void) = default;
protected:
	double reltol;
};

template <class T>
ALinearEqSolver<T>::ALinearEqSolver(void) : reltol(DEFAULT_RELTOL_SOLV) { }

template <class T>
void ALinearEqSolver<T>::setRelTol(double reltol) {
	this->reltol = reltol;
}