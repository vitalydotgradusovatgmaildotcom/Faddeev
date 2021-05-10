#pragma once

#include "AnOperator.h"

template <class T>
class AnEigenProblem {
public:
	AnOperator<T> *op;
	Vector<Complex> eval;
	std::vector<shared_ptr<AFunction<T>>> evec;
	AnEigenProblem();
	AnEigenProblem(const AnEigenProblem &rhs) = delete;
	AnEigenProblem(AnEigenProblem &&rhs) = delete;
	AnEigenProblem & operator=(const AnEigenProblem &rhs) = delete;
	AnEigenProblem & operator=(AnEigenProblem &&rhs) = delete;
	virtual Vector<Complex> getEval() const = 0;
	virtual shared_ptr<AFunction<T>> getEvec(const Int i) = 0;
	virtual ~AnEigenProblem();
};

template <class T>
AnEigenProblem<T>::AnEigenProblem() : op(nullptr) {}

template <class T>
AnEigenProblem<T>::~AnEigenProblem() {
	if (op != nullptr)
		delete op;
	//for (Int k = 0; k < evec.size(); k++)
	//	if (evec[k] != nullptr)
	//		delete evec[k];
}