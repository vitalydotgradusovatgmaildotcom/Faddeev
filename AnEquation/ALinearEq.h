#pragma once

#include "AnOperator.h"
#include "AFunction.h"


template <class T>
class ALinearEq {
public:
	AnOperator<T> *op;
	AFunction<T> *rhs;
	shared_ptr<AFunction<T>> sol;
	bool solved = false;
	ALinearEq();
	ALinearEq(const ALinearEq &rhs) = delete;
	ALinearEq(ALinearEq &&rhs) = delete;
	ALinearEq & operator=(const ALinearEq &rhs) = delete;
	ALinearEq & operator=(ALinearEq &&rhs) = delete;
	void getResidue(AFunction<T> &residue) const;
	bool isSolved() const;
	virtual bool needsSolution() = 0;
	virtual ~ALinearEq(void);
};

template <class T>
ALinearEq<T>::ALinearEq() : op(nullptr), rhs(nullptr), sol(nullptr) { }

template <class T>
void  ALinearEq<T>::getResidue(AFunction<T> &residue) const {
	op->times(*sol, residue);
	residue -= *rhs; residue *= -1.0;
}

template <class T>
bool ALinearEq<T>::isSolved() const {
	return solved;
}

template <class T>
ALinearEq<T>::~ALinearEq(void) {
	if (op != nullptr)
		delete op;
	if (rhs != nullptr)
		delete rhs;
	//if (sol != nullptr)
	//	delete sol;
}