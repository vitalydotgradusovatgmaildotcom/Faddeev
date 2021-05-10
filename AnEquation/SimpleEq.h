#pragma once

#include "ALinearEq.h"

template <class T>
class SimpleEq :
	public ALinearEq<T> {
public:
	SimpleEq(AnOperator<T> *op, \
		AFunction<T> *rhs, const shared_ptr<AFunction<T>> &sol);
	bool needsSolution() override;
	~SimpleEq();
	using ALinearEq<T>::solved;
};

template <class T>
SimpleEq<T>::SimpleEq(AnOperator<T> *op, \
	AFunction<T> *rhs, const shared_ptr<AFunction<T>> &sol) {
	this->op = op;
	this->rhs = rhs;
	this->sol = sol;
}

template <class T>
bool SimpleEq<T>::needsSolution() {
	return !solved;
}

template <class T>
SimpleEq<T>::~SimpleEq() {
	this->op = nullptr;
	this->rhs = nullptr;
	//this->sol = nullptr;
}