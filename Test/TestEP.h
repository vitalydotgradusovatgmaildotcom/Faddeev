#pragma once

#define TEST_EP_SPARSE_BIG

#include "AnEigenProblem.h"
#include "TestOperator.h"
#include "TestOpSparse.h"

#ifdef TEST_EP_SPARSE_BIG

class TestEP :
	public AnEigenProblem<TEST_OP_SP_TYPE> {
public:
	TestEP();
	~TestEP();
};

#else

class TestEP :
	public AnEigenProblem<TEST_OP_TYPE> {
public:
	TestEP();
	~TestEP();
};

#endif