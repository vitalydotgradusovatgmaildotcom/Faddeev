#pragma once

#include "ALinearEq.h"
#include "TestOperator.h"

class TestEq :
	public ALinearEq<TEST_OP_TYPE> {
public:
	TestEq(void);
	~TestEq(void);
};

