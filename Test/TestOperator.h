#pragma once

#include "Operator.h"

#define TEST_OP_TYPE double //Complex

class TestOperator :
	public Operator<TEST_OP_TYPE, 1>
{
public:
	GenMatrix<TEST_OP_TYPE> matr;
	TestOperator(const shared_ptr<const ADiscretization<TEST_OP_TYPE, 1>> &discr);
	Function<TEST_OP_TYPE, 1> operator*(const Function<TEST_OP_TYPE, 1> &u) const;
	void times(const AFunction<TEST_OP_TYPE> &u, AFunction<TEST_OP_TYPE> &res) override;
	void solve(AFunction<Complex> &rhssol) override;
	~TestOperator(void);
protected:
	//Function<TEST_OP_TYPE, 1> & resFun;
};