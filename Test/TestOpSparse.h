#pragma once

#include "Operator.h"
#include "SparseMatr.h"
//#include "SparseMatr_old.h"

#define TEST_OP_SP_TYPE  double //Complex

class TestOpSparse :
	public Operator<TEST_OP_SP_TYPE, 1> {
public:
	SparseMatr<TEST_OP_SP_TYPE> matr;
	TestOpSparse(const shared_ptr<const ADiscretization<TEST_OP_SP_TYPE, 1>> &discr);
	Function<TEST_OP_SP_TYPE, 1> operator*(const Function<TEST_OP_SP_TYPE, 1> &u) const;
	void times(const AFunction<TEST_OP_SP_TYPE> &u, AFunction<TEST_OP_SP_TYPE> &res) override;
	void solve(AFunction<Complex> &rhssol);
	~TestOpSparse(void);
};

