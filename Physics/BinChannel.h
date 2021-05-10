#pragma once

#include "cph.h"
#include "Function.h"

class BinChannel {
public:
	Int alpha;
	Int n, l, m;
	double E;
	shared_ptr<Function<Complex, 1>> radWF;
	BinChannel() : alpha(0), n(0), l(0), m(0), E(0.0) { }
	BinChannel(Int alpha, Int n, Int l, Int m, double E, \
		shared_ptr<Function<Complex, 1>> radWF) \
		: alpha(alpha), n(n), l(l), m(m), E(E), radWF(radWF) { }
	BinChannel(const BinChannel &rhs) = default;
	BinChannel(BinChannel &&rhs) = default;
	BinChannel & operator=(const BinChannel &rhs) = default;
	BinChannel & operator=(BinChannel &&rhs) = default;
	void print() {
		cout << "pair " << alpha << " n = " << n \
			<< " l = " << l << " energy = " \
			<< double135<double> << E << endl;
	}
	bool operator==(const BinChannel &rhs) const {
		return (alpha == rhs.alpha) && (n == rhs.n) && (l == rhs.l) && (m == rhs.m);
	}
	~BinChannel(void) = default;
};

