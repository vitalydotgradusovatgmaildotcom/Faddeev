#pragma once

#include "cph.h"

double pLegendre(const Int l, const Int m, const double x);

double pLegendre(const Int l, const double x);
double dpLegendre(const Int l, const double x);
double ddpLegendre(const Int l, const double x);
void pLegendre(const Int l, const double x, std::vector<double> &vals);
void dpLegendre(const Int l, const double x, std::vector<double> &vals);
void ddpLegendre(const Int l, const double x, std::vector<double> &vals);

void gauleg(const double x1, const double x2, \
	vector<double> &x, vector<double> &w);

Complex wignerD(const Int j, const Int m, const Int mbar, \
	const double alpha, const double beta, const double gamma);

double wignerDSmall(const Int j, const Int m, const Int mbar, \
	const double beta);

//void pLegZeros(const Int m, const Int l, vector<double> &zeros);

inline double factorial(const Int n) {
	double res = 1.0;
	for (Int i = 2; i <= n; i++)
		res *= i;
	return res;
}

inline double double_factorial(const Int n) {
	double res = 1;
	for (Int i = n; i > 1; i-=2)
		res *= i;
	return res;
}

//inline double gleg(const Int l, const Int m, const double x) {
//	assert(x > -1.0 && x < 1.0);
//	return sqrt(2*PI)*pLegendre(l, m, x) / pow(1.0-x*x, 0.5*m);
//}