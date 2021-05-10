#pragma once

#include "Grid.h"
#include "MathAn.h"

inline void KvitsHu_log(Grid &r) {
	//mapping from Kvitsinsky, Hu, PRA, 47(2), 994 (1993)
	//close to equidistant in [0, 0.5*rmax], sparse near 0 and rmax
	double rmax = r.getRightmostPoint();
	double h = r[1] - r[0];
	for (Int i = 0; i < r.getNPoints(); i++)
		r[i] = rmax * log(1.0 - r[i] / (rmax + h)) / log(1.0 - rmax / (rmax + h));
}

inline void KvitsHu_pow(Grid &r, const double h) {
	//mapping from Kvitsinsky, Hu, PRA, 47(2), 994 (1993)
	double rmax = r.getRightmostPoint();
	for (Int i = 0; i < r.getNPoints(); i++)
		r[i] = h * (pow(1.0 + rmax / h, r[i] / rmax) - 1.0);
}

inline void KvitsHu_pow_my(Grid &r, const double alpha) {
	//the bigger is alpha > 1,
	//the denser is the grid near 0 and sparser near xmax
	//(mapping looks more like _|)
	assert(alpha > 1.0);
	double rmax = r.getRightmostPoint();
	for (Int i = 0; i < r.getNPoints(); i++)
		r[i] = rmax * (( pow(alpha, r[i] / rmax)-1.0 )/(alpha-1.0));
}

inline void Roudnev_pow(Grid &r, const double alpha) {
	//from VA java code
	//with alpha -> +infty, rmax fixed, u=r/rmax > 0
	//behaves as (rmax+1)^u - 1
	//with rmax -> +infty, alpha fixed, u > 0 as rmax*u^alpha
	//rmax * pow(r[i], alpha) - alternative?
	double rmax = r.getRightmostPoint();
	double k = pow(rmax + 1.0, 1.0 / alpha) - 1.0;
	for (Int i = 0; i < r.getNPoints(); i++) {
		r[i] = pow((k * (r[i] / rmax) + 1.0), alpha) - 1.0;
	}
}

inline void Kornev_sin(Grid &z) {
	//by Alexey Kornev
	for (Int i = 0; i < z.getNPoints(); i++)
		z[i] = sin(0.5 * z[i] * PI);
}

inline void Roudnev_tanh(Grid &z, const double alpha) {
	//from VA java code
	for (Int i = 0; i < z.getNPoints(); i++)
		z[i] = tanh(alpha*z[i]) / tanh(alpha);
}

inline void Roudnev_cubic(Grid &z, const double alpha) {
	//from VA java code
	for (Int i = 0; i < z.getNPoints(); i++)
		z[i] = (z[i] + alpha * z[i] * z[i] * z[i]) / (1.0 + alpha);
}

inline void mapChi(Grid &r, const vector<double> &xs, \
									const vector<double> &chi) {
	//use piecewise linear interpolation - monotonicity!
	assert(xs.size() == chi.size());
	Int n = r.getNPoints();
	double rmax = r.getRightmostPoint();
	for (Int i = 0; i < n; i++) {
		r[i] = rmax * piecewiseLinear(xs, chi, (double)(i) / (n - 1));
	}
}