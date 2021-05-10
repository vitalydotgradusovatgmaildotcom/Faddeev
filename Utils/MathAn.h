#pragma once

#include "cph.h"
#include "Algorithms.h"

template <class T>
void firstDeriv(vector<T> &f, const double h) {
	assert(f.size() >= 3);
	Int n = f.size();
	T last = 0.5*(f[n - 3] - 4.0 * f[n - 2] + 3.0 * f[n-1]) / h;
	T prev = f[0];
	f[0] = 0.5*(-3.0 * f[0] + 4.0 * f[1] - f[2]) / h;
	T curr;
	for (Int i = 1; i < n - 1; i++) {
		curr = 0.5*(f[i + 1] - prev) / h;
		prev = f[i];
		f[i] = curr;
	}
	f[n - 1] = last;
}

template <class T>
void secondDeriv(vector<T> &f, const double h) {
	assert(f.size() >= 4);
	Int n = f.size();
	double h2 = h * h;
	T last = (-f[n - 4] + 4.0*f[n - 3] - 5.0*f[n - 2] + 2.0*f[n - 1]) / h2;
	T prev = f[0];
	f[0] = (2.0*f[0] - 5.0*f[1] + 4.0*f[2] - f[3]) / h2;
	T curr;
	for (Int i = 1; i < n - 1; i++) {
		curr = (f[i + 1] - 2.0*f[i] + prev) / h2;
		prev = f[i];
		f[i] = curr;
	}
	f[n - 1] = last;
}

template <class T>
void secondDeriv(vector<double> &x, vector<T> &f) {
	assert(x.size() == f.size());
	assert(f.size() >= 4);
	Int n = f.size();
	double h, f_, g, alpha, beta, gamma, coef;

	T last[3];
	for (Int i = n - 3; i < n; i++) {
		h = x[i - 1] - x[i];
		g = x[i - 2] - x[i];
		f_ = x[i - 3] - x[i];
		alpha = -f_ * g*(g + f_) / ((h + g + f_) * (h - g) * (h - f_));
		beta = (f_ - alpha * (h - f_)) / (g - f_);
		gamma = -1.0 - alpha - beta;
		coef = 2.0 / (alpha*h*h + beta * g*g + gamma * f_*f_);
		last[i-(n-3)] = (f[i] + alpha * f[i - 1] + beta * f[i - 2] + gamma * f[i - 3]) / coef;
	}

	for (Int i = 0; i < n-3; i++) {
		h = x[i + 1] - x[i];
		g = x[i + 2] - x[i];
		f_ = x[i + 3] - x[i];
		alpha = -f_ * g*(g + f_) / ((h + g + f_) * (h - g) * (h - f_));
		beta = (f_ - alpha * (h - f_)) / (g - f_);
		gamma = -1.0 - alpha - beta;
		coef = 2.0 / (alpha*h*h + beta * g*g + gamma * f_*f_);
		f[i] = (f[i] + alpha * f[i + 1] + beta * f[i + 2] + gamma * f[i + 3]) / coef;
	}

	f[n - 3] = last[0]; f[n - 2] = last[1]; f[n - 1] = last[2];
}

template <class T>
T simpson(const vector<T> &f, const double h) {
	Int n = f.size();
	//NB! n must be odd!
	assert(n % 2 == 1);
	T sum = T();
	for (Int i = 1; i < n - 3; i += 2) {
		sum += f[i]*4.0*h / 3.0;
		sum += f[i+1]*2.0*h / 3.0;
	}
	sum += f[n -2]*4.0*h / 3.0;
	sum += (f[0] + f[n -1])*h / 3.0;
	return sum;
}

//NB! nodes must be sorted
inline double piecewiseLinear(const vector<double> &nodes, \
	const vector<double> &f, const double x) {
	Int ind1, ind2;
	findIntervSorted(nodes.begin(), nodes.end(), x, ind1, ind2);
	double h = nodes[ind2] - nodes[ind1];
	double dx = x - nodes[ind1];
	return f[ind1] + (f[ind2] - f[ind1])*dx / h;
}