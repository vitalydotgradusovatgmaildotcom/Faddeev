#pragma once

#include "ADiscretization.h"

//T defines whether the solution is real/complex valued
//TBas defines whether basis functions are real/complex valued
template <class T, Int d, class TBas>
class AProjDiscr :
	public ADiscretization<T, d> {
public:
	std::array<shared_ptr<const ABasis<TBas>>, d> bases;
	AProjDiscr(array<shared_ptr<const ABasis<TBas>>, d> &bases);
	T f(const std::array<double, d> &xs, const AFunction<T> &func) const override;
	T df(const std::array<double, d> &xs, const AFunction<T> &func) const override;
	T d2f(const std::array<double, d> &xs, const AFunction<T> &func) const override;
	using ADiscretization<T, d>::getRaw;
	~AProjDiscr(void) = default;
protected:
	using ADiscretization<T, d>::n;
	using ADiscretization<T, d>::ns;
};

template <class T, Int d, class TBas>
AProjDiscr<T, d, TBas>::AProjDiscr(array<shared_ptr<const ABasis<TBas>>, d> &bases) {
	for (Int k = 0; k < bases.size(); k++)
		this->bases[k] = bases[k];
	n = 1;
	for (Int k = 0; k < bases.size(); k++) {
		//this->bases[k] = std::move(bases[k]);
		ns[k] = this->bases[k]->getNCoef();
		n *= ns[k];
	}
}

template <>
double AProjDiscr<double, 1, double>::f(const std::array<double, 1> &xs, \
			const AFunction<double> &func) const {
	return bases[0]->f(xs[0], func);
}

template <>
Complex AProjDiscr<Complex, 1, double>::f(const std::array<double, 1> &xs, \
			const AFunction<Complex> &func) const {
	return bases[0]->f(xs[0], func);
}

template <>
Complex AProjDiscr<Complex, 1, Complex>::f(const std::array<double, 1> &xs, \
	const AFunction<Complex> &func) const {
	return bases[0]->f(xs[0], func);
}

template <class T, Int d, class TBas>
T AProjDiscr<T, d, TBas>::f(const std::array<double, d> &xs, \
			const AFunction<T> &func) const {
//NB!!! only zero boundary conditions if d > 1 !!!
	array<NumbersList, d> nls;
	array<NumbersList::const_iterator, d> nums;
	Int nnonz = 1;
	array<vector<TBas>, d> fs;
	for (Int k = 0; k < d; k++) {
		bases[k]->getNonzero(xs[k], nls[k]);
		if ( nls[k].empty() ) return T();
		nums[k] = nls[k].begin();
		nnonz *= nls[k].size();
		fs[k].reserve(nls[k].size());
		for (Int i = 0; i < nls[k].size(); i++)
			fs[k].push_back(bases[k]->f(xs[k], nls[k][i]));
	}
	nls[0].push_back(-1); nums[0] = nls[0].begin();

	Int dim = d-1;
	array<T, d> res; res.fill(T());
	array<Int, d> is;
	for (Int k = 0 ; k < d; k++)
		is[k] = *(nums[k]);
	for (Int i = 0; i < nnonz; i++) {
		//dim == d - 1 atm
		res[dim] += func[getRaw(is)]*fs[dim][nums[dim]-nls[dim].begin()];
		while (nums[dim] == nls[dim].end()-1) {
			nums[dim] = nls[dim].begin();
			is[dim--] = *(nums[dim]);
			res[dim] += res[dim+1]*fs[dim][nums[dim]-nls[dim].begin()];
			res[dim+1] = T();
		}
		//nums[dim] != nums[dim].end() now
		nums[dim]++; is[dim] = *(nums[dim]);
		dim = d - 1;
	}

	return res[0];
}

template <>
double AProjDiscr<double, 1, double>::df(const std::array<double, 1> &xs, \
			const AFunction<double> &func) const {
	return bases[0]->d(xs[0], func);
}

template <>
Complex AProjDiscr<Complex, 1, double>::df(const std::array<double, 1> &xs, \
			const AFunction<Complex> &func) const {
	return bases[0]->d(xs[0], func);
}

template <>
Complex AProjDiscr<Complex, 1, Complex>::df(const std::array<double, 1> &xs, \
	const AFunction<Complex> &func) const {
	return bases[0]->d(xs[0], func);
}

template <class T, Int d, class TBas>
T AProjDiscr<T, d, TBas>::df(const std::array<double, d> &xs, const AFunction<T> &func) const {
	//NB!!! only zero boundary conditions !!!
	array<NumbersList, d> nls;
	array<NumbersList::const_iterator, d> nums;
	Int nnonz = 1;
	array<vector<TBas>, d> fs;
	for (Int k = 0; k < d; k++) {
		bases[k]->getNonzero(xs[k], nls[k]);
		nums[k] = nls[k].begin();
		nnonz *= nls[k].size();
		fs[k].reserve(nls[k].size());
		for (Int i = 0; i < nls[k].size(); i++)
			fs[k].push_back(bases[k]->d(xs[k], nls[k][i]));
	}
	nls[0].push_back(-1);

	Int dim = d-1;
	array<T, d> res; res.fill(T());
	array<Int, d> is;
	for (Int k = 0 ; k < d; k++)
		is[k] = *(nums[k]);
	for (Int i = 0; i < nnonz; i++) {
		//dim == d - 1 atm
		res[dim] += func[getRaw(is)]*fs[dim][nums[dim]-nls[dim].begin()];
		while (nums[dim] == nls[dim].end()) {
			nums[dim] = nls[dim].begin();
			is[dim--] = *(nums[dim]);
			res[dim] += res[dim+1]*fs[dim][nums[dim]-nls[dim].begin()];
			res[dim+1] = T();
		}
		//nums[dim] != nums[dim].end() now
		nums[dim]++; is[dim] = *(nums[dim]);
	}

	return res[0];
}

template <>
double AProjDiscr<double, 1, double>::d2f(const std::array<double, 1> &xs, \
			const AFunction<double> &func) const {
	return bases[0]->dd(xs[0], func);
}

template <>
Complex AProjDiscr<Complex, 1, double>::d2f(const std::array<double, 1> &xs, \
			const AFunction<Complex> &func) const {
	return bases[0]->dd(xs[0], func);
}

template <>
Complex AProjDiscr<Complex, 1, Complex>::d2f(const std::array<double, 1> &xs, \
	const AFunction<Complex> &func) const {
	return bases[0]->dd(xs[0], func);
}

template <class T, Int d, class TBas>
T AProjDiscr<T, d, TBas>::d2f(const std::array<double, d> &xs, const AFunction<T> &func) const {
	//NB!!! only zero boundary conditions !!!
	array<NumbersList, d> nls;
	array<NumbersList::const_iterator, d> nums;
	Int nnonz = 1;
	array<vector<TBas>, d> fs;
	for (Int k = 0; k < d; k++) {
		bases[k]->getNonzero(xs[k], nls[k]);
		nums[k] = nls[k].begin();
		nnonz *= nls[k].size();
		fs[k].reserve(nls[k].size());
		for (Int i = 0; i < nls[k].size(); i++)
			fs[k].push_back(bases[k]->dd(xs[k], nls[k][i]));
	}
	nls[0].push_back(-1);

	Int dim = d-1;
	array<T, d> res; res.fill(T());
	array<Int, d> is;
	for (Int k = 0 ; k < d; k++)
		is[k] = *(nums[k]);
	for (Int i = 0; i < nnonz; i++) {
		//dim == d - 1 atm
		res[dim] += func[getRaw(is)]*fs[dim][nums[dim]-nls[dim].begin()];
		while (nums[dim] == nls[dim].end()) {
			nums[dim] = nls[dim].begin();
			is[dim--] = *(nums[dim]);
			res[dim] += res[dim+1]*fs[dim][nums[dim]-nls[dim].begin()];
			res[dim+1] = T();
		}
		//nums[dim] != nums[dim].end() now
		nums[dim]++; is[dim] = *(nums[dim]);
	}

	return res[0];
}
