#pragma once

template <class T>
class MerkurievCutoff {
public:
	MerkurievCutoff(const double nu, const double x0, const double y0);
	MerkurievCutoff(const MerkurievCutoff &rhs) = delete;
	MerkurievCutoff(const MerkurievCutoff &&rhs) = delete;
	MerkurievCutoff & operator=(const MerkurievCutoff &rhs) = delete;
	MerkurievCutoff & operator=(const MerkurievCutoff &&rhs) = delete;
	inline T operator()(const T x, const T y) const;
	inline T tail(const T x, const T y) const;
	~MerkurievCutoff(void) = default;
protected:
	double nu, x0, y0;
};

template <class T>
MerkurievCutoff<T>::MerkurievCutoff( \
	const double nu, const double x0, const double y0) \
		: nu(nu), x0(x0), y0(y0) { }

template<>
inline double MerkurievCutoff<double>::operator()(const double x, const double y) const {
	//use cutoff function from
	//Kvitsinsky, Carbonell, Gignoux, PRA, 46, 1310 (1992)
	return 2.0 / (1.0 + exp(pow(x / x0, nu) / (1.0 + y / y0)));
}

template<>
inline Complex MerkurievCutoff<Complex>::operator()(const Complex x, const Complex y) const {
	//use cutoff function from
	//Kvitsinsky, Carbonell, Gignoux, PRA, 46, 1310 (1992)
	Complex pwr = pow(x / x0, nu) / (1.0 + y / y0);
	if (pwr.real() > 700.0)
		return zzero;
	return 2.0 / (1.0 + exp(pwr));
}

template<class T>
inline T MerkurievCutoff<T>::tail(const T x, const T y) const {
	return 1.0 - (*this)(x, y);
}