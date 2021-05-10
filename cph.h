#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <cctype>
#include <functional>
#include <memory>
#include <numeric>

using namespace std;

//#define NDEBUG //cancels assert
#include <cassert>
/*#define error(msg) {std::cout << "ERROR: " << msg << " in file " << \
	__FILE__ << " at line " << __LINE__ << std::endl; throw(1);}
*/

typedef complex<double> Complex;

#define MKL_Complex16 Complex
#include "mkl.h"
#include "mkl_lapacke.h"
#include "mkl_spblas.h"

typedef MKL_INT Int;

#define MALLOC(size) mkl_malloc(size, 64)
#define FREE(p) {mkl_free(p); p = nullptr;}

static const double NaN = numeric_limits<double>::quiet_NaN();
static const double Inf = numeric_limits<double>::infinity();
static const double eps = numeric_limits<double>::epsilon();
static const double maxval = numeric_limits<double>::max();
static const Int maxInt = numeric_limits<Int>::max();
static const Complex zzero = Complex(0.0, 0.0);
static const Complex zone = Complex(1.0, 0.0);
static const Complex zi = Complex(0.0, 1.0);
static const double PI = 3.14159265358979323846;
static const double bytes2Gbytes = 1.0 / (1 << 30);

/*inline double MAX(const double &x, const double &y) { \
	return x > y ? x : y;}
inline double MIN(const double &x, const double &y) { \
	return x < y ? x : y;}
inline Int MAX(const Int &x, const Int &y) { \
	return x > y ? x : y;}
inline Int MIN(const Int &x, const Int &y) { \
	return x < y ? x : y;}*/

template <class T>
inline void SWAP(T &a, T &b) {
	T tmp = b; b = a; a = tmp;
}

inline bool IS_EPS(const double x) {
	return abs(x) < 1e-15;
}

inline double CONJ(const double x) {
	return x;
}
inline Complex CONJ(const Complex z) {
	return conj(z);
}

template <class T>
inline Int SGN(const T x) {
	return (x > 0) - (x < 0);
}

inline Int DELTA(const Int n, const Int m) {
	return n == m ? 1 : 0;
}

inline double ABS(const double x) {
	return abs(x);
}

inline double ABS(const Complex z) {
	return LAPACKE_dlapy2(z.real(), z.imag());
}

template <class T> ostream& double52(ostream &stream);
template <class T> ostream& double135(ostream &stream);
template <class T> ostream& double2315(ostream &stream);

//NB! fixed notation, more than 5 positions possible
template <>
inline ostream& double52<double>(ostream &stream) {
	stream << fixed << setw(5) << setprecision(2);
	return stream;
}

//NB! fixed notation, more than 13 positions possible
template <>
inline ostream& double52<Complex>(ostream &stream) {
	stream << fixed << setw(13) << setprecision(2);
	return stream;
}

template<>
inline ostream& double135<double>(ostream &stream) {
	stream << scientific << setw(13) << setprecision(5);
	return stream;
}

template<>
inline ostream& double135<Complex>(ostream &stream) {
	stream << scientific << setw(29) << setprecision(5);
	return stream;
}

template<>
inline ostream& double2315<double>(ostream &stream) {
	stream << scientific << setw(23) << setprecision(15);
	return stream;
}

template<>
inline ostream& double2315<Complex>(ostream &stream) {
	stream << scientific << setw(49) << setprecision(15);
	return stream;
}

/*inline ostream& operator<<(ostream &stream, const Complex &z) {
	stream << "(" << double52<T> << z.real() << ", ";
	stream << double52<T> << z.imag() << ")";
	return stream;
}*/

inline void Header(const std::string &str) {
	cout << endl;
	cout << \
		"========================================" << endl;
	cout << "||   " << setw(33) << std::left << str << "||" << endl;
	cout << \
		"========================================" << endl;
	cout << endl;
}

inline void header(const std::string &str) {
	cout << endl;
	cout << str << endl;
	cout << \
		"----------------------------------------" << endl;
}

//===============================================
//STRING UTILS

//String split
//http://www.martinbroadhurst.com/how-to-split-a-string-in-c.html
template <class Container>
void split(const std::string& str, Container& cont,
	const std::string& delims = " ") {
	std::size_t current, previous = 0;
	current = str.find_first_of(delims);
	while (current != std::string::npos) {
		cont.push_back(str.substr(previous, current - previous));
		previous = current + 1;
		current = str.find_first_of(delims, previous);
	}
	cont.push_back(str.substr(previous, current - previous));
}

//trim string
// trim from start
static inline std::string & ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
		[](Int c) {return !std::isspace(c); }));
    return s;
}

// trim from end
static inline std::string & rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
		[](Int c) {return !std::isspace(c); }).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string & trim(std::string &s) {
    return ltrim(rtrim(s));
}

//===============================================
// note: this implementation does not disable this overload for array types
/*template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
template<class T, class Arg>
unique_ptr<T> make_unique(Arg&& arg) {
	return unique_ptr<T>(new T(std::forward<Arg>(arg)));
}*/