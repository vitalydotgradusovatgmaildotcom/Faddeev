#include "PermutationM.h"

template class PermutationM<double>;
template class PermutationM<Complex>;

template <class T>
PermutationM<T>::PermutationM(void) { }

template <class T>
PermutationM<T>::PermutationM(std::vector<Int> &&perm) \
			: p(std::move(perm)) {
//The permutation is specified by its full form
//(1 2 ... n)
//(   perm  )
	this->n = p.size();
	this->m = n;
}

template <class T>
PermutationM<T> & PermutationM<T>::operator=(const PermutationM<T> &rhs) {
	if (this != &rhs) {
		n = rhs.n; m = rhs.m;
		p = rhs.p;
	}
	return *this;
}

template <class T>
PermutationM<T> & PermutationM<T>::operator=(const AMatrix<T> &rhs) {
	const PermutationM<T> & rhs_ = dynamic_cast<const PermutationM<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
PermutationM<T> & PermutationM<T>::operator=(PermutationM<T> &&rhs) {
	n = rhs.n; rhs.n = 0;
	m = rhs.m; rhs.m = 0;
	p = std::move(rhs.p);
	return *this;
}

template <class T>
PermutationM<T> & PermutationM<T>::operator=(AMatrix<T> &&rhs) {
	PermutationM<T> &&rhs_ = dynamic_cast<PermutationM<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void PermutationM<T>::fill(const T val) {
	assert(false);
}

template <class T>
void PermutationM<T>::resize(const Int nnew) {
	assert(false);
}

template <class T>
void PermutationM<T>::solve(Vector<T> & rhssol) {
	Int size = rhssol.size();
	assert(size == n);
//https://blogs.msdn.microsoft.com/oldnewthing/20170102-00/?p=95095
	std::vector<Int> p_copy(p);
	Int curr, next;
	for (Int i = 0; i < n; i++) {
		curr = i;
		while (i != p_copy[curr]) {
			next = p_copy[curr];
			SWAP(rhssol[curr], rhssol[next]);
			p_copy[curr] = curr;
			curr = next;
		}
		p_copy[curr] = curr;
	}
}

template <class T>
void PermutationM<T>::print() const {
	cout << "Permutation matrix (permutation in cyclic form):" << endl;
	for (Int i = 0; i < n; i++)
		cout << p[i] << endl;
}

template <class T>
void PermutationM<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	f << "Permutation matrix (permutation in cyclic form):" << endl;
	for (Int i = 0; i < n; i++)
		f << p[i] << endl;
	f.close();
}

template <class T>
PermutationM<T> & PermutationM<T>::operator*=(const T c) {
	assert(false);
	return *this;
}

template <class T>
Vector<T> PermutationM<T>::operator*(const Vector<T> &vec) const {
	Int size = vec.size();
	assert(size == n);
	Vector<T> y(vec.size());
	for (Int i = 0; i < n; i++)
		y[p[i]] = vec[i];
	return y;
}

template <class T>
double PermutationM<T>::sizeGb() const {
	return p.size() * sizeof(Int)*bytes2Gbytes;
}