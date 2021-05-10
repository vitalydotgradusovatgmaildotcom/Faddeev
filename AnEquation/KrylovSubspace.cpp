#include "KrylovSubspace.h"

template class KrylovSubspace<double>;
template class KrylovSubspace<Complex>;

template <class T>
//KrylovSubspace<T>::KrylovSubspace(const shared_ptr<const AnOperator<T>> &op, \
		const AFunction<T> &v0, const Int &nmax) : \
					nmax(nmax), op(op), Hbar(nmax+1, nmax) {
KrylovSubspace<T>::KrylovSubspace(AnOperator<T> &op, \
	const AFunction<T> &v0, const Int &nmax) : \
					nmax(nmax), op(op), Hbar(nmax, nmax-1) {
	V.reserve(nmax);
	V.push_back(v0.clone());
	n = 1;
	V[0]->coef.normalize();
	Hbar.fill(T());
}

template <class T>
void KrylovSubspace<T>::makeGramSchmidtStep() {
	assert(n <= nmax);

	V.push_back(V.back()->clone());
	op.times(**(V.end() - 2), *(V.back()));
	double avnorm = V.back()->coef.norm();
	for (Int i = 0; i < V.size()-1; i++) {
		Hbar[i][n-1] = V[i]->coef.scal(V.back()->coef);
		//Hbar[i][n - 1] = V.back()->coef.scal(V[i]->coef);
		unique_ptr<AFunction<T>> tmp = V[i]->clone();
		*tmp *= Hbar[i][n-1];
		*(V.back()) -= *tmp;
	}	

//reorthogonalisation
//http://web.mit.edu/ehliu/Public/sunhu/chapter5-6.pdf
	if ( V.back()->coef.norm() / avnorm <= REORTHOG_THRESH ) {
		T scal;
		for (Int i = 0; i < V.size()-1; i++) {
			scal = V[i]->coef.scal(V.back()->coef);
			//scal = V.back()->coef.scal(V[i]->coef);
			unique_ptr<AFunction<T>> tmp = V[i]->clone();
			*tmp *= scal;
			*(V.back()) -= *tmp;
			Hbar[i][n-1] += scal;
		}
	}

	Hbar[n][n-1] = V.back()->coef.norm();
	if (Hbar[n][n-1] != 0.0)
		*(V.back()) *= 1.0/Hbar[n][n-1];

	n++;
}

template <class T>
Int KrylovSubspace<T>::dim() {
	return n;
}

template <class T>
const AFunction<T> & KrylovSubspace<T>::operator[](const Int i) const {
	return *(V[i]);
}

template <class T>
const GenMatrix<T> KrylovSubspace<T>::getHbar() const {
	assert(n > 1);
	return GenMatrix<T>(Hbar, 0, n-1, 0, n-2);
}

template <class T>
void KrylovSubspace<T>::renew(const AFunction<T> &v0) {
	V.clear();
	V.push_back(v0.clone());
	n = 1;
	V[0]->coef.normalize();
}

template <class T>
void KrylovSubspace<T>::reduce(const Int newdim) {
	V.resize(newdim);
	n = newdim;
}

template <class T>
KrylovSubspace<T>::~KrylovSubspace(void) { }
