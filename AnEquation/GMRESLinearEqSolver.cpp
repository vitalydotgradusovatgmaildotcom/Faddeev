#include "GMRESLinearEqSolver.h"

template class GMRESLinearEqSolver<double>;
template class GMRESLinearEqSolver<Complex>;

template <class T>
GMRESLinearEqSolver<T>::GMRESLinearEqSolver(void) : ks(nullptr) { }

template <class T>
void GMRESLinearEqSolver<T>::calculateSolution(ALinearEq<T> &eq) {
	header("Start Arnoldi GMRES process ...");
	cout << "Restart    iter    rel. residue" << endl;
	
	//using initial guess x0 = eq.sol
	AnOperator<T> &op = *(eq.op);
	const AFunction<T> &rhs = *(eq.rhs);
	double resnorm, rhsnorm = rhs.coef.norm();
	unique_ptr<AFunction<T>> residue = eq.sol->clone();

	dimmax = min(op.getRank()+1, (Int)KRYLOV_DIM_MAX_GMRES);
	eq.getResidue(*residue);
	resnorm = residue->coef.norm();

	Int m, nrestart = 0;
	si.reserve(dimmax); ci.reserve(dimmax);
	Hmrot.resize(dimmax); gmrot.resize(dimmax);
	ks = new KrylovSubspace<T>(op, *residue, dimmax);
	Vector<T> ym;
	bool condRestart = resnorm/rhsnorm > reltol && \
				nrestart <= RESTART_NUM_MAX_GMRES;
	while(condRestart) {
		//Restarted GMRES
		ks->renew(*residue); m = ks->dim()-1;
		gmrot[0] = resnorm;
		while (resnorm/rhsnorm > reltol \
					&& m < dimmax-1) {
			ks->makeGramSchmidtStep(); m = ks->dim()-1;
			makeRot();
			resnorm = abs(gmrot[m]);
			cout << setw(11) << nrestart << setw(8) \
				<< m  << double135<double> \
											<< resnorm/rhsnorm;
			if (m % 3 == 0)
				cout << endl;
		}
		ym.resize(m);
		for (Int i = 0; i < m; i++)
			ym[i] = gmrot[i];
		Hmrot.solve(ym);
		//eq.sol = x0 !!!
		for (Int i = m-1; i >= 0; i--) {
			unique_ptr<AFunction<T>> tmp = (*ks)[i].clone();
			*tmp *= ym[i];
			*(eq.sol) += *tmp;
		}
		nrestart++;
		condRestart = resnorm/rhsnorm > reltol && \
				nrestart <= RESTART_NUM_MAX_GMRES;
		if (condRestart) {
			eq.getResidue(*residue);
			resnorm = residue->coef.norm();
			si.resize(0); ci.resize(0);
		}
	}
	cout << endl << endl;
	assert(nrestart <= RESTART_NUM_MAX_GMRES);
	eq.solved = true;
	Hmrot.resize(0); gmrot.resize(0);
	si.clear(); ci.clear();
	delete ks;
}

template <class T>
void GMRESLinearEqSolver<T>::makeRot() {
	Int m = ks->dim()-1;
	Vector<T> hmbarcol(m+1);
	const GenMatrix<T> & Hmbar = ks->getHbar();
	for (Int i = 0; i < m+1; i++)
		hmbarcol[i] = Hmbar[i][m-1];
	T h1, h2;
	for (Int i = 0; i < si.size(); i++) {
		h1 = hmbarcol[i]; h2 = hmbarcol[i+1];
		hmbarcol[i] = CONJ(ci[i])*h1 + si[i]*h2;
		hmbarcol[i+1] = -si[i]*h1 + ci[i]*h2;
	}
	double d = sqrt(norm(hmbarcol[m - 1]) + norm(hmbarcol[m]));
	si.push_back(hmbarcol[m] / d);
	ci.push_back(hmbarcol[m-1] / d);
	h1 = hmbarcol[m-1]; h2 = hmbarcol[m];
	hmbarcol[m-1] = CONJ(ci.back())*h1 + si.back()*h2;
	hmbarcol[m] = -si.back()*h1 + ci.back()*h2;
	for (Int i = 0; i < m; i++)
		Hmrot.set(i, m-1, hmbarcol[i]);
	T g = gmrot[m-1];
	gmrot[m-1] = CONJ(ci.back())*g;
	gmrot[m] = -si.back()*g;
}

template <class T>
GMRESLinearEqSolver<T>::~GMRESLinearEqSolver(void) { }
