
#include "FaddeevConfigurator.h"
#include "System3Body.h"
#include "FMEq.h"
#include "FMEigenProb.h"
#include "RadEigenProb.h"
#include "BoundGrCalculator.h"
#include "GMRESLinearEqSolver.h"
#include "IRAMEigenSolver.h"
#include "Test.h"
#include "RadEq.h"
#include <chrono>

int main() {
	//mkl_set_num_threads(14);
	test();

	const FaddeevConfigurator &fc = \
		FaddeevConfigurator::getInstance();

	System3Body &s3b = *fc.s3b;
	Header("System");
	s3b.printSystem();
	std::vector<double> energies = fc.energies;
	switch (fc.problem) {
	case scatter: {
		//make and solve equation
		GMRESLinearEqSolver<Complex> solv;

		for (Int k = 0; k < energies.size(); k++) {
			//renew energy and solve
			Header("ENERGY = " + std::to_string(energies[k]));

			auto start = std::chrono::high_resolution_clock::now();
			FMEq fe(s3b, energies[k]);
			auto finish = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_op = finish - start;

			std::chrono::duration<double> elapsed_sv(0.0);
			std::chrono::duration<double> elapsed_post(0.0);
			while (fe.needsSolution()) {
				start = std::chrono::high_resolution_clock::now();
				solv.calculateSolution(fe);
				finish = std::chrono::high_resolution_clock::now();
				elapsed_sv += finish - start;

				start = std::chrono::high_resolution_clock::now();
				fe.postprocess();
				finish = std::chrono::high_resolution_clock::now();
				elapsed_post += finish - start;
			}

			start = std::chrono::high_resolution_clock::now();
			fe.calcCrossSect();
			finish = std::chrono::high_resolution_clock::now();
			elapsed_post += finish - start;


			header("Solution time");
			cout << "Build matrices, preconditioner and rhs time: " << \
				fixed << elapsed_op.count() << " sec" << endl;
			cout << "Solve equations time: " << \
				fixed << elapsed_sv.count() << " sec" << endl;
			cout << "Postprocessing time: " << \
				fixed << elapsed_post.count() << " sec" << endl;
			cout << endl;
		}
		break;
	}
	case bound: {
		Int nEig = fc.nEig;
		//make and solve eigen problem
		auto start = std::chrono::high_resolution_clock::now();
		FMEigenProb fp(s3b, energies[0]);
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_op = finish - start;

		Header("Calculate eigenvalues");
		IRAMEigenSolver<Complex> es;
		es.setRelTol(1.0e-12);
		start = std::chrono::high_resolution_clock::now();
		es.getEEV(fp, nEig, LarM, false);
		finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_sv = finish - start;
		header("Eigenvalues");
		Vector<Complex> ev(fp.getEval());
		for (Int k = 0; k < ev.size(); k++)
			ev[k] = energies[0] + zone / ev[k];
		sort(ev.v, ev.v + ev.size(), [&ev](Complex a, Complex b) { return a.real() < b.real(); });
		ev.print();
		//for (Int k = 0; k < ev.size(); k++)
		//	cout << k << " " << double2315<Complex> << ev[k] << endl;
		//cout << endl;
		header("Solution time");
		cout << "Build matrices and preconditioner time: " << \
			fixed << elapsed_op.count() << " sec" << endl;
		cout << "Find eigenvalues time: " << \
			fixed << elapsed_sv.count() << " sec" << endl;
		cout << endl;
		break;
	}
	case scatter_rad: {
		Pair pair = s3b.getPair(fc.numPair);
		//pair.alpha = 0;
		Int l = fc.J;
		for (Int k = 0; k < energies.size(); k++) {
			RadEq re(pair, l, energies[k]);
			re.writeSol();
		}
		break;
	}
	case bound_rad: {
		Pair pair = s3b.getPair(fc.numPair);
		//pair.alpha = 0;
		Int l = fc.J;
		RadEigenProb tbep(pair, l);
		//tbep.writeSolution();
		header("Eigenvalues");
		Vector<Complex> ev = tbep.getEval();
		for (Int k = 0; k < tbep.numberOfBound(); k++)
			cout << k << " " << double2315<Complex> << ev[k] << endl;
		cout << endl;
		break;
	}
	case bound_gr: {
		Pair pair = s3b.getPair(fc.numPair);
		//pair.alpha = 0;
		BoundGrCalculator bgc(pair);
	}
	}
}