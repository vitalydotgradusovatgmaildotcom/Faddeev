
#include "cph.h"
#include "GMRESLinearEqSolver.h"
#include "TestEq.h"
#include "MatrixOp.h"
#include "Collocation.h"
#include "HSpline5.h"
#include "HSpline3.h"
#include "Function.h"
#include "VectDiscr.h"
#include "VectFunction.h"
#include "MatrixProd.h"
#include "BlockDiagMatr.h"
#include "BlockMatr.h"
#include "UTMatrix.h"
#include "TriDiagM.h"
#include "PermutationM.h"
#include "IRAMEigenSolver.h"
#include "TestEP.h"
#include "TTVector.h"
#include "System3Body.h"
#include "SparseMatr.h"
#include "CWFCalculator.h"
#include "HybridBasis.h"
#include "wignerSymbols-cpp.h"
#include "SimplePoly.h"

void test() {

	//GMRESLinearEqSolver test
	
	/*
	TestEq eq;
	std::shared_ptr<TestOperator> tOp = dynamic_pointer_cast<TestOperator, AnOperator<TEST_OP_TYPE>>(eq.op);

	GMRESLinearEqSolver<TEST_OP_TYPE> solv;
	solv.setRelTol(1e-10);

	eq.sol->coef.fill(0.0);

	solv.calculateSolution(eq);
	Function<TEST_OP_TYPE, 1> resid(tOp->getDiscr());
	eq.getResidue(resid);
	cout << "Residue: " << resid.coef.norm() << endl;

	cout << "Solution:" << endl;
	eq.sol->coef.print();
	*/
//=============================
	//IRAMEigenSolver test
/*
#ifdef TEST_EP_SPARSE_BIG
	TestEP ep;
	
	IRAMEigenSolver<TEST_OP_SP_TYPE> solv;

	solv.getEEV(ep, 10, LarM, false);
	ep.eval.print();

#else

	TestEP ep;
	std::shared_ptr<TestOperator> tOp = dynamic_pointer_cast<TestOperator, AnOperator<TEST_OP_TYPE>>(ep.op);
	
	GenMatrix<TEST_OP_TYPE> a(tOp->matr.nrows(), tOp->matr.ncols());
	for (Int i = 0; i < a.nrows(); i++)
		for (Int j = 0; j < a.ncols(); j++) {
			a[i][j] = tOp->matr[i][j];
		}
	Vector<Complex> ev;
	GenMatrix<Complex> vm1, v;
	a.getEEV(ev, vm1, v, false, false);
	ev.print();

	IRAMEigenSolver<TEST_OP_TYPE> solv;
	solv.setRelTol(1e-10);

	solv.getEEV(ep, 2, LarR, false);
	ep.eval.print();

#endif
	*/
//=============================
	// Splines tests
	/*
	vector<double> vvv; vvv.resize(4);
	vvv[0] = 0.5; vvv[1] = 0.6; vvv[2] = 0.85; vvv[3] = 1.0;
	Grid gr(vvv);
	HSpline5<double> bas(gr);
	ofstream f; f.open("splines.dat");
	double h = 0.01;
	double x = 0.5;
	while (x <= 1.005) {
		f << x << "  " << bas.f(x, 6) << endl;
		x += h;
	}
	f << endl << endl;

	x = 0.5;
	while (x <= 1.005) {
		f << x << "  " << 10.0*bas.f(x, 7) << endl;
		x += h;
	}
	f << endl << endl;

	x = 0.5;
	while (x <= 1.005) {
		f << x << "  " << 200.0*bas.f(x, 8) << endl;
		x += h;
	}
	f << endl << endl;

	x = 0.5;
	while (x <= 1.005) {
		f << x << "  " << bas.f(x, 9) << endl;
		x += h;
	}
	f << endl << endl;

	x = 0.5;
	while (x <= 1.005) {
		f << x << "  " << 10.0*bas.f(x, 10) << endl;
		x += h;
	}
	f << endl << endl;

	x = 0.5;
	while (x <= 1.005) {
		f << x << "  " << 200.0*bas.f(x, 11) << endl;
		x += h;
	}
	f << endl << endl;


	f.close();
	*/

	/*
	vector<double> vvv; vvv.resize(4);
	vvv[0] = 0.5; vvv[1] = 0.6; vvv[2] = 0.85; vvv[3] = 1.0;
	Grid gr(vvv);

	array<unique_ptr<const ABasis>, 1> barr;
	//barr[0] = make_unique<const HSpline5<double>>(gr, Dir, Neu, 0.5, -1.0);
	barr[0] = make_unique<const HSpline5<double>>(gr, Mix, Mix, 0.5, -1.0);
	shared_ptr<const ADiscretization<double, 1>> discr = \
		make_shared<const Collocation<double, 1>>(barr);

	Function<double, 1> func(discr);
	for (Int i = 0; i < func.size(); i++)
		func[i] = 2.0;

	ofstream f; f.open("splines.dat");
	double h = 0.01;
	array<double, 1> x = { 0.45 };
	while (x[0] <= 1.05) {
		f << x[0] << "  " << discr->f(x, func) << endl;
		x[0] += h;
	}
	f << endl << endl;

	f.close();
	*/
//=============================
	//Discretization and Function (3d) test

	/*

	vector<double> vvv; vvv.resize(6);
	vvv[0] = -7.8; vvv[1] = -4.8; vvv[2] = -1.0;
	vvv[3] = 2.2; vvv[4] = 9.3; vvv[5] = 18.0;
	Grid gr(vvv);
	array<unique_ptr<const ABasis>, 3> barr;
	barr[0] = make_unique<const HSpline5<double>>(gr, Dir, Neu, 1.0, 1.0);
	vvv.resize(4); vvv[0] = -10.0; gr = Grid(vvv);
	barr[1] = make_unique<const HSpline5<double>>(gr, Neu, none, -1.0, 1.0);
	vvv.resize(5); vvv[1] = -2.0; vvv[4] = 5.0; gr = Grid(vvv);
	barr[2] = make_unique<const HSpline5<double>>(gr, none, none, -1.0, 1.0);

	vector<const shared_ptr<const ADiscretization<double, 3>>> discrs;
	discrs.reserve(3);

	shared_ptr<const ADiscretization<double, 3>> discr = \
		make_shared<const Collocation<double, 3>>(barr);
	discrs.push_back(discr);

	vvv.resize(2); vvv[0] = -5.0; gr = Grid(vvv);
	barr[0] = make_unique<const HSpline5<double>>(gr, Dir, Neu, 1.0, 1.0);
	barr[1] = make_unique<const HSpline5<double>>(gr, Neu, none, -1.0, 1.0);
	barr[2] = make_unique<const HSpline5<double>>(gr, none, none, -1.0, 1.0);

	discr = \
		make_shared<const Collocation<double, 3>>(barr);
	discrs.push_back(discr);

	vvv.resize(3); vvv[0] = -3.0; gr = Grid(vvv);
	barr[0] = make_unique<const HSpline5<double>>(gr, Dir, Neu, 1.0, 1.0);
	vvv.resize(2); vvv[0] = -4.0; gr = Grid(vvv);
	barr[1] = make_unique<const HSpline5<double>>(gr, Neu, none, -1.0, 1.0);
	barr[2] = make_unique<const HSpline5<double>>(gr, none, none, -1.0, 1.0);

	discr = \
		make_shared<const Collocation<double, 3>>(barr);
	discrs.push_back(discr);

	vvv.resize(2); vvv[0] = -6.0; gr = Grid(vvv);
	barr[0] = make_unique<const HSpline5<double>>(gr, Dir, Neu, 1.0, 1.0);
	barr[1] = make_unique<const HSpline5<double>>(gr, Neu, none, -1.0, 1.0);
	barr[2] = make_unique<const HSpline5<double>>(gr, none, none, -1.0, 1.0);

	discr = \
		make_shared<const Collocation<double, 3>>(barr);
	discrs.push_back(discr);

	shared_ptr<const VectDiscr<double, 3>> vd = make_shared<const VectDiscr<double, 3>>(discrs);
	cout << "Total size: " << vd->getN() << ", nFunc: " \
		<< vd->getNFunc() << endl;
	Int iFunc; array<Int, 3> inds;
	vd->getInds(3089, iFunc, inds);
	cout << "iFunc = " << iFunc << ", inds = " << inds[0] \
		<< "  " << inds[1] << "  " << inds[2] << endl;
	cout << "i = " << vd->getRaw(iFunc, inds) << endl;
	VectFunction<double, 3> func(vd);
	func[0] = 2.0;
	func[1] = -3.0;
	func *= 10.0;
	for (Int k = 0; k < func.size(); k++)
		func[k] = 1.0;
	Int jjj = 0;
	array<double, 3> xxx = {-5.0, -4.0, -3.0};
	func = func + func*9.0;//func += func*9.0;//func *= 10.0;
	cout << func.getF(3)(xxx) << endl;
	*/

	/*
	vector<double> vvv; vvv.resize(6);
	vvv[0] = -7.8; vvv[1] = -4.8; vvv[2] = -1.0;
	vvv[3] = 2.2; vvv[4] = 9.3; vvv[5] = 18.0;
	Grid gr(vvv);
	
	array<unique_ptr<const ABasis>, 3> barr;
	barr[0] = make_unique<const HSpline5<double>>(gr, Dir, Neu, 1.0, 1.0);
	barr[1] = make_unique<const HSpline5<double>>(gr, Neu, none, -1.0, 1.0);
	barr[2] = make_unique<const HSpline5<double>>(gr, none, none, -1.0, 1.0);
	//barr[0] = make_unique<const HSpline5<double>>(gr, none, none, 1.0, 1.0);
	//barr[1] = make_unique<const HSpline5<double>>(gr, none, none, -1.0, 1.0);

	shared_ptr<const ADiscretization<double, 3>> d = \
		make_shared<const Collocation<double, 3>>(barr);

	ofstream f; f.open("OUTPUT.dat");
	double a, b, h; a = -8.79; //-10.0;
	b = -4.0;//20.0;
	h = 1.0;

	array<double, 3> x = {a, a, a};
	Function<double, 3> func(d);
	for (Int i = 0; i < func.size(); i++)
		func[i] = 0.0;
	//func[0] = 2.0;
	//func[1] = -10.0; func[2] = -1.0;
	//func[3] = -1.0; func[4] = 2.0;
	func[612] = 2.0;
	//func[126] = 3.0;
	//func[func.size()-3] = 10.0;
	//func[func.size()-20] = 1.0;
	//func[func.size()-37] = -2.0;
	while (x[0] <= b) {
		x[1] = a;
		while (x[1] <= b) {
			x[2] = a;
			while (x[2] <= b) {
				f << double135<double> << x[0] << "\t" \
					<< double135<double> << x[1] << "\t" \
					<< double135<double> << x[2] << "\t" \
					<< double135<double> << d->f(x, func) << endl;
				x[2] += h;
			}
			x[1] += h;
		}
		f << endl;
		x[0] += h;
	}

	f << endl << endl;
	f.close();
	*/
	
//==================================================
//TensorProd test
	
	/*
	auto m1 = make_shared<GenMatrix<Complex>>(2);
	auto m2 = make_shared<DiagMatrix<Complex>>(2);
	auto m3 = make_shared<GenMatrix<Complex>>(2);
	(*m1)[0][0] = Complex(2, 1); (*m1)[0][1] = Complex(3, 3);
	(*m1)[1][0] = Complex(1, -1); (*m1)[1][1] = Complex(4, 0);
	(*m2)[0] = Complex(1, 1); (*m2)[1] = Complex(2, -1);
	(*m3)[0][0] = Complex(-1, 2); (*m3)[0][1] = Complex(3, 2);
	(*m3)[1][0] = Complex(1, 0); (*m3)[1][1] = Complex(1, -1);

	std::array<shared_ptr<AMatrix<Complex>>, 3> matrs;
	matrs[0] = m1;
	matrs[1] = m2;
	matrs[2] = m3;
	TensorProd<Complex, 3> m(matrs);

	Vector<Complex> v(8);
	v[0] = Complex(-4, 1); v[1] = -3; v[2] = -2; v[3] = -1;
	v[4] = 0; v[5] = 1; v[6] = 2; v[7] = 3;

	m.solve(v);
	v.print();
	*/

	/*
	auto m1 = make_shared<GenMatrix<double>>(3, 2);
	auto m2 = make_shared<GenMatrix<double>>(4, 3);
	auto m3 = make_shared<GenMatrix<double>>(2, 3);
	(*m1)[0][0] = 1; (*m1)[0][1] = 2;
	(*m1)[1][0] = 3; (*m1)[1][1] = 4;
	(*m1)[2][0] = 2; (*m1)[2][1] = 1;
	
	(*m2)[0][0] = 2; (*m2)[0][1] = 3; (*m2)[0][2] = -5;
	(*m2)[1][0] = 1; (*m2)[1][1] = 4; (*m2)[1][2] = 3;
	(*m2)[2][0] = -1; (*m2)[2][1] = 2; (*m2)[2][2] = -3;
	(*m2)[3][0] = 6; (*m2)[3][1] = 1; (*m2)[3][2] = 1;

	(*m3)[0][0] = -1; (*m3)[0][1] = 3; (*m3)[0][2] = -4;
	(*m3)[1][0] = 1; (*m3)[1][1] = 1; (*m3)[1][2] = 2;
	
	std::array<shared_ptr<AMatrix<double>>, 3> matrs;
	matrs[0] = m1;
	matrs[1] = m2;
	matrs[2] = m3;
	TensorProd<double, 3> m(matrs);
	GenMatrix<double> mg = m.toGenMatrix();
	mg.print();
	Vector<double> v(18);
	v[0] = -4; v[1] = -3; v[2] = -2; v[3] = -1; 
	v[4] = 0; v[5] = 1; v[6] = 2; v[7] = 3;
	v[8] = 0; v[9] = 4; v[10] = 3; v[11] = 2;
	v[12] = -1; v[13] = -2; v[14] = -3; v[15] = -4;
	v[16] = -1; v[17] = 1;

	v = m * v;
	v.print();
	*/

//===================================================
	//MatrixProd test

	/*	
	auto m1 = make_shared<GenMatrix<Complex>>(7);
	auto m2 = make_shared<DiagMatrix<Complex>>(7);
	auto m3 = make_shared<GenMatrix<Complex>>(7);

	(*m1)[0][0] = Complex(10, 11); (*m1)[0][1] = Complex(0.54, 0.95); (*m1)[0][2] = Complex(0.80, 0.69); (*m1)[0][3] = Complex(0.03, 0.84); (*m1)[0][4] = Complex(0.65, 0.25); (*m1)[0][5] = Complex(0.82, 0.91); (*m1)[0][6] = Complex(0.76, 0.05);
	(*m1)[1][0] = Complex(0.90,  0.27); (*m1)[1][1] = Complex(12, 13); (*m1)[1][2] = Complex(0.14, 0.89); (*m1)[1][3] = Complex(0.84, 0.25); (*m1)[1][4] = Complex(0.17, 0.61); (*m1)[1][5] = Complex(0.69, 0.28); (*m1)[1][6] = Complex(0.79, 0.53);
	(*m1)[2][0] = Complex(0.12, 0.67); (*m1)[2][1] = Complex(0.96, 0.58); (*m1)[2][2] = Complex(13, 12); (*m1)[2][3] = Complex(0.93, 0.81); (*m1)[2][4] = Complex(0.70, 0.47); (*m1)[2][5] = Complex(0.31, 0.75); (*m1)[2][6] = Complex(0.18, 0.77);
    (*m1)[3][0] = Complex(0.91, 0.65); (*m1)[3][1] = Complex(0.15, 0.22); (*m1)[3][2] = Complex(0.91, 0.54); (*m1)[3][3] = Complex(10, 9); (*m1)[3][4] = Complex(0.03, 0.35); (*m1)[3][5] = Complex(0.95, 0.75); (*m1)[3][6] = Complex(0.48, 0.93);
    (*m1)[4][0] = Complex(0.63, 0.16); (*m1)[4][1] = Complex(0.97, 0.75); (*m1)[4][2] = Complex(0.79, 0.13); (*m1)[4][3] = Complex(0.75, 0.92); (*m1)[4][4] = Complex(14, 15); (*m1)[4][5] = Complex(0.03, 0.38); (*m1)[4][6] = Complex(0.44, 0.12);
    (*m1)[5][0] = Complex(0.09, 0.11); (*m1)[5][1] = Complex(0.95, 0.25); (*m1)[5][2] = Complex(0.95, 0.14); (*m1)[5][3] = Complex(0.74, 0.34); (*m1)[5][4] = Complex(0.04, 0.58); (*m1)[5][5] = Complex(11, 11); (*m1)[5][6] = Complex(0.64, 0.56);
    (*m1)[6][0] = Complex(0.27, 0.49); (*m1)[6][1] = Complex(0.48, 0.50); (*m1)[6][2] = Complex(0.65, 0.25); (*m1)[6][3] = Complex(0.39, 0.19); (*m1)[6][4] = Complex(0.09, 0.54); (*m1)[6][5] = Complex(0.38, 0.07); (*m1)[6][6] = Complex(15, 16);

	/*
	A1 = [10 + 11i 0.54 + 0.95i 0.80 + 0.69i 0.03 + 0.84i 0.65 + 0.25i 0.82 + 0.91i 0.76 + 0.05i;
    0.90 + 0.27i 12 + 13i 0.14 + 0.89i 0.84 + 0.25i 0.17 + 0.61i 0.69 + 0.28i 0.79 + 0.53i;
    0.12 + 0.67i 0.96 + 0.58i 13 + 12i 0.93 + 0.81i 0.70 + 0.47i 0.31 + 0.75i 0.18 + 0.77i;
    0.91 + 0.65i 0.15 + 0.22i 0.91 + 0.54i 10 + 9i 0.03 + 0.35i 0.95 + 0.75i 0.48 + 0.93i;
    0.63 + 0.16i 0.97 + 0.75i 0.79 + 0.13i 0.75 + 0.92i 14 + 15i 0.03 + 0.38i 0.44 + 0.12i;
    0.09 + 0.11i 0.95 + 0.25i 0.95 + 0.14i 0.74 + 0.34i 0.04 + 0.58i 11 + 11i 0.64 + 0.56i;
    0.27 + 0.49i 0.48 + 0.50i 0.65 + 0.25i 0.39 + 0.19i 0.09 + 0.54i 0.38 + 0.07i 15 + 16i];*/
	/*
	(*m2)[0] = Complex(1.0, 7.0); (*m2)[1] = Complex(2.0, 6.0); (*m2)[2] = Complex(3.0, 5.0); (*m2)[3] = Complex(4.0, 4.0); (*m2)[4] = Complex(5.0, 3.0); (*m2)[5] = Complex(6.0, 2.0); (*m2)[6] = Complex(7.0, 1.0);

	//A2 = diag([1+7i 2+6i 3+5i 4+4i 5+3i 6+2i 7+i]);

	(*m3)[0][0] = 10; (*m3)[0][1] = 0.546881519204984; (*m3)[0][2] = 0.800280468888800; (*m3)[0][3] = 0.0357116785741896; (*m3)[0][4] = 0.655477890177557; (*m3)[0][5] = 0.823457828327293; (*m3)[0][6] = 0.765516788149002;
	(*m3)[1][0] = 0.905791937075619; (*m3)[1][1] = 11; (*m3)[1][2] = 0.141886338627215; (*m3)[1][3] = 0.849129305868777; (*m3)[1][4] = 0.171186687811562; (*m3)[1][5] = 0.694828622975817; (*m3)[1][6] = 0.795199901137063;
	(*m3)[2][0] = 0.126986816293506; (*m3)[2][1] = 0.964888535199277; (*m3)[2][2] = 12; (*m3)[2][3] = 0.933993247757551; (*m3)[2][4] = 0.706046088019609; (*m3)[2][5] = 0.317099480060861; (*m3)[2][6] = 0.186872604554379;
    (*m3)[3][0] = 0.913375856139019; (*m3)[3][1] = 0.157613081677548; (*m3)[3][2] = 0.915735525189067; (*m3)[3][3] = 13; (*m3)[3][4] = 0.0318328463774207; (*m3)[3][5] = 0.950222048838355; (*m3)[3][6] = 0.489764395788231;
    (*m3)[4][0] = 0.632359246225410; (*m3)[4][1] = 0.970592781760616; (*m3)[4][2] = 0.792207329559554; (*m3)[4][3] = 0.757740130578333; (*m3)[4][4] = 14; (*m3)[4][5] = 0.0344460805029088; (*m3)[4][6] = 0.445586200710900;
    (*m3)[5][0] = 0.0975404049994095; (*m3)[5][1] = 0.957166948242946; (*m3)[5][2] = 0.959492426392903; (*m3)[5][3] = 0.743132468124916; (*m3)[5][4] = 0.0461713906311539; (*m3)[5][5] = 15; (*m3)[5][6] = 0.646313010111265;
    (*m3)[6][0] = 0.278498218867048; (*m3)[6][1] = 0.485375648722841; (*m3)[6][2] = 0.655740699156587; (*m3)[6][3] = 0.392227019534168; (*m3)[6][4] = 0.0971317812358475; (*m3)[6][5] = 0.381558457093008; (*m3)[6][6] = 16;

	/*
	A3 = [10 0.546881519204984 0.800280468888800 0.0357116785741896 0.655477890177557 0.823457828327293 0.765516788149002;
    0.905791937075619 11 0.141886338627215 0.849129305868777 0.171186687811562 0.694828622975817 0.795199901137063;
    0.126986816293506 0.964888535199277 12 0.933993247757551 0.706046088019609 0.317099480060861 0.186872604554379;
    0.913375856139019 0.157613081677548 0.915735525189067 13 0.0318328463774207 0.950222048838355 0.489764395788231;
    0.632359246225410 0.970592781760616 0.792207329559554 0.757740130578333 14 0.0344460805029088 0.445586200710900;
    0.0975404049994095 0.957166948242946 0.959492426392903 0.743132468124916 0.0461713906311539 15 0.646313010111265;
    0.278498218867048 0.485375648722841 0.655740699156587 0.392227019534168 0.0971317812358475 0.381558457093008 16];
	*/

	/*Vector<Complex> v(7);
	v[0] = Complex(0.01, 0.60); v[1] = Complex(0.33, 0.26); v[2] = Complex(0.16, 0.65); v[3] = Complex(0.79, 0.68); v[4] = Complex(0.31, 0.74); v[5] = Complex(0.52, 0.45); v[6] = Complex(0.16, 0.08);  
	//f = [0.01 + 0.60i;0.33 + 0.26i;0.16 + 0.65i;0.79 + 0.68i;0.31 + 0.74i;0.52 + 0.45i;0.16 + 0.08i];

	m1->print(); m2->print(); m3->print();

	std::vector<std::shared_ptr<AMatrix<Complex>>> matrs;
	matrs.reserve(3);
	matrs.push_back(m1); matrs.push_back(m2); matrs.push_back(m3);
	std::vector<bool> inverse;
	inverse.reserve(3);
	inverse.push_back(true); inverse.push_back(true); inverse.push_back(true);
	MatrixProd<Complex> mp(matrs, inverse);

	mp *= Complex(2.0, 0.0);
	v = mp * v;
	//mp.solve(v);
	v.print();
	mp.print();*/

//===================================================
	//genEEV and simultDiagonalize test
	/*GenMatrix<Complex> a(7), b(7);

	a[0][0] = Complex(10, 11); a[0][1] = Complex(0.54, 0.95); a[0][2] = Complex(0.80, 0.69); a[0][3] = Complex(0.03, 0.84); a[0][4] = Complex(0.65, 0.25); a[0][5] = Complex(0.82, 0.91); a[0][6] = Complex(0.76, 0.05);
	a[1][0] = Complex(0.90,  0.27); a[1][1] = Complex(12, 13); a[1][2] = Complex(0.14, 0.89); a[1][3] = Complex(0.84, 0.25); a[1][4] = Complex(0.17, 0.61); a[1][5] = Complex(0.69, 0.28); a[1][6] = Complex(0.79, 0.53);
	a[2][0] = Complex(0.12, 0.67); a[2][1] = Complex(0.96, 0.58); a[2][2] = Complex(13, 12); a[2][3] = Complex(0.93, 0.81); a[2][4] = Complex(0.70, 0.47); a[2][5] = Complex(0.31, 0.75); a[2][6] = Complex(0.18, 0.77);
    a[3][0] = Complex(0.91, 0.65); a[3][1] = Complex(0.15, 0.22); a[3][2] = Complex(0.91, 0.54); a[3][3] = Complex(10, 9); a[3][4] = Complex(0.03, 0.35); a[3][5] = Complex(0.95, 0.75); a[3][6] = Complex(0.48, 0.93);
    a[4][0] = Complex(0.63, 0.16); a[4][1] = Complex(0.97, 0.75); a[4][2] = Complex(0.79, 0.13); a[4][3] = Complex(0.75, 0.92); a[4][4] = Complex(14, 15); a[4][5] = Complex(0.03, 0.38); a[4][6] = Complex(0.44, 0.12);
    a[5][0] = Complex(0.09, 0.11); a[5][1] = Complex(0.95, 0.25); a[5][2] = Complex(0.95, 0.14); a[5][3] = Complex(0.74, 0.34); a[5][4] = Complex(0.04, 0.58); a[5][5] = Complex(11, 11); a[5][6] = Complex(0.64, 0.56);
    a[6][0] = Complex(0.27, 0.49); a[6][1] = Complex(0.48, 0.50); a[6][2] = Complex(0.65, 0.25); a[6][3] = Complex(0.39, 0.19); a[6][4] = Complex(0.09, 0.54); a[6][5] = Complex(0.38, 0.07); a[6][6] = Complex(15, 16);
	*/

	/*
	GenMatrix<double> a(7), b(7);

	a[0][0] = 10; a[0][1] = 0.54; a[0][2] = 0.80; a[0][3] = 0.03; a[0][4] = 0.65; a[0][5] = 0.82; a[0][6] = 0.76;
	a[1][0] = 0.90; a[1][1] = 12; a[1][2] = 0.14; a[1][3] = 0.84; a[1][4] = 0.17; a[1][5] = 0.69; a[1][6] = 0.79;
	a[2][0] = 0.12; a[2][1] = 0.96; a[2][2] = 13; a[2][3] = 0.93; a[2][4] = 0.70; a[2][5] = 0.31; a[2][6] = 0.18;
    a[3][0] = 0.91; a[3][1] = 0.15; a[3][2] = 0.91; a[3][3] = 10; a[3][4] = 0.03; a[3][5] = 0.95; a[3][6] = 0.48;
    a[4][0] = 0.63; a[4][1] = 0.97; a[4][2] = 0.79; a[4][3] = 0.75; a[4][4] = 14; a[4][5] = 0.03; a[4][6] = 0.44;
    a[5][0] = 0.09; a[5][1] = 0.95; a[5][2] = 0.95; a[5][3] = 0.74; a[5][4] = 0.04; a[5][5] = 11; a[5][6] = 0.64;
    a[6][0] = 0.27; a[6][1] = 0.48; a[6][2] = 0.65; a[6][3] = 0.39; a[6][4] = 0.09; a[6][5] = 0.38; a[6][6] = 15;
	

	b[0][0] = 10; b[0][1] = 0.546881519204984; b[0][2] = 0.800280468888800; b[0][3] = 0.0357116785741896; b[0][4] = 0.655477890177557; b[0][5] = 0.823457828327293; b[0][6] = 0.765516788149002;
	b[1][0] = 0.905791937075619; b[1][1] = 11; b[1][2] = 0.141886338627215; b[1][3] = 0.849129305868777; b[1][4] = 0.171186687811562; b[1][5] = 0.694828622975817; b[1][6] = 0.795199901137063;
	b[2][0] = 0.126986816293506; b[2][1] = 0.964888535199277; b[2][2] = 12; b[2][3] = 0.933993247757551; b[2][4] = 0.706046088019609; b[2][5] = 0.317099480060861; b[2][6] = 0.186872604554379;
    b[3][0] = 0.913375856139019; b[3][1] = 0.157613081677548; b[3][2] = 0.915735525189067; b[3][3] = 13; b[3][4] = 0.0318328463774207; b[3][5] = 0.950222048838355; b[3][6] = 0.489764395788231;
    b[4][0] = 0.632359246225410; b[4][1] = 0.970592781760616; b[4][2] = 0.792207329559554; b[4][3] = 0.757740130578333; b[4][4] = 14; b[4][5] = 0.0344460805029088; b[4][6] = 0.445586200710900;
    b[5][0] = 0.0975404049994095; b[5][1] = 0.957166948242946; b[5][2] = 0.959492426392903; b[5][3] = 0.743132468124916; b[5][4] = 0.0461713906311539; b[5][5] = 15; b[5][6] = 0.646313010111265;
    b[6][0] = 0.278498218867048; b[6][1] = 0.485375648722841; b[6][2] = 0.655740699156587; b[6][3] = 0.392227019534168; b[6][4] = 0.0971317812358475; b[6][5] = 0.381558457093008; b[6][6] = 16;
	
	GenMatrix<Complex> wbar(0), w(0);
	DiagMatrix<Complex> lambda(0);

	GenMatrix<Complex> a_copy(a.nrows()), b_copy(b.nrows());
	for (Int i = 0; i < a.nrows(); i++)
		for (Int j = 0; j < a.ncols(); j++) {
			a_copy[i][j] = a[i][j];
			b_copy[i][j] = b[i][j];
		}

	simultDiagonalize(a, b, wbar, lambda, w);
	(wbar * a_copy * w).print();
	lambda.print();
	(wbar * b_copy * w).print();

	*/
	/*genEEV(a, b, wl, ev, wr, false, false);
	ev.print();*/
//==========================================
	//BlockDiagMatr test
	/*
	shared_ptr<UTMatrix<double>> ut = \
		make_shared<UTMatrix<double>>(2);
	ut->set(0, 0, 4.0); ut->set(0, 1, 5.0);
	ut->set(1, 1, 6.0); 
	shared_ptr<DiagMatrix<double>> d = \
		make_shared<DiagMatrix<double>>(3);
	(*d)[0] = 1.0; (*d)[1] = 2.0; (*d)[2] = 3.0;
	shared_ptr<GenMatrix<double>> gm = \
		make_shared<GenMatrix<double>>(2);
	(*gm)[0][0] = 7.0; (*gm)[0][1] = 8.0;
	(*gm)[1][0] = 9.0; (*gm)[1][1] = 10.0;
	std::vector<shared_ptr<AMatrix<double>>> blocks;
	blocks.reserve(3);
	blocks.push_back(ut); blocks.push_back(d);
	blocks.push_back(gm);
	vector<bool> invs; invs.resize(3);
	invs[0] = true; invs[1] = false; invs[2] = true;
	BlockDiagMatr<double> bdm(blocks, invs);
	bdm.print();

	Vector<double> vect(7);
	vect[0] = 1.0; vect[1] = 2.0; vect[2] = 3.0; vect[3] = 4.0;
	vect[4] = 5.0; vect[5] = 6.0; vect[6] = 7.0;
	vect.print();

	vect = bdm * vect;
	//bdm.solve(vect);
	vect.print();*/
//==========================================
	//TriDiagM test
	//A = diag(1:9, 1) + diag(10:19) + diag(20:28, -1)
	//v = (1:10)'
	/*TriDiagM<double> tdm(10);
	Vector<double> v(10);
	tdm.set(0, 0, 10.0);
	tdm.set(0, 1, 1.0);
	v[0] = 1.0;
	for (Int i = 1; i < 9; i++) {
		tdm.set(i, i, 10.0+i);
		tdm.set(i, i+1, 1.0+i);
		tdm.set(i, i-1, 19.0+i);
		v[i] = i+1.0;
	}
	tdm.set(9, 9, 19.0);
	tdm.set(9, 8, 28.0);
	v[9] = 10.0;*/

	//A = diag(1:9, 1) + diag(10:19) + diag(20:28, -1) ...
	//    +i*( diag(-9:-1, 1) + diag(-19:-10) + diag(-28:-20, -1) )
	//v = (1:10)'+i*(11:20)'
	
	/*
	TriDiagM<Complex> tdm(10);
	Vector<Complex> v(10);
	tdm.set(0, 0, Complex(10.0, -19.0));
	tdm.set(0, 1, Complex(1.0, -9.0));
	v[0] = Complex(1.0, 11.0);
	for (Int i = 1; i < 9; i++) {
		tdm.set(i, i, Complex(10.0+i, -19.0+i));
		tdm.set(i, i+1, Complex(1.0+i, -9.0+i) );
		tdm.set(i, i-1, Complex(19.0+i, -29.0+i));
		v[i] = Complex(i+1.0, i+11.0);
	}
	tdm.set(9, 9, Complex(19.0, -10.0));
	tdm.set(9, 8, Complex(28.0, -20.0));
	v[9] = Complex(10.0, 20.0);

	tdm.print();
	(tdm*v).print();
	//tdm.solve(v);
	//tdm.solve(v);
	//tdm.solve(v);
	//v.print();
	*/
//==========================================
	//Permutation matrix test
	/*
	std::vector<Int> perm;
	perm.reserve(8);
	perm.push_back(3); perm.push_back(1); perm.push_back(4); 
	perm.push_back(5); perm.push_back(2); perm.push_back(6); 
	perm.push_back(0); perm.push_back(7);
	PermutationM<double> pm(std::move(perm));
	//(1 2 3 4 5 6 7 8)
	//(4 2 5 6 3 7 1 8)
	//inverse:
	//(1 2 3 4 5 6 7 8)
	//(7 2 5 1 3 4 6 8)
	Vector<double> v(8);
	v[0] = 1.0; v[1] = 2.0; v[2] = 3.0; 
	v[3] = 4.0; v[4] = 5.0; v[5] = 6.0; 
	v[6] = 7.0; v[7] = 8.0;

	(pm*v).print();
	pm.solve(v);
	v.print();*/
//==========================================
//Wigner d small test
	/*
	Int j = 1;
	Int m = 1;
	Int mbar = 0;
	double a = 0.0; double b = 3.14;
	const Int n = 100;
	double h = (b - a) / (n - 1);
	double beta[n];
	double summa = 0.0;
	double wds;
	for (Int i = 0; i < n; i++) {
		beta[i] = a + i * h;
		wds = wignerDSmall(j, m, mbar, beta[i]);
		cout << beta[i] << " " << wds << endl;
		summa += wds * wds*sin(beta[i]);
	}
	summa *= h;
	cout << "(j+1/2) Int d cos D^2 = " << summa * (j + 0.5) << endl;
	*/
	//wigner D test
	/*
	Int j = 5;
	Int m = 4;
	Int mbar = 3;
	Int tautilde = 1;
	double x, y, z, phi, theta, vphi;
	x = 2.0; y = 1.5; z = 0.6;

	phi = PI / 4; theta = PI / 3; vphi = 2 * PI / 5; //z = -0.6;
	//phi = PI / 4; theta = PI / 3; vphi = 9 * PI / 10; z = -0.6;
	//phi = PI / 4; theta = PI / 3; vphi = 7 * PI / 5; z = -0.6;
	//phi = PI / 4; theta = PI / 3; vphi = 19 * PI / 10; z = -0.6;
	//phi = 3 * PI / 4; theta = PI / 3; vphi = 2 * PI / 5; z = -0.6;
	//phi = 3 * PI / 4; theta = PI / 3; vphi = 9 * PI / 10; z = -0.6;
	//phi = 3 * PI / 4; theta = PI / 3; vphi = 7 * PI / 5; z = -0.6;
	//phi = 3 * PI / 4; theta = PI / 3; vphi = 19 * PI / 10; z = -0.6;
	//phi = 5 * PI / 4; theta = PI / 3; vphi = 2 * PI / 5; z = -0.6;
	//phi = 5 * PI / 4; theta = PI / 3; vphi = 9 * PI / 10; z = -0.6;
	//phi = 5 * PI / 4; theta = PI / 3; vphi = 7 * PI / 5; z = -0.6;
	//phi = 5 * PI / 4; theta = PI / 3; vphi = 19 * PI / 10; z = -0.6;
	//phi = 7 * PI / 4; theta = PI / 3; vphi = 2 * PI / 5; z = -0.6;
	//phi = 7 * PI / 4; theta = PI / 3; vphi = 9 * PI / 10; z = -0.6;
	//phi = 7 * PI / 4; theta = PI / 3; vphi = 7 * PI / 5; z = -0.6;
	//phi = 7 * PI / 4; theta = PI / 3; vphi = 19 * PI / 10; z = -0.6;

	//phi = PI / 4; theta = 5 * PI / 6; vphi = 2 * PI / 5; z = -0.6;
	//phi = PI / 4; theta = 5 * PI / 6; vphi = 9 * PI / 10; z = -0.6;
	//phi = PI / 4; theta = 5 * PI / 6; vphi = 7 * PI / 5; z = -0.6;
	//phi = PI / 4; theta = 5 * PI / 6; vphi = 19 * PI / 10; z = -0.6;
	//phi = 3 * PI / 4; theta = 5 * PI / 6; vphi = 2 * PI / 5; z = -0.6;
	//phi = 3 * PI / 4; theta = 5 * PI / 6; vphi = 9 * PI / 10; z = -0.6;
	//phi = 3 * PI / 4; theta = 5 * PI / 6; vphi = 7 * PI / 5; z = -0.6;
	//phi = 3 * PI / 4; theta = 5 * PI / 6; vphi = 19 * PI / 10; z = -0.6;
	//phi = 5 * PI / 4; theta = 5 * PI / 6; vphi = 2 * PI / 5; z = -0.6;
	//phi = 5 * PI / 4; theta = 5 * PI / 6; vphi = 9 * PI / 10; z = -0.6;
	//phi = 5 * PI / 4; theta = 5 * PI / 6; vphi = 7 * PI / 5; z = -0.6;
	//phi = 5 * PI / 4; theta = 5 * PI / 6; vphi = 19 * PI / 10; z = -0.6;
	//phi = 7 * PI / 4; theta = 5 * PI / 6; vphi = 2 * PI / 5; z = -0.6;
	//phi = 7 * PI / 4; theta = 5 * PI / 6; vphi = 9 * PI / 10; z = -0.6;
	//phi = 7 * PI / 4; theta = 5 * PI / 6; vphi = 7 * PI / 5; z = -0.6;
	//phi = 7 * PI / 4; theta = 5 * PI / 6; vphi = 19 * PI / 10; z = -0.6;

	//Int alpha = 0; Int beta = 1;
	Int alpha = 1; Int beta = 0;
	
	double xrot, yrot, zrot, wba, phirot, throt, vphirot;
	//System3Body sys("He", {1.0, 1.0, 10.0}, \
					{-1.0, -1.0, 2.0}, 1, 0, -1);
	System3Body sys("He", { 1.0, 1.0, 1000000000.0 }, \
					{-1.0, -1.0, 2.0}, 1, 0, -1);
	jacobiCoo xyza, xyzb;
	xyza.x = x; xyza.y = y; xyza.z = z;
	sys.rotReduced(xyza, alpha, xyzb, beta);
	xrot = xyzb.x; yrot = xyzb.y; zrot = xyzb.z;

	GenMatrix<double> R(3);
	R[0][0] = cos(phi)*cos(theta)*cos(vphi) - sin(phi)*sin(vphi);
	R[0][1] = -cos(phi)*cos(theta)*sin(vphi) - sin(phi)*cos(vphi);
	R[0][2] = cos(phi)*sin(theta);
	R[1][0] = sin(phi)*cos(theta)*cos(vphi)+cos(phi)*sin(vphi);
	R[1][1] = -sin(phi)*cos(theta)*sin(vphi)+cos(phi)*cos(vphi);
	R[1][2] = sin(phi)*sin(theta);
	R[2][0] = -sin(theta)*cos(vphi);
	R[2][1] = sin(theta)*sin(vphi);
	R[2][2] = cos(theta);
	Vector<double> tmp(3);
	tmp[0] = 0; tmp[1] = 0; tmp[2] = 1;
	Vector<double> yvec(R * tmp); yvec *= y;
	tmp[0] = sin(acos(z)); tmp[1] = 0; tmp[2] = z;
	Vector<double> xvec(R * tmp); xvec *= x;
	Vector<double> xvecrot, yvecrot;
	xvecrot = xvec * sys.c[beta][alpha] + yvec * sys.s[beta][alpha];
	yvecrot = yvec * sys.c[beta][alpha]  - xvec * sys.s[beta][alpha];
	
	assert( !(yvecrot[0] == 0.0 && yvecrot[1] == 0.0) );
	if (yvecrot[0] == 0.0) {
		if (yvecrot[1] > 0)
			phirot = PI / 2;
		else //yvecrot[1] > 0
			phirot = 3 * PI / 2;
	} else { //yvecrot[0] != 0.0
		phirot = atan(yvecrot[1] / yvecrot[0]);
		if (yvecrot[0] > 0.0 && yvecrot[1] < 0.0)
			phirot += 2*PI;
		if (yvecrot[0] < 0.0)
			phirot += PI;
	}
	throt = acos(yvecrot[2] / yvecrot.norm());
	//vphirot = \
		acos((-yvecrot[2] / yvecrot.norm() + cos(throt)*zrot) / \
		(sin(throt)*sin(acos(zrot))));
	GenMatrix<double> Rrot2(3), Rrot3(3);
	Rrot2[0][0] = cos(throt); Rrot2[0][1] = 0.0; Rrot2[0][2] = sin(throt);
	Rrot2[1][0] = 0.0; Rrot2[1][1] = 1.0; Rrot2[1][2] = 0.0;
	Rrot2[2][0] = -sin(throt); Rrot2[2][1] = 0.0; Rrot2[2][2] = cos(throt);
	Rrot3[0][0] = cos(phirot); Rrot3[0][1] = -sin(phirot); Rrot3[0][2] = 0.0;
	Rrot3[1][0] = sin(phirot); Rrot3[1][1] = cos(phirot); Rrot3[1][2] = 0.0;
	Rrot3[2][0] = 0.0; Rrot3[2][1] = 0.0; Rrot3[2][2] = 1.0;
	tmp[0] = xvecrot[0]; tmp[1] = xvecrot[1]; tmp[2] = xvecrot[2];
	Rrot3.solve(tmp);
	Rrot2.solve(tmp); //tmp are now coordinates of xvecrot in coo frame 2
	assert(!(tmp[0] == 0.0 && tmp[1] == 0));
	if (tmp[0] == 0.0) {
		if (tmp[1] > 0.0)
			vphirot = PI / 2;
		if (tmp[1] < 0.0)
			vphirot = 3 * PI / 2;
	} else { //tmp[0] != 0.0
		vphirot = atan(tmp[1] / tmp[0]);
		if (tmp[0] > 0.0 && tmp[1] < 0.0)
			vphirot += 2 * PI;
		if (tmp[0] < 0.0)
			vphirot += PI;
	}

	//Check beta coordinates
	GenMatrix<double> Rrot(3);
	Rrot[0][0] = cos(phirot)*cos(throt)*cos(vphirot) - sin(phirot)*sin(vphirot);
	Rrot[0][1] = -cos(phirot)*cos(throt)*sin(vphirot) - sin(phirot)*cos(vphirot);
	Rrot[0][2] = cos(phirot)*sin(throt);
	Rrot[1][0] = sin(phirot)*cos(throt)*cos(vphirot) + cos(phirot)*sin(vphirot);
	Rrot[1][1] = -sin(phirot)*cos(throt)*sin(vphirot) + cos(phirot)*cos(vphirot);
	Rrot[1][2] = sin(phirot)*sin(throt);
	Rrot[2][0] = -sin(throt)*cos(vphirot);
	Rrot[2][1] = sin(throt)*sin(vphirot);
	Rrot[2][2] = cos(throt);
	tmp[0] = 0; tmp[1] = 0; tmp[2] = 1;
	Vector<double> yvecrot_(Rrot * tmp); yvecrot_ *= yrot;
	tmp[0] = sin(acos(zrot)); tmp[1] = 0; tmp[2] = zrot;
	Vector<double> xvecrot_(Rrot * tmp); xvecrot_ *= xrot;
	//cout << (xvecrot - xvecrot_).norm() << " " << (yvecrot - yvecrot_).norm() << endl;

	
	wba = acos((-sys.s[beta][alpha] * x*z + sys.c[beta][alpha] * y) / yrot);
	if ((beta-alpha+1) % 3 != 0)
		wba = 2 * PI - wba;
	//Complex rhs = zzero;
	//for (Int im = -j; im <= j; im++) {
	//	rhs += pow(-1.0, im - mbar) * \
			wignerD(j, m, im, phi, theta, vphi) * \
				wignerD(j, mbar, im, 0.0, wba, 0.0);
	//}
	//Complex lhs = wignerD(j, m, mbar, phirot, throt, vphirot);
	Complex rhs = zzero;
	for (Int im = (1-tautilde)/2; im <= j; im++) {
		rhs += 2.0*pow(-1.0, im - mbar)/ sqrt(2.0 + 2.0*DELTA(mbar, 0)) * \
			((wignerD(j, m, im, phi, theta, vphi) \
				+ tautilde * pow(-1.0, im)*wignerD(j, m, -im, phi, theta, vphi)) \
				/ sqrt(2.0 + 2.0*DELTA(im, 0))) * \
				((wignerD(j, mbar, im, 0.0, wba, 0.0) \
					+ tautilde * pow(-1.0, im)*wignerD(j, mbar, -im, 0.0, wba, 0.0)) \
					/ sqrt(2.0 + 2.0*DELTA(im, 0)));
	}
	Complex lhs = ( wignerD(j, m, mbar, phirot, throt, vphirot) \
		+ tautilde*pow(-1.0, mbar)*wignerD(j, m, -mbar, phirot, throt, vphirot)) \
			/ sqrt(2.0+2.0*DELTA(mbar, 0));
	cout << "lhs = " << lhs << ", rhs = " << rhs << endl;
	*/

//======================================
	//BlockMatr test

	/*
	vector<vector<shared_ptr<AMatrix<double>>>> blocks;
	blocks.resize(3);
	for (Int ib = 0; ib < 3; ib++)
		blocks[ib].reserve(3);

	shared_ptr<UTMatrix<double>> ut = \
		make_shared<UTMatrix<double>>(2);
	ut->set(0, 0, 4.0); ut->set(0, 1, 5.0);
	ut->set(1, 1, 6.0);
	blocks[0].push_back(ut);
	shared_ptr<GenMatrix<double>> gm = \
		make_shared<GenMatrix<double>>(2, 3);
	(*gm)[0][0] = 1.0; (*gm)[0][1] = 0.0; (*gm)[0][2] = 3.0;
	(*gm)[1][0] = 0.0; (*gm)[1][1] = 2.0; (*gm)[1][2] = 4.0;
	blocks[0].push_back(gm);
	gm = make_shared<GenMatrix<double>>(2);
	(*gm)[0][0] = 7.0; (*gm)[0][1] = 8.0;
	(*gm)[1][0] = 9.0; (*gm)[1][1] = 10.0;
	blocks[0].push_back(gm);

	gm = make_shared<GenMatrix<double>>(3, 2);
	(*gm)[0][0] = -1.0; (*gm)[0][1] = 1.0;
	(*gm)[1][0] = 2.0; (*gm)[1][1] = -1.0;
	(*gm)[2][0] = 0.5; (*gm)[2][1] = 0.0;
	blocks[1].push_back(gm);
	gm = make_shared<GenMatrix<double>>(3);
	(*gm)[0][0] = 1.0; (*gm)[0][1] = 2.0;  (*gm)[0][2] = -5.0;
	(*gm)[1][0] = 3.0; (*gm)[1][1] = 4.0; (*gm)[1][2] = -10.0;
	(*gm)[2][0] = 3.5; (*gm)[2][1] = 4.5; (*gm)[2][2] = 5.5;
	blocks[1].push_back(gm);
	gm = make_shared<GenMatrix<double>>(3, 2);
	(*gm)[0][0] = -3.0; (*gm)[0][1] = -1.0;
	(*gm)[1][0] = 2.0; (*gm)[1][1] = -11.0;
	(*gm)[2][0] = 8.0; (*gm)[2][1] = -6.0;
	blocks[1].push_back(gm);

	gm = make_shared<GenMatrix<double>>(2);
	(*gm)[0][0] = -11.0; (*gm)[0][1] = 11.0;
	(*gm)[1][0] = 22.0; (*gm)[1][1] = -11.0;
	blocks[2].push_back(gm);
	gm = make_shared<GenMatrix<double>>(2, 3);
	(*gm)[0][0] = 7.0; (*gm)[0][1] = -8.0; (*gm)[0][2] = 4.0;
	(*gm)[1][0] = 9.0; (*gm)[1][1] = -10.0; (*gm)[1][2] = 2.0;
	blocks[2].push_back(gm);
	gm = make_shared<GenMatrix<double>>(2);
	(*gm)[0][0] = 10.0; (*gm)[0][1] = 20.0;
	(*gm)[1][0] = 30.0; (*gm)[1][1] = 40.0;
	blocks[2].push_back(gm);

	vector<vector<bool>> invs; invs.resize(3);
	for (Int ib = 0; ib < 3; ib++)
		invs[ib].resize(3);
	invs[0][0] = true; invs[0][1] = false; invs[0][2] = true;
	invs[1][0] = false; invs[1][1] = true; invs[1][2] = false;
	invs[2][0] = false; invs[2][1] = false; invs[2][2] = false;
	BlockMatr<double> bm(blocks, invs);
	bm.print();

	Vector<double> vect(7);
	vect[0] = 1.0; vect[1] = 2.0; vect[2] = 3.0; vect[3] = 4.0;
	vect[4] = 5.0; vect[5] = 6.0; vect[6] = 7.0;
	vect.print();

	vect = bm * (bm * vect);
	vect.print();
	*/

//================================
	//SparseMatr test
	/*
	Int n = 5; Int m = 4; long long int nz = 10;

	SparseMatr<double> sm(n, m, nz);
	sm.aa[0] = 1.0; sm.aa[1] = 2.0; sm.aa[2] = 3.0;
	sm.aa[3] = 4.0; sm.aa[4] = 5.0;
	sm.aa[5] = 6.0; sm.aa[6] = 7.0; sm.aa[7] = 8.0;
	sm.aa[8] = 9.0; sm.aa[9] = 10.0;
	sm.ja[0] = 0; sm.ja[1] = 2; sm.ja[2] = 0;
	sm.ja[3] = 1; sm.ja[4] = 3;
	sm.ja[5] = 1; sm.ja[6] = 2; sm.ja[7] = 2;
	sm.ja[8] = 3; sm.ja[9] = 3;
	sm.ia[0] = 0; sm.ia[1] = 2; sm.ia[2] = 5; sm.ia[3] = 7;
	sm.ia[4] = 9; sm.ia[5] = 10;

	Vector<double> vec(4);
	vec[0] = 1.0; vec[1] = 2.0; vec[2] = 3.0; vec[3] = 4.0;
	Vector<double> tmp(sm * vec);
	tmp.print();

	vec[0] = 0.5; vec[1] = -0.2; vec[2] = 2.1; vec[3] = -3.3;
	tmp = sm * vec;
	tmp.print();

	sm.print();
	sm.set(2, 1, 100.0);
	cout << endl;
	sm.print();
	cout << endl;
	cout << sm.get(4, 3) << endl;
	*/
//================================
	//Order and squeeze test
	/*
	long long Int n = 10;
	Int *ja = new Int[n];
	ja[0] = 3; ja[1] = 3; ja[2] = 0; ja[3] = 2; ja[4] = 2;
	ja[5] = 5; ja[6] = 1; ja[7] = 3; ja[8] = 5; ja[9] = 2;
	//ja[0] = 3; ja[1] = 3; ja[2] = 3; ja[3] = 3; ja[4] = 3;
	//ja[5] = 3; ja[6] = 3; ja[7] = 3; ja[8] = 3; ja[9] = 3;
	double *aa = new double[n];
	aa[0] = 1.0; aa[1] = 2.0; aa[2] = 3.0; aa[3] = 4.0; aa[4] = 5.0;
	aa[5] = 6.0; aa[6] = 7.0; aa[7] = 8.0; aa[8] = 9.0; aa[9] = 10.0;
	orderAndSqueeze(n, aa, ja);
	for (long long Int i = 0; i < n; i++)
		cout << ja[i] << " " << aa[i] << endl;
	delete[] aa;
	delete[] ja;
	*/
//==========================================
	//PARDISO test
	/*
#define TEST_SP_MATR_T double //Complex //double
	SparseMatr<TEST_SP_MATR_T> sp(100, 3000);
	//SparseMatr_old<TEST_SP_MATR_T> sp(100, 3000);
	sp.fillRandom();
	Vector<TEST_SP_MATR_T> sol(100), rhs(0);
	sol.fill(1.0);
	sol = rhs = sp * sol;
	sp.solve(sol);
	cout << "Residue 1 = " << (sp * sol - rhs).norm() << endl;
	for (Int k = 0; k < sol.size(); k++)
		sol[k] = k;
	sol = rhs = sp * sol;
	sp.solve(sol);
	cout << "Residue 2 = " << (sp * sol - rhs).norm() << endl;
	for (Int k = 0; k < sol.size(); k++)
		sol[k] = 1.0/(k+1);
	sol = rhs = sp * sol;
	sp.solve(sol);
	cout << "Residue 3 = " << (sp * sol - rhs).norm() << endl;
	*/
//==========================================
	//ILU0 test
	/*
	
	Int n = 4;
	long long Int nz = 9;


	//SparseMatr_old<double> a(n, nz);
	//a.aa[0] = 1.0; a.aa[1] = 2.0; a.aa[2] = 3.0; a.aa[3] = 4.0;
	//a.aa[4] = 5.0; a.aa[5] = 6.0; a.aa[6] = 7.0; a.aa[7] = 8.0;
	//a.aa[8] = 9.0;


	SparseMatr_old<Complex> a(n, nz);
	a.aa[0] = Complex(1.0, -2.0); a.aa[1] = Complex(2.0, 3.0);
	a.aa[2] = Complex(3.0, -5.0); a.aa[3] = Complex(4.0, 0.0);
	a.aa[4] = Complex(5.0, 0.5); a.aa[5] = Complex(6.0, -1.5);
	a.aa[6] = Complex(7.0, -1.0); a.aa[7] = Complex(8.0, 10.0);
	a.aa[8] = Complex(9.0, -9.0);

	a.ja[0] = 0; a.ja[1] = 2; a.ja[2] = 0; a.ja[3] = 1;
	a.ja[4] = 3; a.ja[5] = 1; a.ja[6] = 2; a.ja[7] = 2;
	a.ja[8] = 3;
	a.ia[0] = 0; a.ia[1] = 2; a.ia[2] = 5;
	a.ia[3] = 7; a.ia[4] = 9;

	//a.print();

	auto b = ILU0(a);
	b.print();
	Vector<Complex> vvv(n);
	vvv.fill(zone); vvv[0] *= zi;
	b.solve(vvv);
	vvv.print();


	//n = 30; nz = 200;
	//SparseMatr<Complex> a_(n, nz);
	//a_.fillRandom();
	//a_.write("SPARSE_MATRIX.dat");
	//auto b_ = ILU0(a_);
	//b_.write("ILU0.dat");
	//Vector<Complex> v(n);
	//v.fill(zone); v[0] *= zi;
	//b_.solve(v);
	//v.write("LU_SOLUTION.dat");
		*/
//===================================
	//Coulomb functions test	
	
	/*
	CWFCalculator cwf(3, -2.0);
	CWFCalculator cwf2(0, 1.5);

	cout << sigmaL(3, -2.0) << endl;
	cout << sigmaL(0, 1.5) << endl;
	cout << cwf.ulp(0.5) << endl;
	cout << cwf.dulp(0.5) << endl;
	cout << cwf2.ulp(0.5) << endl;
	cout << cwf2.dulp(0.5) << endl;
	*/
//===================================
	//Hybrid basis test
	/*
	vector<double> vvv; vvv.resize(6);
	vvv[0] = 0.5; vvv[1] = 1.0; vvv[2] = 2.0;
	vvv[3] = 3.2; vvv[4] = 9.3; vvv[5] = 18.0;
	Grid gr(vvv);
	unique_ptr<ABasis<Complex>> bas = \
		make_unique<HSpline5<Complex>>(gr, Dir, none, zzero, zzero);
	//unique_ptr<ABasis<Complex>> bas = \
	//	make_unique<HSpline5<Complex>>(gr, none, none, zzero, zzero);
	//unique_ptr<ABasis<Complex>> bas = \
		make_unique<HSpline3<Complex>>(gr, Dir, none, zzero, zzero);

	double eta = 2.1;
	Int l1 = 0, l2 = 0, l3 = 1;
	vector<shared_ptr<CWFCalculator>> cwfs; cwfs.reserve(3);
	cwfs.push_back(make_shared<CWFCalculator>(l1, eta));
	cwfs.push_back(make_shared<CWFCalculator>(l2, eta));
	cwfs.push_back(make_shared<CWFCalculator>(l3, eta));
	vector<double> pn; pn.resize(3);
	pn[0] = 3.0; pn[1] = 10.0; pn[2] = 3.0;
	vector<basset> bs; bs.resize(3);
	bs[0] = basset{ 0, false }; bs[1] = basset{ 1, true };
	bs[2] = basset{1, false};
	vector<double> rnl; rnl.resize(3);
	rnl[0] = 3.0; rnl[1] = 1.4; rnl[2] = 4.0;

	//cwfs.clear(); cwfs.shrink_to_fit();
	//pn.clear(); pn.shrink_to_fit();
	//rnl.clear(); rnl.shrink_to_fit();
	unique_ptr<HybridBasis<AHermitSpline<Complex>>> hbas = \
		make_unique<HybridBasis<AHermitSpline<Complex>>>( \
			std::move(bas), cwfs, pn, rnl, bs);
	Grid cgr(hbas->collocGrid());
	gr.print();
	cgr.print();

	NumbersList limits;
	hbas->getNonzero(0.6, limits);
	for (auto num : limits)
		cout << num << "  ";
	cout << endl;

	hbas->getNonzero(1.2, limits);
	for (auto num : limits)
		cout << num << "  ";
	cout << endl;

	hbas->getNonzero(3.0, limits);
	for (auto num : limits)
		cout << num << "  ";
	cout << endl;

	hbas->getNonzero(5.0, limits);
	for (auto num : limits)
		cout << num << "  ";
	cout << endl;

	hbas->getNonzero(17.0, limits);
	for (auto num : limits)
		cout << num << "  ";
	cout << endl;

	array<shared_ptr<const ABasis<Complex>>, 1> tmp = \
		array<shared_ptr<const ABasis<Complex>>, 1>({ std::move(hbas) });
	shared_ptr<AProjDiscr<Complex, 1, Complex>> discr = \
		make_shared<AProjDiscr<Complex, 1, Complex>>(tmp);
	Function<Complex, 1> func(discr);
	func.coef.fill(zzero);
	func.coef[0] = zzero; func.coef[1] = zzero;
	func.coef[2] = zone; func.coef[3] = zzero; func.coef[4] = zzero;
	func.coef[5] = zzero; func.coef[6] = zzero; func.coef[7] = zzero;
	func.coef[8] = zzero; func.coef[9] = zzero; func.coef[10] = zzero;
	func.coef[11] = zone; func.coef[12] = zzero; func.coef[13] = zzero;
	func.coef[14] = zzero; func.coef[15] = zzero;
	//func.coef[16] = zzero;
	Int fnum = 15;
	ofstream f;
	f.open("RNL.dat");
	double a = vvv[0];
	double b = vvv[5];
	Int n = 1001;
	double h = (b - a) / (n - 1);
	double x;
	for (Int i = 0; i < n; i++) {
		x = a + i * h;
		//f << x << "  " << discr->bases[0].get()->f(x, fnum).real();
		//f << "  " << discr->bases[0].get()->d(x, fnum).real();
		//f << "  " << discr->bases[0].get()->dd(x, fnum).real() << endl;
		f << x << "  " << discr->bases[0].get()->f(x, fnum).imag();
		f << "  " << discr->bases[0].get()->d(x, fnum).imag();
		f << "  " << discr->bases[0].get()->dd(x, fnum).imag() << endl;
		//f << x << "  " << discr->bases[0].get()->f(x, func).real();
		//f << "  " << discr->bases[0].get()->d(x, func).real();
		//f << "  " << discr->bases[0].get()->dd(x, func).real() << endl;
	}
	f << endl << endl;
	f.close();
	*/
//===================================
	//Algorithms test
/*
	//vector<double> v({1.0, 2.0, 3.0, 4.0});
	vector<double> v({ 1.0, 2.0 });
	Int ind1, ind2;
	double val = 1.0;
	findIntervSorted(v.begin(), v.end(), val, ind1, ind2);
	cout << ind1 << "  " << ind2 << endl;
	Int vvv = 0;
*/
//===================================
	//Legendre and Clebsch-Gordan test
	//cout << pLegendre(1, 1, 0.0) << endl;
	//Int l1 = 5, l2 = 4, l3 = 1, m1 = -3, m2 = 2, m3 = -1;
	//cout << WignerSymbols::clebschGordan(l1, l2, l3, m1, m2, m3) << endl;

//===================================
	//Incoming wave test
	/*
	double x = 1.0, y = 20.0, theta = PI / 5.0;
	double phi, vtheta, vphi;
	phi = PI / 4; vtheta = PI / 3; vphi = 2 * PI / 5;
	double p = 2.0;
	GenMatrix<double> R(3);
	R[0][0] = cos(phi)*cos(vtheta)*cos(vphi) - sin(phi)*sin(vphi);
	R[0][1] = -cos(phi)*cos(vtheta)*sin(vphi) - sin(phi)*cos(vphi);
	R[0][2] = cos(phi)*sin(vtheta);
	R[1][0] = sin(phi)*cos(vtheta)*cos(vphi) + cos(phi)*sin(vphi);
	R[1][1] = -sin(phi)*cos(vtheta)*sin(vphi) + cos(phi)*cos(vphi);
	R[1][2] = sin(phi)*sin(vtheta);
	R[2][0] = -sin(vtheta)*cos(vphi);
	R[2][1] = sin(vtheta)*sin(vphi);
	R[2][2] = cos(vtheta);

	Vector<double> xv, yv, tmp(3);
	tmp[0] = 0.0; tmp[1] = 0.0; tmp[2] = 1.0;
	xv = R * tmp; xv *= x;
	tmp[0] = sin(theta); tmp[1] = 0.0; tmp[2] = cos(theta);
	yv = R * tmp; yv *= y;

	cout << "3D value: " << y * exp(zi*p*yv[2]) / sqrt(4 * PI) << endl;
	
	Int jmax = 100;
	Complex sum = zzero;
	Complex fmmjt, dmmj1, dmmj2;
	for (Int j = 0; j <= jmax; j++) {
		CWFCalculator cwf(j, 0.0);
		for (Int im = 0; im <= j; im++) {
			dmmj1 = wignerDSmall(j, 0, im, vtheta)*exp(-zi*im*vphi);
			dmmj2 = wignerDSmall(j, 0, -im, vtheta)*exp(zi * im*vphi);
			fmmjt = (dmmj1 + pow(-1.0, im)*dmmj2) / sqrt(2.0+2.0*DELTA(im, 0));
			sum += zone / sqrt(2.0 + 2.0*DELTA(im, 0)) * sqrt(2.0*j + 1) * \
				pow(zi, j) * cwf.fl(p*y) / p * pLegendre(j, im, cos(theta)) * \
				WignerSymbols::clebschGordan(j, 0, j, 0, 0, 0) * \
				WignerSymbols::clebschGordan(j, 0, j, im, 0, im) * 2.0 * fmmjt;
		}
		cout << "Sum after j = " << j << ": " << sum << endl;
	}
	cout << "Partial value: " << sum << endl;

	Int lambda = 5, mlam = 1;
	double phiy = atan2(yv[1], yv[0]);
	if (phiy < 0)
		phiy += 2 * PI;
	cout << pLegendre(lambda, mlam, yv[2] / y)*exp(zi*mlam*phiy) << endl;
	Complex term;
	sum = zzero;
	for (Int imbar = -lambda; imbar <= lambda; imbar++) {
		dmmj1 = exp(-zi * mlam*phi)*wignerDSmall(lambda, mlam, imbar, vtheta)*exp(-zi * imbar*vphi);
		//term = pow(-1.0, mlam)*conj(dmmj1)*pow(-1.0, imbar)*pLegendre(lambda, abs(imbar), cos(theta));
		term = pow(-1.0, mlam)*conj(dmmj1)*pLegendre(lambda, abs(imbar), cos(theta));
		if (imbar < 0)
			term *= pow(-1.0, imbar);
		sum += term;
	}
	cout << sum << endl;
	*/
//===================================
	//Legendre polynomials and basis
	/*
	double a = -1.0, b = 1.0;
	//a = -5.0; b = 3.0;
	Grid gr(a, b, 2);
	Int deg = 5;
	
	array<shared_ptr<const ABasis<double>>, 1> barr;
	barr[0] = make_shared<const SimplePoly<double>>(gr, deg);
	shared_ptr<const ADiscretization<double, 1>> discr = \
		make_shared<const Collocation<double, 1, double>>(barr);

	Int l = 4;
	double x = 0.7;
	cout << "P_l = " << barr[0]->f(x, l) << ", dP_l = " << barr[0]->d(x, l);
	cout << ", ddP_l = " << barr[0]->dd(x, l) << endl;
	cout << "P_l = " << pLegendre(l, x) << ", dP_l = " << dpLegendre(l, x);
	cout << ", ddP_l = " << ddpLegendre(l, x) << endl;

	Function<double, 1> func(discr);
	func.coef.fill(0.0);
	func.coef[0] = 1.0;
	func.coef[2] = -2.0;
	func.coef[3] = 3.0;
	func.coef[deg] = -1.0;

	cout << "val = " << barr[0]->f(x, func) << ", dval = " << barr[0]->d(x, func);
	cout << ", ddval = " << barr[0]->dd(x, func) << endl;

	GenMatrix<double> ovlp(barr[0]->getOverlap());
	ovlp.print(); cout << endl;
	GenMatrix<double> ovlpw(barr[0]->getOverlapW([](double x) {return x * x; }, 2));
	ovlpw.print();

	barr[0]->collocGrid().print();
	*/
}