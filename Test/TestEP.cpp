#include "TestEP.h"
#include "HSpline5.h"
#include "Collocation.h"

TestEP::TestEP() {
#ifdef TEST_EP_SPARSE_BIG
	//Create discretization

	Int npoi = 200;
	Grid gr(-1.0, 1.0, npoi);

	array<shared_ptr<const ABasis<double>>, 1> bases;
	bases[0] = make_shared<const HSpline5<double>>(gr, Dir, Dir, 1.0, 1.0);

	shared_ptr<const ADiscretization<TEST_OP_SP_TYPE, 1>> d = \
		make_shared<const Collocation<TEST_OP_SP_TYPE, 1, double>>(bases);

	//============================================================

	op = new TestOpSparse(d);
	TestOpSparse &op_ = dynamic_cast<TestOpSparse &>(*op);

	op_.matr.fillRandom();
	//op_.matr.write("SPARSE_MATR.dat");

	evec.reserve(1);
	evec.push_back( \
		make_unique<Function<TEST_OP_SP_TYPE, 1>>(d));
	evec[0]->coef.setZero();

#else

	//Create discretization

	vector<double> poi; poi.resize(3);//poi.resize(6);
	poi[0] = -7.8; poi[1] = -4.6; poi[2] = -1.0;
	//poi[3] = 2.2; poi[4] = 9.3; poi[5] = 18.0;
	//Grid gr(6, poi);
	Grid gr(poi);

	array<shared_ptr<const ABasis>, 1> bases;
	bases[0] = make_shared<const HSpline5<double>>(gr, Dir, Dir, 1.0, 1.0);

	shared_ptr<const ADiscretization<TEST_OP_TYPE, 1>> d = \
		make_shared<const Collocation<TEST_OP_TYPE, 1>>(bases);

	//============================================================

	op = make_shared<TestOperator>(d);
	shared_ptr<TestOperator> op_ = dynamic_pointer_cast<TestOperator, AnOperator<TEST_OP_TYPE>>(op);

	evec.reserve(1);
	evec.push_back(make_shared<Function<TEST_OP_TYPE, 1>>(d));

	//double gen problem

	op_->matr[0][0] = 10; op_->matr[0][1] = 0.546881519204984; op_->matr[0][2] = 0.800280468888800; op_->matr[0][3] = 0.0357116785741896; op_->matr[0][4] = 0.655477890177557; op_->matr[0][5] = 0.823457828327293; op_->matr[0][6] = 0.765516788149002;
	op_->matr[1][0] = 0.905791937075619; op_->matr[1][1] = 11; op_->matr[1][2] = 0.141886338627215; op_->matr[1][3] = 0.849129305868777; op_->matr[1][4] = 0.171186687811562; op_->matr[1][5] = 0.694828622975817; op_->matr[1][6] = 0.795199901137063;
	op_->matr[2][0] = 0.126986816293506; op_->matr[2][1] = 0.964888535199277; op_->matr[2][2] = 12; op_->matr[2][3] = 0.933993247757551; op_->matr[2][4] = 0.706046088019609; op_->matr[2][5] = 0.317099480060861; op_->matr[2][6] = 0.186872604554379;
	op_->matr[3][0] = 0.913375856139019; op_->matr[3][1] = 0.157613081677548; op_->matr[3][2] = 0.915735525189067; op_->matr[3][3] = 13; op_->matr[3][4] = 0.0318328463774207; op_->matr[3][5] = 0.950222048838355; op_->matr[3][6] = 0.489764395788231;
	op_->matr[4][0] = 0.632359246225410; op_->matr[4][1] = 0.970592781760616; op_->matr[4][2] = 0.792207329559554; op_->matr[4][3] = 0.757740130578333; op_->matr[4][4] = 14; op_->matr[4][5] = 0.0344460805029088; op_->matr[4][6] = 0.445586200710900;
	op_->matr[5][0] = 0.0975404049994095; op_->matr[5][1] = 0.957166948242946; op_->matr[5][2] = 0.959492426392903; op_->matr[5][3] = 0.743132468124916; op_->matr[5][4] = 0.0461713906311539; op_->matr[5][5] = 15; op_->matr[5][6] = 0.646313010111265;
	op_->matr[6][0] = 0.278498218867048; op_->matr[6][1] = 0.485375648722841; op_->matr[6][2] = 0.655740699156587; op_->matr[6][3] = 0.392227019534168; op_->matr[6][4] = 0.0971317812358475; op_->matr[6][5] = 0.381558457093008; op_->matr[6][6] = 16;


	/*
	A = [10 0.546881519204984 0.800280468888800 0.0357116785741896 0.655477890177557 0.823457828327293 0.765516788149002;
		0.905791937075619 11 0.141886338627215 0.849129305868777 0.171186687811562 0.694828622975817 0.795199901137063;
		0.126986816293506 0.964888535199277 12 0.933993247757551 0.706046088019609 0.317099480060861 0.186872604554379;
		0.913375856139019 0.157613081677548 0.915735525189067 13 0.0318328463774207 0.950222048838355 0.489764395788231;
		0.632359246225410 0.970592781760616 0.792207329559554 0.757740130578333 14 0.0344460805029088 0.445586200710900;
		0.0975404049994095 0.957166948242946 0.959492426392903 0.743132468124916 0.0461713906311539 15 0.646313010111265;
		0.278498218867048 0.485375648722841 0.655740699156587 0.392227019534168 0.0971317812358475 0.381558457093008 16];
	f = [0.754686681982361;0.276025076998578;0.679702676853675;0.655098003973841;0.162611735194631;0.118997681558377;0.498364051982143];
	*/


	//complex gen problem

		/*
		op_->matr[0][0] = Complex(10, 11); op_->matr[0][1] = Complex(0.54, 0.95); op_->matr[0][2] = Complex(0.80, 0.69); op_->matr[0][3] = Complex(0.03, 0.84); op_->matr[0][4] = Complex(0.65, 0.25); op_->matr[0][5] = Complex(0.82, 0.91); op_->matr[0][6] = Complex(0.76, 0.05);
		op_->matr[1][0] = Complex(0.90,  0.27); op_->matr[1][1] = Complex(12, 13); op_->matr[1][2] = Complex(0.14, 0.89); op_->matr[1][3] = Complex(0.84, 0.25); op_->matr[1][4] = Complex(0.17, 0.61); op_->matr[1][5] = Complex(0.69, 0.28); op_->matr[1][6] = Complex(0.79, 0.53);
		op_->matr[2][0] = Complex(0.12, 0.67); op_->matr[2][1] = Complex(0.96, 0.58); op_->matr[2][2] = Complex(13, 12); op_->matr[2][3] = Complex(0.93, 0.81); op_->matr[2][4] = Complex(0.70, 0.47); op_->matr[2][5] = Complex(0.31, 0.75); op_->matr[2][6] = Complex(0.18, 0.77);
		op_->matr[3][0] = Complex(0.91, 0.65); op_->matr[3][1] = Complex(0.15, 0.22); op_->matr[3][2] = Complex(0.91, 0.54); op_->matr[3][3] = Complex(10, 9); op_->matr[3][4] = Complex(0.03, 0.35); op_->matr[3][5] = Complex(0.95, 0.75); op_->matr[3][6] = Complex(0.48, 0.93);
		op_->matr[4][0] = Complex(0.63, 0.16); op_->matr[4][1] = Complex(0.97, 0.75); op_->matr[4][2] = Complex(0.79, 0.13); op_->matr[4][3] = Complex(0.75, 0.92); op_->matr[4][4] = Complex(14, 15); op_->matr[4][5] = Complex(0.03, 0.38); op_->matr[4][6] = Complex(0.44, 0.12);
		op_->matr[5][0] = Complex(0.09, 0.11); op_->matr[5][1] = Complex(0.95, 0.25); op_->matr[5][2] = Complex(0.95, 0.14); op_->matr[5][3] = Complex(0.74, 0.34); op_->matr[5][4] = Complex(0.04, 0.58); op_->matr[5][5] = Complex(11, 11); op_->matr[5][6] = Complex(0.64, 0.56);
		op_->matr[6][0] = Complex(0.27, 0.49); op_->matr[6][1] = Complex(0.48, 0.50); op_->matr[6][2] = Complex(0.65, 0.25); op_->matr[6][3] = Complex(0.39, 0.19); op_->matr[6][4] = Complex(0.09, 0.54); op_->matr[6][5] = Complex(0.38, 0.07); op_->matr[6][6] = Complex(15, 16);
		
		*/
		/*
		A = [10 + 11i 0.54 + 0.95i 0.80 + 0.69i 0.03 + 0.84i 0.65 + 0.25i 0.82 + 0.91i 0.76 + 0.05i;
			0.90 + 0.27i 12 + 13i 0.14 + 0.89i 0.84 + 0.25i 0.17 + 0.61i 0.69 + 0.28i 0.79 + 0.53i;
			0.12 + 0.67i 0.96 + 0.58i 13 + 12i 0.93 + 0.81i 0.70 + 0.47i 0.31 + 0.75i 0.18 + 0.77i;
			0.91 + 0.65i 0.15 + 0.22i 0.91 + 0.54i 10 + 9i 0.03 + 0.35i 0.95 + 0.75i 0.48 + 0.93i;
			0.63 + 0.16i 0.97 + 0.75i 0.79 + 0.13i 0.75 + 0.92i 14 + 15i 0.03 + 0.38i 0.44 + 0.12i;
			0.09 + 0.11i 0.95 + 0.25i 0.95 + 0.14i 0.74 + 0.34i 0.04 + 0.58i 11 + 11i 0.64 + 0.56i;
			0.27 + 0.49i 0.48 + 0.50i 0.65 + 0.25i 0.39 + 0.19i 0.09 + 0.54i 0.38 + 0.07i 15 + 16i];
		f = [0.01 + 0.60i;0.33 + 0.26i;0.16 + 0.65i;0.79 + 0.68i;0.31 + 0.74i;0.52 + 0.45i;0.16 + 0.08i];
		*/
#endif
}


TestEP::~TestEP() { }
