/*
 * formMatrix_test.cpp
 *
 *  Created on: Jun 18, 2014
 *      Author: pxiang
 */
#include "gtest/gtest.h"
#include "formMatrix.h"

TEST(SetUpInteraction, CheckMagnitude) {
	extern Interaction *pInteraction;
	LatticeShape lattice1D(1);
	int xmax = 200;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	double onSiteE = 1.0;
	double hop = 2.0;
	double dyn = 3.0;
	bool randomOnsite = true;
	bool randomHop = false;
	bool randomDyn = false;
	int maxDistance = 10;
	unsigned seed = 22667;
	bool longRangeHop = true;
	bool longRangeDyn = true;

	InteractionData interactionData = {onSiteE,hop,dyn,randomOnsite,
			                           randomHop,randomDyn,maxDistance,seed,
			                           longRangeHop,longRangeDyn};
//	generateIndexMatrix(lattice1D);

	setInteractions(lattice1D, interactionData);

	// check whether the values of hop and dyn are in the right range
	for (int i=0; i<=xmax; ++i) {
		for (int j=0; j<=xmax; ++j) {
			double t_element = tMatrix(*pInteraction, i,j);
			double d_element = dMatrix(*pInteraction, i,j);
			if (i==j) {
				EXPECT_DOUBLE_EQ(t_element, 0.0);
				EXPECT_DOUBLE_EQ(d_element, 0.0);
			} else {
				// check the range of the matrix elements
				EXPECT_TRUE(t_element>=0.0 && t_element<=hop);
				EXPECT_TRUE(d_element>=0.0 && d_element<=dyn);
				// check the symmetry
				EXPECT_DOUBLE_EQ(t_element, tMatrix(*pInteraction, j,i));
				EXPECT_DOUBLE_EQ(d_element, dMatrix(*pInteraction, j,i));

				// check the long range behavior
				if (longRangeHop) {

					if (std::abs(i-j)<=maxDistance) {
						double max = hop/std::pow(std::abs(i-j), 3.0);
						EXPECT_TRUE(t_element>=0.0 && t_element<=max);
					} else {
						EXPECT_DOUBLE_EQ(t_element, 0.0);
					}
				} else {
					if (std::abs(i-j)==1) {
						EXPECT_TRUE(t_element>=0.0 && t_element<=hop);
					} else {
						EXPECT_DOUBLE_EQ(t_element, 0.0);
					}
				}

				if (longRangeDyn) {

					if (std::abs(i-j)<=maxDistance) {
						double max = dyn/std::pow(std::abs(i-j), 3.0);
						EXPECT_TRUE(t_element>=0.0 && t_element<=max);
					} else {
						EXPECT_DOUBLE_EQ(d_element, 0.0);
					}
				} else {
					if (std::abs(i-j)==1) {
						EXPECT_TRUE(t_element>=0.0 && t_element<=dyn);
					} else {
						EXPECT_DOUBLE_EQ(d_element, 0.0);
					}
				}


			}
		}
	}

	// check whether the values of onSiteE are in the right range
	for (int i=0; i<=xmax; ++i) {
		double e_element = eVector(*pInteraction, i);
		if (randomOnsite) {
			EXPECT_TRUE(e_element>=-onSiteE/2 && e_element<=onSiteE/2);
		} else {
			EXPECT_DOUBLE_EQ(e_element, onSiteE);
		}
	}

}


TEST(SetUpLongRange, Hop) {
	extern Interaction *pInteraction;
	LatticeShape lattice1D(1);
	int xmax = 200;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,5,230,true,true};
	generateIndexMatrix(lattice1D);

	setInteractions(lattice1D, interactionData);
	EXPECT_TRUE(pInteraction->getMaxDistanceHop()==5);
	EXPECT_TRUE(pInteraction->getMaxDistanceDyn()==5);

	interactionData.longRangeHop = true;
	interactionData.longRangeDyn = false;
	setInteractions(lattice1D, interactionData);
	EXPECT_TRUE(pInteraction->getMaxDistanceHop()==5);
	EXPECT_TRUE(pInteraction->getMaxDistanceDyn()==1);

	interactionData.longRangeHop = false;
	interactionData.longRangeDyn = true;
	setInteractions(lattice1D, interactionData);
	EXPECT_TRUE(pInteraction->getMaxDistanceHop()==1);
	EXPECT_TRUE(pInteraction->getMaxDistanceDyn()==5);

	interactionData.longRangeHop = false;
	interactionData.longRangeDyn = false;
	setInteractions(lattice1D, interactionData);
	EXPECT_TRUE(pInteraction->getMaxDistanceHop()==1);
	EXPECT_TRUE(pInteraction->getMaxDistanceDyn()==1);
}


TEST(FormMatrixZ, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 200;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,5,230,true,true};
	generateIndexMatrix(lattice1D);
	setInteractions(lattice1D, interactionData);
	int K = 24;
	dcomplex energy = dcomplex(10.0, 0.001);
	CDMatrix ZK;
	formMatrixZ(K, energy, ZK);
	EXPECT_TRUE(true);
}


TEST(FormMatrixM, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 200;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,5,230,true,true};
	generateIndexMatrix(lattice1D);
	setInteractions(lattice1D, interactionData);
	CDMatrix MKKp;
	for (int K=1; K<=xmax+xmax-1; ++K) {
		for (int Kp=1; Kp<=xmax+xmax-1; ++Kp) {
//			std::cout << "K = " << K << "\t " << "Kp = " << Kp << std::endl;
			formMatrixM(K, Kp, MKKp);
		}
	}

	EXPECT_TRUE(true);
}


TEST(FormMatrixW, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 200;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	generateIndexMatrix(lattice1D);


	dcomplex energy = dcomplex(10.0, 0.001);

	for (int K=1; K<=xmax+xmax-1; ++K) {
		for (int maxDistance=1; maxDistance<=5; ++maxDistance) {
			InteractionData interactionData = {1.0,1.0,1.0,true,false,true,
					                           maxDistance,230,true,true};
			setInteractions(lattice1D, interactionData);
			CDMatrix WK;
			formMatrixW(K, energy, WK);
		}
	}


	EXPECT_TRUE(true);
}


TEST(FormMatrixAlpha, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 200;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	generateIndexMatrix(lattice1D);



	for (int maxDistance=1; maxDistance<=5; ++maxDistance) {
		//std::cout <<"numNeighbor = " << numNeighbor << std::endl;
		InteractionData interactionData = {1.0,1.0,1.0,true,false,true,
				                           maxDistance,230,true,true};
		setInteractions(lattice1D, interactionData);
		CDMatrix Alpha;
		// alpha_K: K starts from 2, K=1 is meaningless because V_0 doesn't exist
		for (int K=2; K<=xmax+xmax-1; ++K) {
			if (K-maxDistance>=1) {
				//std::cout << "K="<< K << " Kp=" << K-numNeighbor << std::endl;
				formMatrixAlpha(K, Alpha);
			}
		}
	}

	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,1,230,true,true};
	int rows, cols;
	int rows_expected, cols_expected;
	int K;
	int Kmax = xmax + xmax - 1;
	CDMatrix Alpha;

	K = xmax+xmax-1; // this is Kmax
	interactionData.maxDistance = 1;
	setInteractions(lattice1D, interactionData);
	formMatrixAlpha(K, Alpha);
	rows = Alpha.rows();
	cols = Alpha.cols();
	rows_expected = DimsOfV[K];
	cols_expected = DimsOfV[K-1];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	interactionData.maxDistance = 2;
	setInteractions(lattice1D, interactionData);
	formMatrixAlpha(K-1, Alpha);
	rows = Alpha.rows();
	cols = Alpha.cols();
	rows_expected = DimsOfV[K-1] + DimsOfV[K];
	cols_expected = DimsOfV[K-3] + DimsOfV[K-2];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	interactionData.maxDistance = 3;
	setInteractions(lattice1D, interactionData);
	formMatrixAlpha(K-2, Alpha);
	rows = Alpha.rows();
	cols = Alpha.cols();
	rows_expected = DimsOfV[K-2] + DimsOfV[K-1] + DimsOfV[K];
	cols_expected = DimsOfV[K-5] + DimsOfV[K-4] + DimsOfV[K-3];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	interactionData.maxDistance = 4;
	setInteractions(lattice1D, interactionData);
	formMatrixAlpha(K-3, Alpha);
	rows = Alpha.rows();
	cols = Alpha.cols();
	rows_expected = DimsOfV[K-3] + DimsOfV[K-2] + DimsOfV[K-1] + DimsOfV[K];
	cols_expected = DimsOfV[K-7] + DimsOfV[K-6] + DimsOfV[K-5] + DimsOfV[K-4];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);


	interactionData.maxDistance = 4;
	setInteractions(lattice1D, interactionData);
	K = Kmax-1;
	formMatrixAlpha(K, Alpha);
	rows = Alpha.rows();
	cols = Alpha.cols();
	rows_expected = DimsOfV[K] + DimsOfV[K+1];
	cols_expected = DimsOfV[K-4] + DimsOfV[K-3] + DimsOfV[K-2] + DimsOfV[K-1];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

}



TEST(FormMatrixBeta, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 200;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	generateIndexMatrix(lattice1D);
	int Kmax = xmax+xmax-1;

	for (int maxDistance=1; maxDistance<=5; ++maxDistance) {
		InteractionData interactionData = {1.0,1.0,1.0,true,false,true,
				                           maxDistance,230,true,true};
		setInteractions(lattice1D, interactionData);
		CDMatrix Beta;
		// Beta_K: K is in range [1, Kmax-1] Beta_{Kmax} is meaningless
		for (int K=1; K<=Kmax-1; ++K) {
			if (K+maxDistance<=Kmax) {
				//std::cout << "K="<<K<< " Kp=" << K+numNeighbor << std::endl;
				formMatrixBeta(K, Beta);
			}
		}
	}

//	std::cout << "\n\n" << std::endl;
	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,1,230,true,true};
	int rows, cols;
	int rows_expected, cols_expected;
	int K = 1;
	CDMatrix Beta;
	interactionData.maxDistance = 1;
	setInteractions(lattice1D, interactionData);
	formMatrixBeta(K, Beta);
	rows = Beta.rows();
	cols = Beta.cols();
	rows_expected = DimsOfV[K];
	cols_expected = DimsOfV[K+1];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	K = 1;
	interactionData.maxDistance = 2;
	setInteractions(lattice1D, interactionData);
	formMatrixBeta(K, Beta);
	rows = Beta.rows();
	cols = Beta.cols();
	rows_expected = DimsOfV[K] + DimsOfV[K+1];
	cols_expected = DimsOfV[K+2] + DimsOfV[K+1+2];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	interactionData.maxDistance = 3;
	setInteractions(lattice1D, interactionData);
	formMatrixBeta(K, Beta);
	rows = Beta.rows();
	cols = Beta.cols();
	rows_expected = DimsOfV[K] + DimsOfV[K+1] + DimsOfV[K+2];
	cols_expected = DimsOfV[K+3] + DimsOfV[K+1+3] + DimsOfV[K+2+3];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	interactionData.maxDistance = 4;
	setInteractions(lattice1D, interactionData);
	formMatrixBeta(K, Beta);
	rows = Beta.rows();
	cols = Beta.cols();
	rows_expected = DimsOfV[K] + DimsOfV[K+1] + DimsOfV[K+2] + DimsOfV[K+3];
	cols_expected = DimsOfV[K+4] + DimsOfV[K+1+4] + DimsOfV[K+2+4]+ DimsOfV[K+3+4];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	interactionData.maxDistance = 5;
	setInteractions(lattice1D, interactionData);
	formMatrixBeta(K,  Beta);
	rows = Beta.rows();
	cols = Beta.cols();
	rows_expected = DimsOfV[K] + DimsOfV[K+1] + DimsOfV[K+2] + DimsOfV[K+3]
	                + DimsOfV[K+4];
	cols_expected = DimsOfV[K+5] + DimsOfV[K+1+5] + DimsOfV[K+2+5]+ DimsOfV[K+3+5]
	                + DimsOfV[K+4 + 5];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	interactionData.maxDistance = 6;
	setInteractions(lattice1D, interactionData);
	formMatrixBeta(K,  Beta);
	rows = Beta.rows();
	cols = Beta.cols();
	rows_expected = DimsOfV[K] + DimsOfV[K+1] + DimsOfV[K+2] + DimsOfV[K+3]
	                + DimsOfV[K+4] + DimsOfV[K+5];
	cols_expected = DimsOfV[K+6] + DimsOfV[K+1+6] + DimsOfV[K+2+6]+ DimsOfV[K+3+6]
	                + DimsOfV[K+4+6] + DimsOfV[K+5+6];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

}

