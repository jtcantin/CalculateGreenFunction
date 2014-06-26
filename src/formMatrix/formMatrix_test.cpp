/*
 * formMatrix_test.cpp
 *
 *  Created on: Jun 18, 2014
 *      Author: pxiang
 */
#include "gtest/gtest.h"
#include "formMatrix.h"


TEST(FormMatrixZ, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 200;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,5,230};
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
	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,5,230};
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
		for (int maxDistance=1; maxDistance<=5 && K+maxDistance<=xmax+xmax-1; ++maxDistance) {
			InteractionData interactionData = {1.0,1.0,1.0,true,false,true,maxDistance,230};
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
		InteractionData interactionData = {1.0,1.0,1.0,true,false,true,maxDistance,230};
		setInteractions(lattice1D, interactionData);
		CDMatrix Alpha;
		for (int K=1; K<=xmax+xmax-1; ++K) {
			if (K + maxDistance-1<=xmax+xmax-1 && K-maxDistance>=1) {
				//std::cout << "K="<< K << " Kp=" << K-numNeighbor << std::endl;
				formMatrixAlpha(K, Alpha);
			}
		}
	}

	EXPECT_TRUE(true);
}



TEST(FormMatrixBeta, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 200;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	generateIndexMatrix(lattice1D);


	for (int maxDistance=1; maxDistance<=5; ++maxDistance) {
		InteractionData interactionData = {1.0,1.0,1.0,true,false,true,maxDistance,230};
		setInteractions(lattice1D, interactionData);
		CDMatrix Beta;
		for (int K=1; K<=xmax+xmax-1; ++K) {
			if (K + 2*maxDistance-1<=xmax+xmax-1) {
				//std::cout << "K="<<K<< " Kp=" << K+numNeighbor << std::endl;
				formMatrixBeta(K, Beta);
			}
		}
	}

//	std::cout << "\n\n" << std::endl;
	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,1,230};
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

	EXPECT_TRUE(true);
}

