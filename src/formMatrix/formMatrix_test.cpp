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
	Parameters parameters = {1.0,1.0,1.0,true,false,true,5,230};
	generateIndexMatrix(lattice1D);
	setInteractions(lattice1D, parameters);
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
	Parameters parameters = {1.0,1.0,1.0,true,false,true,5,230};
	generateIndexMatrix(lattice1D);
	setInteractions(lattice1D, parameters);
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
	Parameters parameters = {1.0,1.0,1.0,true,false,true,5,230};
	generateIndexMatrix(lattice1D);
	setInteractions(lattice1D, parameters);
	CDMatrix WK;
	dcomplex energy = dcomplex(10.0, 0.001);

	for (int K=1; K<=xmax+xmax-1; ++K) {
		for (int numNeighbor=1; numNeighbor<=5 && K+numNeighbor<=xmax+xmax-1; ++numNeighbor) {
			formMatrixW(K, numNeighbor, energy, WK);
		}
	}


	int K = 1;
	int numNeighbor = 2;
	formMatrixW(K, numNeighbor, energy, WK);

	numNeighbor = 3;
	formMatrixW(K, numNeighbor, energy, WK);

	numNeighbor = 4;
	formMatrixW(K, numNeighbor, energy, WK);

	numNeighbor = 5;
	formMatrixW(K, numNeighbor, energy, WK);

	numNeighbor = 6;
	formMatrixW(K, numNeighbor, energy, WK);

	EXPECT_TRUE(true);
}


TEST(FormMatrixAlpha, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 200;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Parameters parameters = {1.0,1.0,1.0,true,false,true,5,230};
	generateIndexMatrix(lattice1D);
	setInteractions(lattice1D, parameters);
	CDMatrix Alpha;

	for (int numNeighbor=1; numNeighbor<=5; ++numNeighbor) {
		//std::cout <<"numNeighbor = " << numNeighbor << std::endl;
		for (int K=1; K<=xmax+xmax-1; ++K) {
			if (K + numNeighbor-1<=xmax+xmax-1 && K-numNeighbor>=1) {
				//std::cout << "K="<< K << " Kp=" << K-numNeighbor << std::endl;
				formMatrixAlpha(K, numNeighbor, Alpha);
			}
		}
	}

	EXPECT_TRUE(true);
}



TEST(FormMatrixBeta, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 200;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Parameters parameters = {1.0,1.0,1.0,true,false,true,5,230};
	generateIndexMatrix(lattice1D);
	setInteractions(lattice1D, parameters);
	CDMatrix Beta;

	for (int numNeighbor=1; numNeighbor<=5; ++numNeighbor) {
		for (int K=1; K<=xmax+xmax-1; ++K) {
			if (K + 2*numNeighbor-1<=xmax+xmax-1) {
				//std::cout << "K="<<K<< " Kp=" << K+numNeighbor << std::endl;
				formMatrixBeta(K, numNeighbor, Beta);
			}
		}
	}

//	std::cout << "\n\n" << std::endl;
	int rows, cols;
	int rows_expected, cols_expected;
	int K = 1;
	int numNeighbor = 1;
	formMatrixBeta(K, numNeighbor, Beta);
	rows = Beta.rows();
	cols = Beta.cols();
	rows_expected = DimsOfV[K];
	cols_expected = DimsOfV[K+1];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	K = 1;
	numNeighbor = 2;
	formMatrixBeta(K, numNeighbor, Beta);
	rows = Beta.rows();
	cols = Beta.cols();
	rows_expected = DimsOfV[K] + DimsOfV[K+1];
	cols_expected = DimsOfV[K+2] + DimsOfV[K+1+2];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	numNeighbor = 3;
	formMatrixBeta(K, numNeighbor, Beta);
	rows = Beta.rows();
	cols = Beta.cols();
	rows_expected = DimsOfV[K] + DimsOfV[K+1] + DimsOfV[K+2];
	cols_expected = DimsOfV[K+3] + DimsOfV[K+1+3] + DimsOfV[K+2+3];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	numNeighbor = 4;
	formMatrixBeta(K, numNeighbor, Beta);
	rows = Beta.rows();
	cols = Beta.cols();
	rows_expected = DimsOfV[K] + DimsOfV[K+1] + DimsOfV[K+2] + DimsOfV[K+3];
	cols_expected = DimsOfV[K+4] + DimsOfV[K+1+4] + DimsOfV[K+2+4]+ DimsOfV[K+3+4];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	numNeighbor = 5;
	formMatrixBeta(K, numNeighbor, Beta);
	rows = Beta.rows();
	cols = Beta.cols();
	rows_expected = DimsOfV[K] + DimsOfV[K+1] + DimsOfV[K+2] + DimsOfV[K+3]
	                + DimsOfV[K+4];
	cols_expected = DimsOfV[K+5] + DimsOfV[K+1+5] + DimsOfV[K+2+5]+ DimsOfV[K+3+5]
	                + DimsOfV[K+4 + 5];
	EXPECT_EQ(rows, rows_expected);
	EXPECT_EQ(cols, cols_expected);

	numNeighbor = 6;
	formMatrixBeta(K, numNeighbor, Beta);
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

