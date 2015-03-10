/*
 * binaryIO_test.cpp
 *
 *  Created on: Jan 2, 2014
 *      Author: pxiang
 */
#include "gtest/gtest.h"
#include "binaryIO.h"

TEST(BinaryIO, DISABLED_TwoByTwoMatirx) {
	CDMatrix cm(2,2);
	cm(0,0) = dcomplex(1.0, 0.5);

	cm(0,1) = dcomplex(1.3, 0.3);

	cm(1,0) = dcomplex(1.89, 0.35);

	cm(1,1) = dcomplex(1.07899, 0.1135);

	saveMatrixBin("cm.bin",cm);
	CDMatrix cm2;
	loadMatrixBin("cm.bin",cm2);
	EXPECT_EQ(cm.rows(), cm2.rows());
	EXPECT_EQ(cm.cols(), cm2.cols());
	for (int i=0; i<cm.rows(); ++i) {
		for (int j=0; j<cm.cols(); ++j) {
			dcomplex c1 = cm(i,j);
			dcomplex c2 = cm2(i,j);
			EXPECT_DOUBLE_EQ(c1.real(), c2.real());
			EXPECT_DOUBLE_EQ(c1.imag(), c2.imag());
		}
	}
}


TEST(BinaryIO, DISABLED_TwoByOneMatirx) {
	CDMatrix cm(2,1);
	cm(0,0) = dcomplex(1.0, 0.5);
	cm(1,0) = dcomplex(1.89, 0.35);
	saveMatrixBin("cm.bin",cm);
	CDMatrix cm2;
	loadMatrixBin("cm.bin",cm2);
	EXPECT_EQ(cm.rows(), cm2.rows());
	EXPECT_EQ(cm.cols(), cm2.cols());
	for (int i=0; i<cm.rows(); ++i) {
		for (int j=0; j<cm.cols(); ++j) {
			dcomplex c1 = cm(i,j);
			dcomplex c2 = cm2(i,j);
			EXPECT_DOUBLE_EQ(c1.real(), c2.real());
			EXPECT_DOUBLE_EQ(c1.imag(), c2.imag());
		}
	}
}


TEST(BinaryIO, DISABLED_OneByTwoMatirx) {
	CDMatrix cm(1,2);
	cm(0,0) = dcomplex(1.0, 0.5);
	cm(0,1) = dcomplex(1.3, 0.3);
	saveMatrixBin("cm.bin",cm);
	CDMatrix cm2;
	loadMatrixBin("cm.bin",cm2);
	EXPECT_EQ(cm.rows(), cm2.rows());
	EXPECT_EQ(cm.cols(), cm2.cols());
	for (int i=0; i<cm.rows(); ++i) {
		for (int j=0; j<cm.cols(); ++j) {
			dcomplex c1 = cm(i,j);
			dcomplex c2 = cm2(i,j);
			EXPECT_DOUBLE_EQ(c1.real(), c2.real());
			EXPECT_DOUBLE_EQ(c1.imag(), c2.imag());
		}
	}
}


TEST(BinaryIO, DISABLED_OneByOneMatirx) {
	CDMatrix cm(1,1);
	cm(0,0) = dcomplex(1.0, 0.5);
	saveMatrixBin("cm.bin",cm);
	CDMatrix cm2;
	loadMatrixBin("cm.bin",cm2);
	EXPECT_EQ(cm.rows(), cm2.rows());
	EXPECT_EQ(cm.cols(), cm2.cols());
	for (int i=0; i<cm.rows(); ++i) {
		for (int j=0; j<cm.cols(); ++j) {
			dcomplex c1 = cm(i,j);
			dcomplex c2 = cm2(i,j);
			EXPECT_DOUBLE_EQ(c1.real(), c2.real());
			EXPECT_DOUBLE_EQ(c1.imag(), c2.imag());
		}
	}
}


TEST(BinaryIO, DISABLED_ThousandByThousandMatirx) {
	CDMatrix cm(1000,1000);
	RandomNumberGenerator random;
	for (int i=0; i<cm.rows(); ++i) {
		for (int j=0; j<cm.cols(); ++j) {
			cm(i,j) = random.randomComplex();
		}
	}
	saveMatrixBin("cm.bin",cm);
	CDMatrix cm2;
	loadMatrixBin("cm.bin",cm2);
	EXPECT_EQ(cm.rows(), cm2.rows());
	EXPECT_EQ(cm.cols(), cm2.cols());
	for (int i=0; i<cm.rows(); ++i) {
		for (int j=0; j<cm.cols(); ++j) {
			dcomplex c1 = cm(i,j);
			dcomplex c2 = cm2(i,j);
			EXPECT_DOUBLE_EQ(c1.real(), c2.real());
			EXPECT_DOUBLE_EQ(c1.imag(), c2.imag());
		}
	}
}


TEST(TestIOSpeed, DISABLED_Normal) {
	CDMatrix cm = CDMatrix::Random(10000,10000);
	saveMatrixBin("cm.bin",cm);
	CDMatrix cm2;
	loadMatrixBin("cm.bin",cm2);
	EXPECT_TRUE(true);
}


TEST(TestIOSpeed, DISABLED_Tmpfs) {
	CDMatrix cm = CDMatrix::Random(10000,10000);
	saveMatrixBin("/dev/shm/cm.bin",cm);
	CDMatrix cm2;
	loadMatrixBin("/dev/shm/cm.bin",cm2);
	EXPECT_TRUE(true);
}
