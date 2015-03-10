/*
 * textIO_test.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: pxiang
 */

#include "gtest/gtest.h"
#include "textIO.h"

TEST(TextIO, DISABLED_TwoByTwoMatirx) {
	CDMatrix cm(2,2);
	cm(0,0) = dcomplex(1.0, 0.5);

	cm(0,1) = dcomplex(1.3, 0.3);

	cm(1,0) = dcomplex(1.89, 0.35);

	cm(1,1) = dcomplex(1.07899, 0.1135);

	saveMatrixText("cm.txt",cm, true);
	CDMatrix cm2;
	loadMatrixText("cm.txt",cm2);
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


TEST(TextIO, DISABLED_TwoByOneMatirx) {
	CDMatrix cm(2,1);
	cm(0,0) = dcomplex(1.0, 0.5);
	cm(1,0) = dcomplex(1.89, 0.35);
	saveMatrixText("cm.txt",cm, true);
	CDMatrix cm2;
	loadMatrixText("cm.txt",cm2);
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


TEST(TextIO, DISABLED_OneByTwoMatirx) {
	CDMatrix cm(1,2);
	cm(0,0) = dcomplex(1.0, 0.5);
	cm(0,1) = dcomplex(1.3, 0.3);
	saveMatrixText("cm.txt",cm, true);
	CDMatrix cm2;
	loadMatrixText("cm.txt",cm2);
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


TEST(TextIO, DISABLED_OneByOneMatirx) {
	CDMatrix cm(1,1);
	cm(0,0) = dcomplex(1.0, 0.5);
	saveMatrixText("cm.txt",cm, true);
	CDMatrix cm2;
	loadMatrixText("cm.txt",cm2);
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


TEST(TextIO, DISABLED_ThousandByThousandMatirx) {
	CDMatrix cm(1000,1000);
	RandomNumberGenerator random;
	for (int i=0; i<cm.rows(); ++i) {
		for (int j=0; j<cm.cols(); ++j) {
			cm(i,j) = random.randomComplex();
		}
	}
	saveMatrixText("cm.txt",cm, true);
	CDMatrix cm2;
	loadMatrixText("cm.txt",cm2);
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



