/*
 * initTesting_test.cpp
 *
 *  Created on: Jun 30, 2014
 *      Author: pxiang
 *  Edits:
 *  2015/03/11 - Joshua T Cantin - Generated this file from comparison_test.cpp, using the "CheckAllMatrixElements" test
 */
#include "../directCalculation/direct_calculation.h"
#include "../recursiveCalculation/recursiveCalculation.h"
#include "gtest/gtest.h"
#include "../Utility/misc.h"

/**
 * check all matrix elements of the Green's operator
 */
TEST(ComparisonTest, DISABLED_CheckAllMatrixElements) {
	LatticeShape lattice1D(1);
	int xmax = 101;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	int ni = xmax/2;
	int mi = ni + 1;
	Basis initialSites(ni,mi);

	/********** Set the stage for the calculations**********/
	int maxDistance = 10;
	InteractionData interactionData = {1.0,1.0,1.0,true,
			false,false,maxDistance,230,true,true};
	setUpIndexInteractions(lattice1D, interactionData);


	/*********** Recursive calculation *****************/
	std::vector<dcomplex> zList;
	dcomplex z;
	z.real() = 0.0;
	z.imag() = 0.01;
	zList.push_back(z);

	std::vector< std::string > fileList_recur;
	fileList_recur.push_back("GF_recur.bin");

	calculateAllGreenFunc(lattice1D,  initialSites,
			                   interactionData, zList, fileList_recur);

	// extract Green's function from the file
	CDMatrix gf_recur;
	loadMatrix(fileList_recur[0], gf_recur);



	/*********** Direct calculation *****************/
	// no need to set up the interaction matrix
	// because it is already done in the recursive calculation
	// setInteractions(lattice1D, interactionData);

	std::vector< std::string > fileList_direct;
	fileList_direct.push_back("GF_direct.bin");

	calculateAllGreenFunc_direct(lattice1D,  initialSites,
			                          zList, fileList_direct);

	// extract Green's function from the file
	CDMatrix gf_direct;
	loadMatrix(fileList_direct[0], gf_direct);

	// test equality
	EXPECT_EQ(gf_recur.rows(), gf_direct.rows());
	EXPECT_EQ(gf_recur.cols(), gf_direct.cols());

	double abs_error = 1.e-9;
	for (int i=0; i<gf_recur.rows(); ++i) {
		for (int j=0; j<gf_recur.cols(); ++j) {
			dcomplex element_recur = gf_recur(i,j);
			dcomplex element_direct = gf_direct(i,j);
			EXPECT_NEAR(element_recur.real(), element_direct.real(), abs_error);
			EXPECT_NEAR(element_recur.imag(), element_direct.imag(), abs_error);
		}
	}

}
