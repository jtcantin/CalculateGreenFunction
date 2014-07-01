/*
 * comparison_test.cpp
 *
 *  Created on: Jun 30, 2014
 *      Author: pxiang
 */
#include "../directCalculation/direct_calculation.h"
#include "../recursiveCalculation/recursiveCalculation.h"
#include "gtest/gtest.h"
#include "../Utility/misc.h"

/**
 * test the results from recursive and direct calculations
 */

TEST(ComparisonTest, DISABLED_CheckDOS) {
	LatticeShape lattice1D(1);
	int xmax = 100;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);
	int size = 201;
	//calculate the indexMatrix
	generateIndexMatrix(lattice1D);

	for (int maxDistance=1; maxDistance<=5; ++maxDistance) {
		InteractionData interactionData = {0.0,5.0,15.0,false,
				                           false,false,maxDistance,230};
		/*********** Recursive calculation *****************/
		int zsize = 101;
		std::vector<dcomplex> zList(zsize);
		std::vector<double> zRealList = linspace(-50,50,zsize);
		for (int i=0; i<zList.size(); ++i) {
			zList[i].real() = zRealList[i];
			zList[i].imag() = 0.1;
		}
		std::vector<double> rhoList;

		calculateDensityOfState(lattice1D, initialSites,
				                 interactionData,
				                 zList, rhoList);
		std::string file1 = "rho_xmax100_recursive_maxDistance_"
				            + itos(maxDistance)+".txt";
		save_two_arrays(file1, zRealList, rhoList);

		/*********** Direct calculation *****************/
		// no need to set up the interaction matrix
		// because it is already done in the recursive calculation
		// setInteractions(lattice1D, interactionData);

		densityOfState_direct(lattice1D, initialSites, zList, rhoList);
		std::string file2 = "rho_xmax100_direct_maxDistance_"
				            + itos(maxDistance)+".txt";
		save_two_arrays(file2, zRealList, rhoList);
	}

	EXPECT_TRUE(true);
}



TEST(ComparisonTest, CheckOffDiagonalOrdered) {
	LatticeShape lattice1D(1);
	int xmax = 100;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);
	Basis finalSites(10, 71);

	//calculate the indexMatrix
	generateIndexMatrix(lattice1D);

	int maxDistance = 5;

	InteractionData interactionData = {0.0,5.0,15.0,false,
			false,false,maxDistance,230};
	/*********** Recursive calculation *****************/
	int zsize = 51;
	std::vector<dcomplex> zList(zsize);
	std::vector<double> zRealList = linspace(-50,50,zsize);
	for (int i=0; i<zList.size(); ++i) {
		zList[i].real() = zRealList[i];
		zList[i].imag() = 0.1;
	}
	std::vector<dcomplex> gfList;
	std::vector<double> gf_real;
	std::vector<double> gf_imag;

	calculateGreenFunc(lattice1D, finalSites, initialSites,
			           interactionData, zList, gfList);

	for (int i=0; i<gfList.size(); ++i) {
		gf_real.push_back(gfList[i].real());
		gf_imag.push_back(gfList[i].imag());
	}

	std::string file1 = "GF_minus40_20_recursive_real.txt";
	save_two_arrays(file1, zRealList, gf_real);
	file1 = "GF_minus40_20_recursive_imag.txt";
	save_two_arrays(file1, zRealList, gf_imag);

	/*********** Direct calculation *****************/
	// no need to set up the interaction matrix
	// because it is already done in the recursive calculation
	// setInteractions(lattice1D, interactionData);

	gf_real.clear();
	gf_imag.clear();
	greenFunc_direct(lattice1D, finalSites, initialSites, zList, gfList);

	for (int i=0; i<gfList.size(); ++i) {
		gf_real.push_back(gfList[i].real());
		gf_imag.push_back(gfList[i].imag());
	}

	file1 = "GF_minus40_20_real.txt";
	save_two_arrays(file1, zRealList, gf_real);
	file1 = "GF_minus40_20_imag.txt";
	save_two_arrays(file1, zRealList, gf_imag);


	EXPECT_TRUE(true);
}



TEST(ComparisonTest, CheckOffDiagonalDisordered) {
	LatticeShape lattice1D(1);
	int xmax = 100;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);
	Basis finalSites(33, 90);

	//calculate the indexMatrix
	generateIndexMatrix(lattice1D);

	int maxDistance = 5;

	InteractionData interactionData = {1.0,5.0,15.0,true,
			true,true,maxDistance,230};
	/*********** Recursive calculation *****************/
	int zsize = 51;
	std::vector<dcomplex> zList(zsize);
	std::vector<double> zRealList = linspace(-50,50,zsize);
	for (int i=0; i<zList.size(); ++i) {
		zList[i].real() = zRealList[i];
		zList[i].imag() = 0.1;
	}
	std::vector<dcomplex> gfList;
	std::vector<double> gf_real;
	std::vector<double> gf_imag;

	calculateGreenFunc(lattice1D, finalSites, initialSites,
			           interactionData, zList, gfList);

	for (int i=0; i<gfList.size(); ++i) {
		gf_real.push_back(gfList[i].real());
		gf_imag.push_back(gfList[i].imag());
	}

	std::string file1 = "GF_minus17_39_recursive_real.txt";
	save_two_arrays(file1, zRealList, gf_real);
	file1 = "GF_minus17_39_recursive_imag.txt";
	save_two_arrays(file1, zRealList, gf_imag);

	/*********** Direct calculation *****************/
	// no need to set up the interaction matrix
	// because it is already done in the recursive calculation
	// setInteractions(lattice1D, interactionData);

	gf_real.clear();
	gf_imag.clear();
	greenFunc_direct(lattice1D, finalSites, initialSites, zList, gfList);

	for (int i=0; i<gfList.size(); ++i) {
		gf_real.push_back(gfList[i].real());
		gf_imag.push_back(gfList[i].imag());
	}

	file1 = "GF_minus17_39_real.txt";
	save_two_arrays(file1, zRealList, gf_real);
	file1 = "GF_minus17_39_imag.txt";
	save_two_arrays(file1, zRealList, gf_imag);


	EXPECT_TRUE(true);
}
