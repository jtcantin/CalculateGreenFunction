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

TEST(ComparisonTest, CheckDOS) {
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
