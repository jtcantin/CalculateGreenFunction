/*
 * direct_calculation_test.cpp
 *
 *  Created on: Jan 12, 2014
 *      Author: pxiang
 */

#include "direct_calculation.h"
#include "gtest/gtest.h"
#include "../Utility/misc.h"

TEST(DirectCalculationTest, CheckDOS) {
	LatticeShape lattice1D(1);
	int xmax = 100;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	InteractionData interactionData = {0.0,5.0,15.0,false,false,false,2,230,true,true};
	// calculate the indexMatrix and set up the interaction matrix
	generateIndexMatrix(lattice1D);
	setInteractions(lattice1D, interactionData);

	Basis initialSites(xmax/2,xmax/2+1);

	int size = 201;
	std::vector<dcomplex > zList(size);
	std::vector<double> zRealList = linspace(-50,50,size);

	for (int i=0; i<zList.size(); ++i) {
		zList[i] =dcomplex( zRealList[i], 0.1);
	}

	std::vector<double> rhoList;
	densityOfState_direct(lattice1D, initialSites, zList, rhoList);
	save_two_arrays("rho_vs_energy_direct.txt", zRealList, rhoList);
	EXPECT_TRUE(true);
}



TEST(DirectCalculationTest, DISABLED_CheckOffDiagonal) {
	LatticeShape lattice1D(1);
	int xmax = 100;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	InteractionData interactionData = {1.0,1.0,1.0,false,false,false,1,230,true,true};
	// calculate the indexMatrix and set up the interaction matrix
	generateIndexMatrix(lattice1D);
	setInteractions(lattice1D, interactionData);

	int n1i = lattice1D.getXmax()/2;
	int n2i = n1i + 1;
	Basis initialSites(n1i, n2i);

	std::vector<dcomplex > zList(49);
	std::vector<double> zRealList = linspace(-24,24,49);

	for (int i=0; i<zList.size(); ++i) {
		zList[i] =dcomplex( zRealList[i], 0.1);
	}

	std::vector<dcomplex> rhoList;
	int n1f = n1i-17;
	int n2f = n2i + 39;
	Basis finalSites(n1f, n2f);
	greenFunc_direct(lattice1D, finalSites, initialSites, zList, rhoList);
	save_two_arrays("GF_minus17_39.txt", zRealList, rhoList);
	EXPECT_EQ(2,2);
}
