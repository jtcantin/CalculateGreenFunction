/*
 * recursiveCalculation_test.cpp
 *
 *  Created on: Jun 25, 2014
 *      Author: pxiang
 */
#include "gtest/gtest.h"
#include "recursiveCalculation.h"

TEST(SetUpRecursion, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 120;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,5,230};
	int maxDistance = interactionData.maxDistance;
	Basis initialSites(60,61);
	RecursionData rd;

	setUpRecursion(lattice1D,  interactionData, initialSites, rd);

	std::cout<<"xmax = " << xmax << std::endl;
	std::cout<<"x1 = " << initialSites[0]
	         <<"\t x2 = " << initialSites[1]<< std::endl;
	std::cout<<"maxDistance = " << maxDistance << std::endl;
	std::cout<<"KLeftStart = " << rd.KLeftStart << std::endl;
	std::cout<<"KLeftStop = " << rd.KLeftStop << std::endl;
	std::cout<<"KCenter = " << rd.KCenter << std::endl;
	std::cout<<"KRightStop = " << rd.KRightStop << std::endl;
	std::cout<<"KRightStart = " << rd.KRightStart << std::endl;
}


TEST(FromRightToCenter, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 120;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,2,230};
	Basis initialSites(60,61);
	RecursionData recursionData;


	setUpRecursion(lattice1D,  interactionData, initialSites,recursionData);

	dcomplex z = dcomplex(10, 0.01);
	CDMatrix AKRightStop;

	deleteMatrixFiles("A[0-9]*.bin");
	fromRightToCenter(recursionData, z, AKRightStop, true);

	EXPECT_TRUE(true);
}


TEST(FromLeftToCenter, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 120;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,2,230};
	Basis initialSites(60,61);
	RecursionData recursionData;

	setUpRecursion(lattice1D,  interactionData, initialSites, recursionData);


	dcomplex z = dcomplex(10, 0.01);
	CDMatrix ATildeKLeftStop;

	deleteMatrixFiles("ATilde[0-9]*.bin");
	fromLeftToCenter(recursionData, z, ATildeKLeftStop, true);

	EXPECT_TRUE(true);
}


TEST(SolveVKCenter, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 121;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	InteractionData interactionData = {1.0,1.0,1.0,true,false,true,1,230};
	Basis initialSites(60,61);
	RecursionData recursionData;

	setUpRecursion(lattice1D,  interactionData, initialSites,recursionData);
	std::cout<<"KLeftStart = " << recursionData.KLeftStart << std::endl;
	std::cout<<"KLeftStop = " << recursionData.KLeftStop << std::endl;
	std::cout<<"KCenter = " << recursionData.KCenter << std::endl;
	std::cout<<"KRightStop = " << recursionData.KRightStop << std::endl;
	std::cout<<"KRightStart = " << recursionData.KRightStart << std::endl;

	dcomplex z = dcomplex(10, 0.01);
	CDMatrix ATildeKLeftStop;

	deleteMatrixFiles("ATilde[0-9]*.bin");
	fromLeftToCenter(recursionData, z, ATildeKLeftStop, true);

	CDMatrix AKRightStop;

	deleteMatrixFiles("A[0-9]*.bin");
	fromRightToCenter(recursionData, z, AKRightStop, true);

	CDMatrix VKCenter;
	solveVKCenter(recursionData, z, ATildeKLeftStop, AKRightStop, VKCenter);

	EXPECT_TRUE(true);
}


TEST(GeneratingDensityOfStates, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 100;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	InteractionData interactionData = {0.0,5.0,15.0,false,false,false,2,230};
	Basis initialSites(xmax/2, xmax/2 + 1);
	int zsize = 201;
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
	save_two_arrays("rho_vs_energy.txt", zRealList, rhoList);
	EXPECT_TRUE(true);
}
