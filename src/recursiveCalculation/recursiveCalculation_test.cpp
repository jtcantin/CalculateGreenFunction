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
	Parameters parameters = {1.0,1.0,1.0,true,false,true,5,230};
	int maxDistance = 2;
	Basis initialSites(60,61);
	RecursionData rd;

	setUpRecursion(lattice1D,  parameters, initialSites, maxDistance,rd);

	std::cout<<"xmax = " << xmax << std::endl;
	std::cout<<"x1 = " << initialSites[0]
	         <<"\t x2 = " << initialSites[1]<< std::endl;
	std::cout<<"maxDistance = " << maxDistance << std::endl;
	std::cout<<"KLeftStart = " << rd.KLeftStart << std::endl;
	std::cout<<"KLeftStop = " << rd.KLeftStop << std::endl;
	std::cout<<"KCritial = " << rd.KCenter << std::endl;
	std::cout<<"KRightStop = " << rd.KRightStop << std::endl;
	std::cout<<"KRightStart = " << rd.KRightStart << std::endl;
}


TEST(FromRightToCenter, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 120;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Parameters parameters = {1.0,1.0,1.0,true,false,true,5,230};
	int maxDistance = 2;
	Basis initialSites(60,61);
	RecursionData rd;


	setUpRecursion(lattice1D,  parameters, initialSites, maxDistance,rd);

	dcomplex z = dcomplex(10, 0.01);
	CDMatrix AKRightStop;

	deleteMatrixFiles("A[0-9]*.bin");
	fromRightToCenter(rd.KRightStart, rd.KRightStop, maxDistance,
			z, AKRightStop, true);

	EXPECT_TRUE(true);
}


TEST(FromLeftToCenter, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 120;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Parameters parameters = {1.0,1.0,1.0,true,false,true,5,230};
	int maxDistance = 2;
	Basis initialSites(60,61);
	RecursionData rd;

	setUpRecursion(lattice1D,  parameters, initialSites, maxDistance,rd);


	dcomplex z = dcomplex(10, 0.01);
	CDMatrix ATildeKLeftStop;

	deleteMatrixFiles("ATilde[0-9]*.bin");
	fromLeftToCenter(rd.KLeftStart, rd.KLeftStop, maxDistance,
			z, ATildeKLeftStop, true);

	EXPECT_TRUE(true);
}


TEST(SolveVKCenter, RunningOK) {
	LatticeShape lattice1D(1);
	int xmax = 120;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Parameters parameters = {1.0,1.0,1.0,true,false,true,5,230};
	int maxDistance = 5;
	Basis initialSites(60,61);
	RecursionData rd;

	setUpRecursion(lattice1D,  parameters, initialSites, maxDistance,rd);


	dcomplex z = dcomplex(10, 0.01);
	CDMatrix ATildeKLeftStop;

	deleteMatrixFiles("ATilde[0-9]*.bin");
	fromLeftToCenter(rd.KLeftStart, rd.KLeftStop, maxDistance,
			z, ATildeKLeftStop, true);

	CDMatrix AKRightStop;

	deleteMatrixFiles("A[0-9]*.bin");
	fromRightToCenter(rd.KRightStart, rd.KRightStop, maxDistance,
			z, AKRightStop, true);

	CDMatrix VKCenter;
	solveVKCenter(VKCenter, ATildeKLeftStop, AKRightStop,
			       rd.KCenter, lattice1D, initialSites,
			       maxDistance, z);

	EXPECT_TRUE(true);
}

