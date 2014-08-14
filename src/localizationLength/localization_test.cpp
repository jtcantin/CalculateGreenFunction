/*
 * localization_test.cpp
 *
 *  Created on: Aug 14, 2014
 *      Author: pxiang
 */

#include "../directCalculation/direct_calculation.h"
#include "../recursiveCalculation/recursiveCalculation.h"
#include "gtest/gtest.h"
#include "../Utility/misc.h"


/**
 * 1D crystal
 * Randomize the on-site energies to be between 0 an 1
 * Fix the hopping terms to be 1 for the nearest neighbours and decaying as 1/r^3
 * Fix the dynamical terms to be 1 for the nearest neighbours and decaying as 1/r^3
 *
 * calculate
 * 1/two_particle_localization_length
 * = - lim_{|n-m|--> infty} << ln|<m, m+a|G(z)|n, n+a>| >> / |n-m|
 * where << ... >> represents the average over different configurations of disorder
 *
 * In the current case, we let a = 1 and
 * we want z.real to be at the center of the band
 */
TEST(LocalizationTest, RandomOnSiteEnergy) {
	LatticeShape lattice1D(1);
	int xmax = 401;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	int n = xmax/2;
	int a = 1;
	int nPlusa = n+a;
	Basis initialSites(n, nPlusa);

	int maxDistance = 10;

	// only random onsite energy
	unsigned seed = 567;
	InteractionData interactionData = {1.0,1.0,1.0,true,
			false,false,maxDistance,seed};
	setUpIndexInteractions(lattice1D, interactionData);

	std::vector< dcomplex > zList;
	zList.push_back(dcomplex(0, 0.01));

	std::vector<std::string> fileList;
	fileList.push_back("GF_random_onsite_energy_recursive_E_0.bin");

	calculateAllGreenFunc(lattice1D,  initialSites,
			              interactionData, zList,
	                       fileList);
	// extract Green's function from the binary file
	CDMatrix gf;
	loadMatrix(fileList[0], gf);
	int m = n;
	int mPlusa = m + a;
	for (;mPlusa<=xmax;++m, ++mPlusa ) {
		dcomplex element = gf(m, mPlusa);
		double localization_length = -1.0/( log(std::abs( element ))/abs(n-m) );
		std::cout << localization_length << std::endl;
	}

	EXPECT_TRUE(true);
}




