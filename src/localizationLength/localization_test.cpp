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
#include "../IO/readInput.h"


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
	int xmax;
	/* read xmax */
	InputVariable xmax_input("int", xstr(xmax), &xmax);
	xmax_input.read("input.txt");
	std::cout << "xmax is " << xmax << std::endl;

	lattice1D.setXmax(xmax); //xsite = xmax + 1
	int n = xmax/2;
	int a = 1;
	int nPlusa = n+a;
	Basis initialSites(n, nPlusa);

	/* read max distance */
	int maxDistance;
	InputVariable maxDistance_input("int", xstr(maxDistance), &maxDistance);
	maxDistance_input.read("input.txt");
	std::cout << "maxDistance is " << maxDistance << std::endl;

	// only random onsite energy
	unsigned seed = 567;
	InteractionData interactionData = {1.0,1.0,1.0,true,
			false,false,maxDistance,seed};

	/* update seed */
	InputVariable inputVar("unsigned", xstr(interactionData.seed),
			                &(interactionData.seed));
	inputVar.read("input.txt");
	std::cout << "The new seed is " << interactionData.seed << std::endl;

	setUpIndexInteractions(lattice1D, interactionData);

	std::vector< dcomplex > zList;
	double E = 0.5;
	InputVariable E_input("double", xstr(E), &E);
	E_input.read("input.txt");
	std::cout << "E is " << E << std::endl;

	double eta = 0.01;
	InputVariable eta_input("double", xstr(eta), &eta);
	eta_input.read("input.txt");
	std::cout << "eta is " << eta << std::endl;

	zList.push_back(dcomplex(E, eta));

	std::vector<std::string> fileList;
	fileList.push_back("GF_random_onsite_energy_recursive_E_0.bin");

	calculateAllGreenFunc(lattice1D,  initialSites,
			              interactionData, zList,
	                       fileList);
	// extract Green's function from the binary file
	CDMatrix gf;
	loadMatrix(fileList[0], gf);
	std::ofstream myFile;

	std::string output = "localization_length.txt";

	myFile.open(output.c_str());
	int m = n;
	int mPlusa = m + a;
	for (;mPlusa<=xmax;++m, ++mPlusa ) {
		dcomplex element = gf(m, mPlusa);
		double localization_length = -1.0/( log(std::abs( element ))/abs(n-m) );
		myFile << abs(n-m) << "\t"<< localization_length << std::endl;
	}
	myFile.close();
	EXPECT_TRUE(true);
}




