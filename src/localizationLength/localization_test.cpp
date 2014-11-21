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
	READ_INPUT("input.txt", "int", xmax);


	lattice1D.setXmax(xmax); //xsite = xmax + 1
	int initialSeperation=1; //default value is 1

	/* set up interaction data */
	double onsiteE, hop, dyn;
	bool randomOnsite, randomHop, randomDyn;
	int maxDistance;
	unsigned seed;
	bool longRangeHop, longRangeDyn;
	READ_INPUT("input.txt", "double", onsiteE);
	READ_INPUT("input.txt", "double", hop);
	READ_INPUT("input.txt", "double", dyn);
	READ_INPUT("input.txt", "bool", randomOnsite);
	READ_INPUT("input.txt", "bool", randomHop);
	READ_INPUT("input.txt", "bool", randomDyn);
	READ_INPUT("input.txt", "int", maxDistance);
	READ_INPUT("input.txt", "unsigned", seed);
	READ_INPUT("input.txt", "bool", longRangeHop);
	READ_INPUT("input.txt", "bool", longRangeDyn);
	READ_INPUT("input.txt", "int", initialSeperation);

	int n = xmax/2;
	int nPlusa = n + initialSeperation;
	Basis initialSites(n, nPlusa);

	InteractionData interactionData = { onsiteE, hop, dyn,
			                            randomOnsite, randomHop, randomDyn,
			                            maxDistance, seed,
	                                    longRangeHop, longRangeDyn};

	setUpIndexInteractions(lattice1D, interactionData);

	std::vector< dcomplex > zList;

	double E;
	READ_INPUT("input.txt", "double", E);


	double eta;
	READ_INPUT("input.txt", "double", eta);


	zList.push_back(dcomplex(E, eta));

	std::vector<std::string> fileList;
	//fileList.push_back("GF.txt");
	fileList.push_back("GF.bin");

	calculateAllGreenFunc(lattice1D,  initialSites,
			              interactionData, zList,
	                       fileList);
	// extract Green's function from the file
	CDMatrix gf;
	loadMatrix(fileList[0], gf);

	std::ofstream myFile;
	std::string output = "a_" + itos(initialSeperation)+".txt";

	myFile.open(output.c_str());
	int m = n;
	int mPlusa = m + initialSeperation;
	for ( ; mPlusa<=xmax; ++m, mPlusa=m+initialSeperation ) {
		dcomplex element = gf(m, mPlusa);
		double G_n_m = std::log( std::abs( element ) );
		myFile << (m-n) << "\t"<< G_n_m << std::endl;
	}
	myFile.close();

	EXPECT_TRUE(true);
}




TEST(LocalizationTest, RandomOnSiteEnergyWithNoDisorderRange) {
	LatticeShape lattice1D(1);
	int xmax;
	READ_INPUT("input.txt", "int", xmax);


	lattice1D.setXmax(xmax); //xsite = xmax + 1
	int initialSeperation=1; //default value is 1

	/* set up interaction data */
	double onsiteE, hop, dyn;
	bool randomOnsite, randomHop, randomDyn;
	int maxDistance;
	unsigned seed;
	bool longRangeHop, longRangeDyn;
	READ_INPUT("input.txt", "double", onsiteE);
	READ_INPUT("input.txt", "double", hop);
	READ_INPUT("input.txt", "double", dyn);
	READ_INPUT("input.txt", "bool", randomOnsite);
	READ_INPUT("input.txt", "bool", randomHop);
	READ_INPUT("input.txt", "bool", randomDyn);
	READ_INPUT("input.txt", "int", maxDistance);
	READ_INPUT("input.txt", "unsigned", seed);
	READ_INPUT("input.txt", "bool", longRangeHop);
	READ_INPUT("input.txt", "bool", longRangeDyn);
	READ_INPUT("input.txt", "int", initialSeperation);

	int n = xmax/2;
	int nPlusa = n + initialSeperation;
	Basis initialSites(n, nPlusa);

	InteractionData interactionData = { onsiteE, hop, dyn,
			                            randomOnsite, randomHop, randomDyn,
			                            maxDistance, seed,
	                                    longRangeHop, longRangeDyn};

	int radius = 200;
	setUpIndexInteractions_test(lattice1D, interactionData, radius);

	std::vector< dcomplex > zList;

	double E;
	READ_INPUT("input.txt", "double", E);


	double eta;
	READ_INPUT("input.txt", "double", eta);


	zList.push_back(dcomplex(E, eta));

	std::vector<std::string> fileList;
	//fileList.push_back("GF.txt");
	fileList.push_back("GF.bin");

	calculateAllGreenFunc(lattice1D,  initialSites,
			              interactionData, zList,
	                       fileList);
	// extract Green's function from the file
	CDMatrix gf;
	loadMatrix(fileList[0], gf);

	std::ofstream myFile;
	std::string output = "size_" + itos(xmax+1)
			              +"_a_" + itos(initialSeperation)+".txt";

	myFile.open(output.c_str());
	int center = xmax/2;
	int start = center - radius;
	int end = center + radius-initialSeperation;

	for ( int i=start; i<=end; ++i ) {
		dcomplex element = gf(i, i+initialSeperation);
		//double G_n_m = std::log( std::abs( element ) );
		myFile << (i-center) << "\t" << (i+initialSeperation-center)
				<<"\t"<< element.real() <<"\t" << element.imag() << std::endl;
	}
	myFile.close();

	EXPECT_TRUE(true);
}



