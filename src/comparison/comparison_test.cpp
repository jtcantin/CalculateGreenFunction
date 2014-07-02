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
	Basis initialSites(xmax/2,xmax/2+10);
	int size = 201;
	//calculate the indexMatrix
	generateIndexMatrix(lattice1D);

	for (int maxDistance=1; maxDistance<=5; ++maxDistance) {
		InteractionData interactionData = {0.0,5.0,15.0,false,
				                           false,false,maxDistance,230};
		setUpIndexInteractions(lattice1D, interactionData);
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

		densityOfState_direct(lattice1D, initialSites, zList, rhoList);
		std::string file2 = "rho_xmax100_direct_maxDistance_"
				            + itos(maxDistance)+".txt";
		save_two_arrays(file2, zRealList, rhoList);
	}

	EXPECT_TRUE(true);
}



TEST(ComparisonTest, DISABLED_CheckOffDiagonalOrdered) {
	LatticeShape lattice1D(1);
	int xmax = 100;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);
	Basis finalSites(10, 71);



	int maxDistance = 5;

	InteractionData interactionData = {0.0,5.0,15.0,false,
			false,false,maxDistance,230};
	setUpIndexInteractions(lattice1D, interactionData);

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



TEST(ComparisonTest, DISABLED_CheckOffDiagonalDisordered) {
	LatticeShape lattice1D(1);
	int xmax = 100;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);
	Basis finalSites(33, 90);


	int maxDistance = 5;

	InteractionData interactionData = {1.0,5.0,15.0,true,
			true,true,maxDistance,230};
	setUpIndexInteractions(lattice1D, interactionData);

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




TEST(ComparisonTest, RandomOnSiteEnergy) {
	LatticeShape lattice1D(1);
	int xmax = 50;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);

	int maxDistance = 5;

	// only random onsite energy
	InteractionData interactionData = {1.0,1.0,1.0,false,
			false,true,maxDistance,230};
	setUpIndexInteractions(lattice1D, interactionData);

	std::vector< dcomplex > zList;
	zList.push_back(dcomplex(0, 0.01));
	zList.push_back(dcomplex(1, 0.01));
	std::vector<std::string> fileList;

	/********** Recursive calculations of Green's function******************/
	fileList.push_back("GF_random_onsite_energy_recursive_E_0.txt");
	fileList.push_back("GF_random_onsite_energy_recursive_E_1.txt");
	calculateAllGreenFunc(lattice1D,  initialSites,
			              interactionData, zList,
	                       fileList);

	/********** Direct calculations of Green's function*********************/
	fileList.clear();
	fileList.push_back("GF_random_onsite_energy_direct_E_0.txt");
	fileList.push_back("GF_random_onsite_energy_direct_E_1.txt");
	calculateAllGreenFunc_direct(lattice1D,  initialSites,
			                      zList, fileList);

	/********** Recursive calculations of density of state ******************/
	fileList.clear();
	fileList.push_back("DOS_random_onsite_energy_recursive_E_0.txt");
	fileList.push_back("DOS_random_onsite_energy_recursive_E_1.txt");
	calculateDensityOfStateAll(lattice1D, interactionData, zList, fileList);

	/********** Direct calculations of density of state *********************/
	fileList.clear();
	fileList.push_back("DOS_random_onsite_energy_direct_E_0.txt");
	fileList.push_back("DOS_random_onsite_energy_direct_E_1.txt");
	densityOfStateAll_direct(lattice1D, zList, fileList);

	EXPECT_TRUE(true);
}




TEST(ComparisonTest, RandomHopping) {
	LatticeShape lattice1D(1);
	int xmax = 50;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);

	int maxDistance = 5;

	// only random hopping interaction
	InteractionData interactionData = {1.0,1.0,1.0,true,
			false,false,maxDistance,230};
	setUpIndexInteractions(lattice1D, interactionData);

	std::vector< dcomplex > zList;
	zList.push_back(dcomplex(0, 0.01));
	zList.push_back(dcomplex(1, 0.01));
	std::vector<std::string> fileList;

	/********** Recursive calculations of Green's function******************/
	fileList.push_back("GF_random_hopping_recursive_E_0.txt");
	fileList.push_back("GF_random_hopping_recursive_E_1.txt");
	calculateAllGreenFunc(lattice1D,  initialSites,
			              interactionData, zList,
	                       fileList);

	/********** Direct calculations of Green's function*********************/
	fileList.clear();
	fileList.push_back("GF_random_hopping_direct_E_0.txt");
	fileList.push_back("GF_random_hopping_direct_E_1.txt");
	calculateAllGreenFunc_direct(lattice1D,  initialSites,
			                      zList, fileList);

	/********** Recursive calculations of density of state ******************/
	fileList.clear();
	fileList.push_back("DOS_random_hopping_recursive_E_0.txt");
	fileList.push_back("DOS_random_hopping_recursive_E_1.txt");
	calculateDensityOfStateAll(lattice1D, interactionData, zList, fileList);

	/********** Direct calculations of density of state *********************/
	fileList.clear();
	fileList.push_back("DOS_random_hopping_direct_E_0.txt");
	fileList.push_back("DOS_random_hopping_direct_E_1.txt");
	densityOfStateAll_direct(lattice1D, zList, fileList);

	EXPECT_TRUE(true);
}


TEST(ComparisonTest, NoDisorderNoDyn) {
	LatticeShape lattice1D(1);
	int xmax = 50;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);

	int maxDistance = 5;

	// no disorder and no dynamic interaction
	InteractionData interactionData = {1.0,0.0,1.0,false,
			false,false,maxDistance,230};
	setUpIndexInteractions(lattice1D, interactionData);

	std::vector< dcomplex > zList;
	zList.push_back(dcomplex(0, 0.01));
	zList.push_back(dcomplex(1, 0.01));
	std::vector<std::string> fileList;

	/********** Recursive calculations of Green's function******************/
	fileList.push_back("GF_ordered_no_dyn_recursive_E_0.txt");
	fileList.push_back("GF_ordered_no_dyn_recursive_E_1.txt");
	calculateAllGreenFunc(lattice1D,  initialSites,
			              interactionData, zList,
	                       fileList);

	/********** Direct calculations of Green's function*********************/
	fileList.clear();
	fileList.push_back("GF_ordered_no_dyn_direct_E_0.txt");
	fileList.push_back("GF_ordered_no_dyn_direct_E_1.txt");
	calculateAllGreenFunc_direct(lattice1D,  initialSites,
			                      zList, fileList);

	/********** Recursive calculations of density of state ******************/
	fileList.clear();
	fileList.push_back("DOS_ordered_no_dyn_recursive_E_0.txt");
	fileList.push_back("DOS_ordered_no_dyn_recursive_E_1.txt");
	calculateDensityOfStateAll(lattice1D, interactionData, zList, fileList);

	/********** Direct calculations of density of state *********************/
	fileList.clear();
	fileList.push_back("DOS_ordered_no_dyn_direct_E_0.txt");
	fileList.push_back("DOS_ordered_no_dyn_direct_E_1.txt");
	densityOfStateAll_direct(lattice1D, zList, fileList);

	EXPECT_TRUE(true);
}
