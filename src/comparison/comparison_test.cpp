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
 * check the diagonal elements of the Green operator
 * from recursive and direct calculations
 */

TEST(ComparisonTest, CheckDOS) {
	LatticeShape lattice1D(1);
	int xmax = 51;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+10);

	//calculate the indexMatrix
	generateIndexMatrix(lattice1D);

	for (int maxDistance=1; maxDistance<=10; ++maxDistance) {
		/********** Set the stage for the calculations**********/
		InteractionData interactionData = {0.0,5.0,15.0,false,
				                           false,false,maxDistance,230,true, true};
		setUpIndexInteractions(lattice1D, interactionData);


		/*********** Recursive calculation *****************/
		int zsize = 101;
		std::vector<dcomplex> zList(zsize);
		std::vector<double> zRealList = linspace(-50,50,zsize);
		for (int i=0; i<zList.size(); ++i) {
			zList[i].real() = zRealList[i];
			zList[i].imag() = 0.1;
		}
		std::vector<double> rhoList_recur;
		calculateDensityOfState(lattice1D, initialSites,
				                 interactionData,
				                 zList, rhoList_recur);
		std::string file1 = "rho_xmax"+itos(xmax)+"_recursive_maxDistance_"
				            + itos(maxDistance)+".txt";
		save_two_arrays(file1, zRealList, rhoList_recur);


		/*********** Direct calculation *****************/
		std::vector<double> rhoList_direct;
		densityOfState_direct(lattice1D, initialSites, zList, rhoList_direct);
		std::string file2 = "rho_xmax"+itos(xmax)+"_direct_maxDistance_"
				            + itos(maxDistance)+".txt";
		save_two_arrays(file2, zRealList, rhoList_direct);

		EXPECT_EQ(rhoList_recur.size(), rhoList_direct.size());
		double abs_error = 1.e-9;
		for (int j=0; j<rhoList_recur.size(); ++j) {
			EXPECT_NEAR(rhoList_recur[j], rhoList_direct[j], abs_error);
		}
	}

	//EXPECT_TRUE(true);
}


/**
 * check the off-diagonal elements of the Green's operator
 */
TEST(ComparisonTest, CheckOffDiagonalOrdered) {
	LatticeShape lattice1D(1);
	int xmax = 51;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	int ni = xmax/2;
	int mi = ni + 1;
	int n_difference = -15;
	int m_difference = 10;
	int nf = ni + n_difference;
	int mf = mi + m_difference;

	Basis initialSites(ni,mi);
	Basis finalSites(nf, mf);



	/********** Set the stage for the calculations**********/
	int maxDistance = 30;
	InteractionData interactionData = {0.0,5.0,5.0,false,
			false,false,maxDistance,230,true,true};
	setUpIndexInteractions(lattice1D, interactionData);


	/*********** Recursive calculation *****************/
	int zsize = 51;
	std::vector<dcomplex> zList(zsize);
	std::vector<double> zRealList = linspace(-25,25,zsize);
	for (int i=0; i<zList.size(); ++i) {
		zList[i].real() = zRealList[i];
		zList[i].imag() = 0.01;
	}
	std::vector<dcomplex> gfList;
	std::vector<double> gf_real_recur;
	std::vector<double> gf_imag_recur;

	calculateGreenFunc(lattice1D, finalSites, initialSites,
			           interactionData, zList, gfList);

	for (int i=0; i<gfList.size(); ++i) {
		gf_real_recur.push_back(gfList[i].real());
		gf_imag_recur.push_back(gfList[i].imag());
	}

	std::string file1 = "GF_"+ itos(nf)+"," + itos(mf) + "," +
			               itos(ni)+"," + itos(mi) + "_recursive_real.txt";
	save_two_arrays(file1, zRealList, gf_real_recur);

	file1 = "GF_"+ itos(nf)+"," + itos(mf) + "," +
            itos(ni)+"," + itos(mi) + "_recursive_imag.txt";
	save_two_arrays(file1, zRealList, gf_imag_recur);



	/*********** Direct calculation *****************/
	// no need to set up the interaction matrix
	// because it is already done in the recursive calculation
	// setInteractions(lattice1D, interactionData);

	std::vector<double> gf_real_direct;
	std::vector<double> gf_imag_direct;
	greenFunc_direct(lattice1D, finalSites, initialSites, zList, gfList);

	for (int i=0; i<gfList.size(); ++i) {
		gf_real_direct.push_back(gfList[i].real());
		gf_imag_direct.push_back(gfList[i].imag());
	}

	file1 = "GF_"+ itos(nf)+"," + itos(mf) + "," +
            itos(ni)+"," + itos(mi) + "_direct_real.txt";
	save_two_arrays(file1, zRealList, gf_real_direct);

	file1 = "GF_"+ itos(nf)+"," + itos(mf) + "," +
            itos(ni)+"," + itos(mi) + "_direct_imag.txt";
	save_two_arrays(file1, zRealList, gf_imag_direct);


	// test equality
	EXPECT_EQ(gf_real_recur.size(), gf_real_direct.size());
	EXPECT_EQ(gf_imag_recur.size(), gf_imag_direct.size());

	double abs_error = 1.e-9;
	for (int j=0; j<gf_real_recur.size(); ++j) {
		EXPECT_NEAR(gf_real_recur[j], gf_real_direct[j], abs_error);
		EXPECT_NEAR(gf_imag_recur[j], gf_imag_direct[j], abs_error);
	}

}



/**
 * check all matrix elements of the Green's operator
 */
TEST(ComparisonTest, CheckAllMatrixElements) {
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








TEST(ComparisonTest, DISABLED_RandomOnSiteEnergy) {
	LatticeShape lattice1D(1);
	int xmax = 101;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);

	int maxDistance = 10;

	// only random onsite energy
	InteractionData interactionData = {1.0,1.0,1.0,true,
			false,false,maxDistance,567,true,true};
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




TEST(ComparisonTest, DISABLED_RandomHopping) {
	LatticeShape lattice1D(1);
	int xmax = 101;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);

	int maxDistance = 10;

	// only random hopping interaction
	InteractionData interactionData = {1.0,1.0,1.0,false,
			true,false,maxDistance,1234,true,true};
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


TEST(ComparisonTest, DISABLED_NoDisorderNoDyn) {
	LatticeShape lattice1D(1);
	int xmax = 101;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);

	int maxDistance = 10;

	// no disorder and no dynamic interaction
	InteractionData interactionData = {1.0,1.0,0.0,false,
			false,false,maxDistance,230,true,true};
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



TEST(ComparisonTest, DISABLED_RandomHoppingLargeDyn) {
	LatticeShape lattice1D(1);
	int xmax = 101;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);

	int maxDistance = 10;

	// no disorder and no dynamic interaction
	InteractionData interactionData = {1.0,1.0,4.0,false,
			true,false,maxDistance,1234,true,true};
	setUpIndexInteractions(lattice1D, interactionData);

	std::vector< dcomplex > zList;
	zList.push_back(dcomplex(0, 0.01));
	zList.push_back(dcomplex(1, 0.01));
	std::vector<std::string> fileList;

	/********** Recursive calculations of Green's function******************/
	fileList.push_back("GF_random_hopping_large_dyn_recursive_E_0.txt");
	fileList.push_back("GF_random_hopping_large_dyn_recursive_E_1.txt");
	calculateAllGreenFunc(lattice1D,  initialSites,
			              interactionData, zList,
	                       fileList);

	/********** Direct calculations of Green's function*********************/
	fileList.clear();
	fileList.push_back("GF_random_hopping_large_dyn_direct_E_0.txt");
	fileList.push_back("GF_random_hopping_large_dyn_direct_E_1.txt");
	calculateAllGreenFunc_direct(lattice1D,  initialSites,
			                      zList, fileList);

	/********** Recursive calculations of density of state ******************/
	fileList.clear();
	fileList.push_back("DOS_random_hopping_large_dyn_recursive_E_0.txt");
	fileList.push_back("DOS_random_hopping_large_dyn_recursive_E_1.txt");
	calculateDensityOfStateAll(lattice1D, interactionData, zList, fileList);

	/********** Direct calculations of density of state *********************/
	fileList.clear();
	fileList.push_back("DOS_random_hopping_large_dyn_direct_E_0.txt");
	fileList.push_back("DOS_random_hopping_large_dyn_direct_E_1.txt");
	densityOfStateAll_direct(lattice1D, zList, fileList);

	EXPECT_TRUE(true);
}



TEST(ComparisonTest, DISABLED_NoDisorderLargeDyn) {
	LatticeShape lattice1D(1);
	int xmax = 101;
	lattice1D.setXmax(xmax); //xsite = xmax + 1
	Basis initialSites(xmax/2,xmax/2+1);

	int maxDistance = 10;

	// no disorder and no dynamic interaction
	InteractionData interactionData = {1.0,1.0,4.0,false,
			false,false,maxDistance,230,true,true};
	setUpIndexInteractions(lattice1D, interactionData);

	std::vector< dcomplex > zList;
	zList.push_back(dcomplex(0, 0.01));
	zList.push_back(dcomplex(1, 0.01));
	std::vector<std::string> fileList;

	/********** Recursive calculations of Green's function******************/
	fileList.push_back("GF_ordered_large_dyn_recursive_E_0.txt");
	fileList.push_back("GF_ordered_large_dyn_recursive_E_1.txt");
	calculateAllGreenFunc(lattice1D,  initialSites,
			              interactionData, zList,
	                       fileList);

	/********** Direct calculations of Green's function*********************/
	fileList.clear();
	fileList.push_back("GF_ordered_large_dyn_direct_E_0.txt");
	fileList.push_back("GF_ordered_large_dyn_direct_E_1.txt");
	calculateAllGreenFunc_direct(lattice1D,  initialSites,
			                      zList, fileList);

	/********** Recursive calculations of density of state ******************/
	fileList.clear();
	fileList.push_back("DOS_ordered_large_dyn_recursive_E_0.txt");
	fileList.push_back("DOS_ordered_large_dyn_recursive_E_1.txt");
	calculateDensityOfStateAll(lattice1D, interactionData, zList, fileList);

	/********** Direct calculations of density of state *********************/
	fileList.clear();
	fileList.push_back("DOS_ordered_large_dyn_direct_E_0.txt");
	fileList.push_back("DOS_ordered_large_dyn_direct_E_1.txt");
	densityOfStateAll_direct(lattice1D, zList, fileList);

	EXPECT_TRUE(true);
}




