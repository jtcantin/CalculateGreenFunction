/*
 * recursiveCalculation.h
 *
 *  Created on: Jun 24, 2014
 *      Author: pxiang
 */

#ifndef RECURSIVECALCULATION_H_
#define RECURSIVECALCULATION_H_

#include "../Utility/types.h"
#include "../Utility/misc.h"
#include "../basis_set/basis.h"
#include "../IO/binaryIO.h"
#include "../formMatrix/formMatrix.h"

typedef struct {
	int KLeftStart;
	int KLeftStop;
	int KCenter;
	int KRightStop;
	int KRightStart;
	int maxDistance; //step size of the Recursion
	int Csize; // the size of the constant vector C in the quation for V_{KCenter}
	int indexForNonzero; // the index for the only nonzero element in vector C
} RecursionData;


void deleteMatrixFiles(std::string files);

void solveDenseLinearEqs(CDMatrix& A, CDMatrix& B, CDMatrix& X);



/**
 * find out which V_{K} the basis belongs to
 */
int findCorrespondingVK(LatticeShape& lattice, int maxDistance, Basis& basis);

/**
 * find out the index for a basis in V_{K}
 */
int getBasisIndexInVK(LatticeShape& lattice, int K, Basis& basis);

/**
 * 	calculate the indexMatrix and set up the interaction matrix
 */
void setUpIndexInteractions(LatticeShape& lattice,
		InteractionData& interactionData);

void setUpRecursion(LatticeShape& lattice, InteractionData& interactionData,
		            Basis& initialSites, RecursionData& rd);


void fromRightToCenter(RecursionData& recursionData,
		dcomplex z, CDMatrix& AKRightStop, bool saveAMatrices=true);

void fromLeftToCenter(RecursionData& recursionData,
		dcomplex z, CDMatrix& ATildeKLeftStop, bool saveAMatrices=true);

void solveVKCenter(RecursionData& recursionData, dcomplex z,
		           CDMatrix& ATildeKLeftStop, CDMatrix& AKRightStop,
		           CDMatrix& VKCenter);

void calculateDensityOfState(LatticeShape& lattice, Basis& initialSites,
		                      InteractionData& interactionData,
		                      const std::vector<dcomplex>& zList,
		                      std::vector<double>& rhoList);

/**
 * calculate density of state at all sites
 *
 * before calling it, you have to call setUpIndexInteractions(lattice, interactionData)
 */
void calculateDensityOfStateAll(LatticeShape& lattice,
		                      InteractionData& interactionData,
		                      const std::vector<dcomplex>& zList,
		                      std::vector<std::string>& fileList);


/**
 * calculate a matrix element of the Green function
 * <final_sites| G(z) |initial_sites> for a list of z values
 */
void calculateGreenFunc(LatticeShape& lattice, Basis& finalSites, Basis& initialSites,
		                InteractionData& interactionData,
                        const std::vector<dcomplex>& zList,
                        std::vector<dcomplex>& gfList);

/**
 * calculate all the matrix elements of the Green function and save them into a text file
 *
 */
void calculateAllGreenFunc(LatticeShape& lattice,  Basis& initialSites,
		                InteractionData& interactionData, std::vector<dcomplex> zList,
                        std::vector< std::string > fileList);


/*
 * extract the matrix element G(n, m, initial_sites) from files stored in disk
 *
 * lattice --- contains the information about the size and shape of the crystal
 * maxDistance --- the range of dipole-dipole interaction
 *                 (in the unit of lattice constant)
 */

dcomplex extractMatrixElement(int n, int m, LatticeShape& lattice, int maxDistance);

#endif /* RECURSIVECALCULATION_HPP_ */
