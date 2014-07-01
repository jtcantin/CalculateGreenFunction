/*
 * recursiveCalculation.h
 *
 *  Created on: Jun 24, 2014
 *      Author: pxiang
 */

#ifndef RECURSIVECALCULATION_H_
#define RECURSIVECALCULATION_H_

#include "../Utility/types.h"
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
 * calculate all the matrix elements of the Green function
 * <final_sites| G(z) |initial_sites> for a list of z values
 */
void calculateGreenFunc(LatticeShape& lattice, Basis& finalSites, Basis& initialSites,
		                InteractionData& interactionData,
                        const std::vector<dcomplex>& zList,
                        std::vector<dcomplex>& gfList);

#endif /* RECURSIVECALCULATION_HPP_ */
