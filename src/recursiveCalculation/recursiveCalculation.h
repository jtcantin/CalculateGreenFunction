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
} RecursionData;


void deleteMatrixFiles(std::string files);

void solveDenseLinearEqs(CDMatrix& A, CDMatrix& B, CDMatrix& X);

void setUpRecursion(LatticeShape& lattice, Parameters& parameters,
		            Basis& initialSites, int maxDistance,RecursionData& rd);


void fromRightToCenter(int KRightStart, int KRightStop, int maxDistance,
		dcomplex z, CDMatrix& AKRightStop, bool saveAMatrices=true);

void fromLeftToCenter(int KLeftStart, int KLeftStop, int maxDistance,
		dcomplex z, CDMatrix& ATildeKLeftStop, bool saveAMatrices=true);

void solveVKCenter(CDMatrix& VKCenter, CDMatrix& ATildeKLeftStop, CDMatrix& AKRightStop,
		           int KCenter, LatticeShape& lattice, Basis& initialSites,
		           int maxDistance, dcomplex z);

#endif /* RECURSIVECALCULATION_HPP_ */
