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

void deleteMatrixFiles(std::string files);

void solveDenseLinearEqs(CDMatrix& A, CDMatrix& B, CDMatrix& X);

void setUpRecursion(LatticeShape& lattice, Parameters& parameters,
		            Basis& initialSites, int maxDistance,
		            int& KLeftStart, int& KLeftStop, int& KCritial,
		            int& KRightStop, int& KRightStart);

void fromRightToCenter(int KRightStart, int KRightStop, int maxDistance,
		dcomplex z, CDMatrix& AKRightStop, bool saveAMatrices=true);

void fromLeftToCenter(int KLeftStart, int KLeftStop, int maxDistance,
		dcomplex z, CDMatrix& ATildeKLeftStop, bool saveAMatrices=true);


#endif /* RECURSIVECALCULATION_HPP_ */
