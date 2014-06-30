/*
 * direct_calculation.h
 *
 *  Created on: Jan 10, 2014
 *      Author: pxiang
 */

#ifndef DIRECT_CALCULATION_H_
#define DIRECT_CALCULATION_H_

#include "../Utility/types.h"
#include "../basis_set/basis.h"
#include "../formMatrix/formMatrix.h"


void formAllBasisSets(LatticeShape& lattice, IMatrix& basisIndex,
		std::vector<Basis>& basisSets);

void formHamiltonianMatrix(LatticeShape& lattice, DMatrix& hamiltonian, IMatrix& basisIndex,
		                   std::vector<Basis>& basisSets);

void obtainEigenVectors(DMatrix& hamiltonian, DVector& eigenValues, DMatrix& eigenVectors);

void numeratorHelper(LatticeShape& lattice, Basis& bra, Basis& ket, IMatrix& basisIndex,
		             DMatrix& eigenVectors, CDArray& numerator);

void denominatorHelper(dcomplex z, DVector& eigenValues, CDArray& oneOverDenominator);

dcomplex greenFuncHelper(CDArray& numerator, CDArray& oneOverDenominator);

void greenFunc_direct(LatticeShape& lattice, Basis& bra, Basis& ket, std::vector<dcomplex >& zList,
		                   std::vector<dcomplex>& gfList);

void densityOfState_direct(LatticeShape& lattice, Basis& basis, std::vector<dcomplex >& zList,
		                   std::vector<double>& dosList);

#endif /* DIRECT_CALCULATION_H_ */
