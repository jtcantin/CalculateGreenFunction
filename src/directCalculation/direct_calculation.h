/*
 * direct_calculation.h
 *
 *  Created on: Jan 10, 2014
 *      Author: pxiang
 */

#ifndef DIRECT_CALCULATION_H_
#define DIRECT_CALCULATION_H_

#include "../Utility/types.h"
#include "../Utility/misc.h"
#include "../basis_set/basis.h"
#include "../formMatrix/formMatrix.h"
#include "../IO/binaryIO.h"
#include "../IO/textIO.h"
#include "../IO/MatrixIO.h"


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

/**
 * calculate all Green's function and save them into files
 */
void calculateAllGreenFunc_direct(LatticeShape& lattice,  Basis& initialSites,
		                          std::vector<dcomplex >& zList,
                                  std::vector<std::string>& fileList);

/**
 * calculate the density of state at all possible sites
 * density_of_state = -Im(<basis | G(z) | basis>)/Pi
 *
 * before calling this, call
 * generateIndexMatrix(lattice1D);
 * setInteractions(lattice1D, interactionData);
 */
void densityOfStateAll_direct(LatticeShape& lattice, std::vector<dcomplex >& zList,
		                   std::vector<std::string>& fileList);

#endif /* DIRECT_CALCULATION_H_ */
