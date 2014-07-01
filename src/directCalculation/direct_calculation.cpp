/*
 * direct_calculation.cpp
 *
 *  Created on: Jan 10, 2014
 *      Author: pxiang
 */

// calculate the Green's function directly from eigenfunctions


#include "direct_calculation.h"

/**
 * forming all the basis sets
 *
 * A two-particle basis can be generally written as
 * basis( index_for_site1, index_for_site2 ) = basis( i, j )
 * In 1D, a site can be labeled by the corresponding lattice coordinate x.
 * In 2D, a site is specified by the lattice coordinate (x,y), and it can be
 * labeled by index = x + y*(xmax + 1). We have a function:
 * void getLatticeIndex(LatticeShape& lattice, Basis& basis, int &site1, int &site2)
 * which can be used to calculate the index for site1 and site2 given a basis set
 *
 * basisIndex(site1, site2) = nth means that the basis set with (site1, site2) is
 * the nth basis set.
 *
 */
void formAllBasisSets(LatticeShape& lattice, IMatrix& basisIndex,
		std::vector<Basis>& basisSets) {
	switch (lattice.getDim()) {
	case 1:
		// generate all basis sets
		int nbasis = 0; //the total number of basis set
		int xmax = lattice.getXmax();
		// to avoid double counting, we only consider (n1, n2) where n1<n2
		for (int n1=0; n1<=xmax-1; ++n1) {
			for (int n2=n1+1; n2<=xmax; ++n2) {
				nbasis++;
			}
		}

		//std::cout<< "OK before resizing basisIndex" << std::endl;
		// allocate memory for basisSets and basisIndex
		basisIndex.resize(xmax+1, xmax+1);
//		std::cout<< "OK after resizing basisIndex" << std::endl;
		basisSets.reserve(nbasis); // allocate the space needed to store all basis sets
		int nth = 0;
		for (int n1=0; n1<=xmax-1; ++n1) {
			for (int n2=n1+1; n2<=xmax; ++n2) {
				Basis basis(n1, n2);
				basisSets.push_back(basis);

				int site1, site2;
				getLatticeIndex(lattice, basis, site1, site2);
				basisIndex(site1, site2) = nth;
				nth++;


			}
		}
		break;

	case 2:
		break;
	case 3:
		break;
	}

}


void formHamiltonianMatrix(LatticeShape& lattice, DMatrix& hamiltonian, IMatrix& basisIndex,
		                   std::vector<Basis>& basisSets) {
	extern Interaction *pInteraction;
	int size = basisSets.size();
	hamiltonian = DMatrix::Zero(size, size);
	// fill in the matrix column by column  (by default, the storage order is column-major)
	for (int col=0; col<size; ++col) {
		Basis ket = basisSets[col];

		// fill in the diagonal part
		hamiltonian(col, col)= pInteraction->onsiteE(ket) + pInteraction->dyn(ket);

		// fill in the off-diagonal part
		int bra_site1;
		int bra_site2;
		int row;
		// find out which basis set will produce non-interacting matrix with ket

		for (int distance=1; distance<=pInteraction->getMaxDistance(); ++distance) {
			Neighbors neighborsAtSomeDistance;
			// positive distance
			generateNeighbors(ket,distance, lattice, neighborsAtSomeDistance);
			for (int i=0; i<neighborsAtSomeDistance.size(); ++i) {
				Basis bra = neighborsAtSomeDistance[i];
				getLatticeIndex(lattice, bra, bra_site1, bra_site2);
				row = basisIndex(bra_site1, bra_site2);
				hamiltonian(row, col) = pInteraction->hop(bra, ket);
			}

			// negative distance
			generateNeighbors(ket, -distance, lattice, neighborsAtSomeDistance);
			for (int i=0; i<neighborsAtSomeDistance.size(); ++i) {
				Basis bra = neighborsAtSomeDistance[i];
				getLatticeIndex(lattice, bra, bra_site1, bra_site2);
				row = basisIndex(bra_site1, bra_site2);
				hamiltonian(row, col) = pInteraction->hop(bra, ket);
			}

		}


	}

}


void obtainEigenVectors(DMatrix& hamiltonian, DVector& eigenValues, DMatrix& eigenVectors) {
	Eigen::SelfAdjointEigenSolver<DMatrix> eigensolver(hamiltonian);
	if (eigensolver.info() == Eigen::Success) {
		eigenValues = eigensolver.eigenvalues();
		eigenVectors = eigensolver.eigenvectors();
	} else {
		printf("Failed to solve the eigen problem!");
	}
}


/**
 * helper functions for calculating the green function <bra | G(z) | ket>
 *
 * for 1D: <bra| = basis(n1f, n2f)    |ket> = basis(n1i, n1i)
 */
/**************************************************************************************/
void numeratorHelper(LatticeShape& lattice, Basis& bra, Basis& ket, IMatrix& basisIndex,
		             DMatrix& eigenVectors, CDArray& numerator) {
	int n1, n2;
	getLatticeIndex(lattice, bra, n1, n2);
	int bra_index = basisIndex(n1, n2);

	getLatticeIndex(lattice, ket, n1, n2);
	int ket_index = basisIndex(n1, n2);

	DVector wavefunc_f = eigenVectors.row(bra_index);
	DVector wavefunc_i = eigenVectors.row(ket_index);
	DArray tmp = wavefunc_f.array()*wavefunc_i.array();
	numerator = tmp.cast< dcomplex >();
}

void denominatorHelper(dcomplex z, DVector& eigenValues, CDArray& oneOverDenominator) {
	int n = eigenValues.rows();
	oneOverDenominator = CDArray(n);
	for (int i=0; i<n; ++i) {
		oneOverDenominator(i) = 1.0/(z - eigenValues(i));
	}

}

dcomplex greenFuncHelper(CDArray& numerator, CDArray& oneOverDenominator) {
	CDArray tmp = numerator*oneOverDenominator;
	return tmp.sum();
}
/**************************************************************************************/





/**
 * calculating the green function <bra | G(z) | ket>
 */
void greenFunc_direct(LatticeShape& lattice, Basis& bra, Basis& ket, std::vector<dcomplex >& zList,
		                   std::vector<dcomplex>& gfList) {
	IMatrix basisIndex;
	std::vector<Basis> basisSets;
	formAllBasisSets(lattice, basisIndex, basisSets);


	DMatrix hamiltonian;
	formHamiltonianMatrix(lattice, hamiltonian, basisIndex, basisSets);


	DVector eigenValues;
	DMatrix eigenVectors;
	obtainEigenVectors(hamiltonian, eigenValues, eigenVectors);


	CDArray numerator;
	numeratorHelper(lattice, bra, ket, basisIndex, eigenVectors, numerator);

	gfList.clear();
	gfList.reserve(zList.size());
	for (int i=0; i<zList.size(); ++i) {
		dcomplex z = zList[i];
		CDArray oneOverDenominator;
		denominatorHelper(z, eigenValues, oneOverDenominator);
		dcomplex gf = greenFuncHelper(numerator, oneOverDenominator);
		gfList[i] = gf; //-gf.imag()/M_PI;
	}

}



/**
 * calculate the density of state at (site1, site2)
 * density_of_state = -Im(<basis | G(z) | basis>)/Pi
 */
void densityOfState_direct(LatticeShape& lattice, Basis& basis, std::vector<dcomplex >& zList,
		                   std::vector<double>& dosList) {
	IMatrix basisIndex;
	std::vector<Basis> basisSets;
	formAllBasisSets(lattice, basisIndex, basisSets);


	DMatrix hamiltonian;
	formHamiltonianMatrix(lattice, hamiltonian, basisIndex, basisSets);


	DVector eigenValues;
	DMatrix eigenVectors;
	obtainEigenVectors(hamiltonian, eigenValues, eigenVectors);


	CDArray numerator;
	numeratorHelper(lattice, basis, basis, basisIndex, eigenVectors, numerator);

	dosList.clear();
	dosList.reserve(zList.size());
	for (int i=0; i<zList.size(); ++i) {
		dcomplex z = zList[i];
		CDArray oneOverDenominator;
		denominatorHelper(z, eigenValues, oneOverDenominator);
		dcomplex gf = greenFuncHelper(numerator, oneOverDenominator);
		dosList[i] = -gf.imag()/M_PI;
	}

}
