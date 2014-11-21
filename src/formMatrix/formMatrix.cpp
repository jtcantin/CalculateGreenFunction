/*
 * formMatrix.cpp
 *
 *  Created on: Jun 17, 2014
 *      Author: pxiang
 */

#include "formMatrix.h"

// the following three functions are defined for testing purposes
double tMatrix(Interaction& interaction, int i, int j) {
	return interaction.t(i,j);
}

double dMatrix(Interaction& interaction, int i, int j) {
	return interaction.d(i,j);
}

double eVector(Interaction& interaction, int i) {
	return interaction.e(i);
}



// form all sorts of matrices




Interaction *pInteraction;
LatticeShape *pLattice;


// call this before forming the matrices
void setInteractions(LatticeShape& lattice, InteractionData& interactionData) {
	// destroy the Interaction object pointed by pInteraction
//	if (pLattice!=NULL) {
//		std::cout<<"Deallocate previously created LatticeShape objects"<<std::endl;
//		pLattice->~LatticeShape();
//	}
	pLattice = &lattice;

	// destroy the Interaction object pointed by pInteraction
	if (pInteraction!=NULL) {
		//std::cout<<"Deallocate previously created Interaction objects"<<std::endl;
		pInteraction->~Interaction();
	}
	pInteraction = new Interaction(lattice, interactionData);
}

// only for testing purpose
void setInteractions_test(LatticeShape& lattice, InteractionData& interactionData,
		                  int radius) {
	pLattice = &lattice;

	// destroy the Interaction object pointed by pInteraction
	if (pInteraction!=NULL) {
		//std::cout<<"Deallocate previously created Interaction objects"<<std::endl;
		pInteraction->~Interaction();
	}
	pInteraction = new Interaction(lattice, interactionData);
	pInteraction->setNoDisorderRange(radius);
}


/**
 * obtain the size of the Matrix M_{K, Kp}
 */
void getMSize(int K, int Kp, int& rows, int& cols) {
	extern std::vector<int> DimsOfV;
	int Kmin = 1;
	int Kmax = DimsOfV.size()-1;
	if (K>=Kmin && K<=Kmax && Kp>=Kmin && Kp<=Kmax) {
		rows = DimsOfV[K];
		cols = DimsOfV[Kp];
	} else {
		std::cout<< "ERROR: K="<<K<< " and Kp=" << Kp<<" are out of range!" << std::endl;
		exit(-1);
	}
}


/**
 * obtain the size of the Matrix Z_{K}
 */
void getZSize(int K, int& rows, int& cols) {
	extern std::vector<int> DimsOfV;
	rows = DimsOfV[K];
	cols = DimsOfV[K];
}

/**
 * before calling the following subroutines that form matrices,
 * make sure that pInteraction points to a valid Interaction object
 * and "void generateIndexMatrix(LatticeShape& lattice)" has been
 * called such that VtoG, DimsOfV, and IndexMatrix have values
 */
void formMatrixZ(int K, dcomplex Energy, CDMatrix& ZK) {
	extern std::vector< std::vector< Basis > > VtoG;
	int rows, cols;
	getZSize(K, rows, cols);
	ZK = CDMatrix::Zero(rows, cols);

	// only diagonal elements are nonzero
	for (int i=0; i<rows; ++i) {
			Basis basis =  VtoG[K][i]; // this is a pointer to a 1D or 2D basis
			ZK(i,i) = Energy - pInteraction->onsiteE(basis)
					   - pInteraction->dyn(basis);
	}
}

/**
 * calculate the matrix M_{K, Kp}
 */
void formMatrixM(int K, int Kp, CDMatrix& MKKp) {
	extern std::vector< std::vector< Basis > > VtoG;
	extern IMatrix IndexMatrix;

	int distance = Kp - K;

	int rows, cols;
	getMSize(K, Kp, rows, cols);
	MKKp = CDMatrix::Zero(rows, cols);

	for (int i=0; i<rows; ++i) {
		Neighbors neighbors;
		Basis basis1 = VtoG[K][i];
		int site1, site2;
		getLatticeIndex(*pLattice, basis1, site1, site2);
//		std::cout << "LEFT site1: " << site1 <<"\t site2: "<< site2 <<std::endl;
		int row = IndexMatrix(site1, site2);

		generateNeighbors( basis1, distance, *pLattice, neighbors);

		for (int j=0; j<neighbors.size(); ++j) {
			Basis basis2 = neighbors[j];
			getLatticeIndex(*pLattice, basis2, site1, site2);
//			std::cout << "OK before calling IndexMatrix" << std::endl;
			int col = IndexMatrix(site1, site2);
//			std::cout << "OK after calling IndexMatrix" << std::endl;
//			std::cout << "RIGHT site1: " << site1 <<"\t site2: "<< site2 <<std::endl;
			MKKp(row, col) = pInteraction->hop(basis1, basis2);
//			std::cout << "OK" << std::endl;
		}

		//freeMemory(neighbors);
		//neighbors.clear();
	}
}


/**
 * form the WK matrix
 */
void formMatrixW(int K, dcomplex energy, CDMatrix& WK) {
	extern Interaction *pInteraction;
	int maxDistance = pInteraction->getMaxDistance();
	// find out the size of WK matrix
	int total_rows=0;
	int total_cols=0;
	// go through the blocks that on the diagonal
	// it may happen that the number of blocks is < maxDistance
	extern std::vector<int> DimsOfV;
	int Kmax = DimsOfV.size()-1;
	int numBlock = min(Kmax-K+1, maxDistance); //the number of blocks
	for (int i=0; i<numBlock; ++i) {
		int rows, cols;
		getZSize(K+i, rows, cols);
		total_rows += rows;
		total_cols += cols;
	}

	//assign memory for WK
	WK = CDMatrix::Zero(total_rows, total_cols);

	int row_start = 0;
	int col_start = 0;
	int row_size;
	int col_size;
	for (int block_row=0; block_row<numBlock; ++block_row) {
		col_start = 0; // rewind block col count
		for (int block_col=0; block_col<numBlock; ++block_col) {
			if (block_row == block_col) { // diagonal blocks
				CDMatrix Z;
				formMatrixZ(K+block_row, energy, Z);
				row_size = Z.rows();
				col_size = Z.cols();
				WK.block(row_start, col_start, row_size, col_size) = Z;
			} else { // off-diagonal blocks
				CDMatrix M;
				formMatrixM(K+block_row, K+block_col, M);
				row_size = M.rows();
				col_size = M.cols();
				WK.block(row_start, col_start, row_size, col_size) = -M;
			}
			//std::cout << "block row: " << block_row << " block col:" << block_col << std::endl;
			col_start += col_size;
		}

		row_start += row_size;
	}
}


/**
 * form the Alpha matrix
 */
void formMatrixAlpha(int K, CDMatrix& AlphaK) {
	extern Interaction *pInteraction;
	int maxDistance = pInteraction->getMaxDistance();
	// find out the size of alpha matrix
	int total_rows=0;
	int total_cols=0;

	// it may happen that the number of blocks is < maxDistance
	extern std::vector<int> DimsOfV;
	int Kmax = DimsOfV.size()-1;
//	std::cout<<"Insider formMatrixAlpha, Kmax = " << Kmax << " K=" <<K << std::endl;
	int numBlockRow = min(Kmax-K+1, maxDistance);
	int numBlockCol = maxDistance;

	// go through the last column block by block
	for (int i=0; i<numBlockRow; ++i) {
		int rows, cols;
//		std::cout << "Insider formMatrixAlpha" << std::endl;
		getMSize(K+i, K-1, rows, cols);
		total_rows += rows;
	}

	// go through the first row block by block
	for (int i=0; i<numBlockCol; ++i) {
		int rows, cols;
//		std::cout << "Insider formMatrixAlpha" << std::endl;
		getMSize(K, K-maxDistance+i, rows, cols);
		total_cols += cols;
	}
//	std::cout<<"Size of AlphaK: " << total_rows << "  " << total_cols << std::endl;
//	std::cout << "before initialize AlphaK" << std::endl;

	//assign memory for AlphaK
	AlphaK = CDMatrix::Zero(total_rows, total_cols);
//	std::cout << "after initialize AlphaK\n" << std::endl;

	int row_start = 0;
	int col_start = 0;
	int row_size;
	int col_size;

	for (int block_row=0; block_row<numBlockRow; ++block_row) {
		// rewind block col count
		col_start = 0;
		// find out the starting column index
		for (int i=0; i<block_row; ++i) {
			int rows, cols;
			getMSize(K+i, K+i-maxDistance, rows, cols);
			col_start += cols;
		}

		for (int block_col=block_row; block_col<numBlockCol; ++block_col) {
			// calculate only the upper triangle part is nonzero
//			std::cout << "block row: " << block_row << " block col: "
//					<< block_col << std::endl;
			CDMatrix M;
//			std::cout << "before calling formMatrixM" << std::endl;
			formMatrixM(K+block_row, K+block_col-maxDistance, M);
//			std::cout << "after calling formMatrixM" << std::endl;
			row_size = M.rows();
			col_size = M.cols();
			AlphaK.block(row_start, col_start, row_size, col_size) = M;
			col_start += col_size;
		} // end of inner for loop

		row_start += row_size;
	} // end of outer for loop
}



/**
 * form the Beta matrix
 */
void formMatrixBeta(int K, CDMatrix& BetaK) {
	extern Interaction *pInteraction;
	int maxDistance = pInteraction->getMaxDistance();
	// find out the size of beta matrix
	int total_rows=0;
	int total_cols=0;

	// go through the blocks that on the diagonal
	// it may happen that the number of blocks is < maxDistance
	extern std::vector<int> DimsOfV;
	int Kmax = DimsOfV.size()-1;
//	std::cout << "Kmax = " << Kmax << std::endl;
	int numBlockRow = min(Kmax-K+1, maxDistance);
//	std::cout << "numBlockRow = " << numBlockRow << std::endl;
	int numBlockCol = min(Kmax-(K+maxDistance)+1, maxDistance); //the number of blocks


	// go through the first column of blocks row by row
	for (int i=0; i<numBlockRow; ++i) {
		int rows, cols;
//		std::cout << "Insider formMatrixBeta" << std::endl;
		getMSize(K+i, K+maxDistance, rows, cols);
		total_rows += rows;
	}

	// go through the column of blocks
	for (int i=0; i<numBlockCol; ++i) {
		int rows, cols;
//		std::cout << "Insider formMatrixBeta" << std::endl;
		getMSize(K+maxDistance-1, K+maxDistance+i, rows, cols);
		total_cols += cols;
	}


	//assign memory for BetaK
	BetaK = CDMatrix::Zero(total_rows, total_cols);

	int row_start = 0;
	int col_start = 0;
	int row_size;
	int col_size;

	for (int block_row=0; block_row<numBlockRow; ++block_row) {
		// rewind block col count
		col_start = 0;

		for (int block_col=0; block_col<=min(block_row,numBlockCol-1); ++block_col) {
			// calculate only the upper triangle part is nonzero
//			std::cout << "block row: " << block_row << " block col: "
//					<< block_col << std::endl;
			CDMatrix M;
//			std::cout << "before calling formMatrixM" << std::endl;
			formMatrixM(K+block_row, K+block_col+maxDistance, M);
//			std::cout << "after calling formMatrixM" << std::endl;
			row_size = M.rows();
			col_size = M.cols();
			BetaK.block(row_start, col_start, row_size, col_size) = M;
			col_start += col_size;
		} // end of inner for loop

		row_start += row_size;
	} // end of outer for loop
}
