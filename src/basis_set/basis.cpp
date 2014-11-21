/*
 * basis.cpp
 *
 *  Created on: Jun 18, 2014
 *      Author: pxiang
 */
#include "basis.h"


/**
 * In a lattice, we give each site an index. For example, in 1D, we can use the
 * index numbers 0 ... xmax to denote each sites. In 2D, we still want to use a
 * single number to denote each lattice site, so we assign a number to each
 * position (x,y). This is a mapping from (x, y) to a nonnegative integer nth.
 *
 * 1D case (x_for_particle1, x_for_particle1):
 *
 * x------x------x------x------...------x------x------x
 * 0      1      2      3      ...   xmax-2  xmax-1  xmax
 *
 *
 *
 * 2D case (x_for_p1, y_for_p1; x_for_p2, y_for_p2):
 *
 *	 x------x------x------x------...------x------x------x
 *	 |      |      |      |      ...      |      |      |
 *	 x------x------x------x------...------x------x------x
 *	 |      |      |      |      ...      |      |      |
 *	 x------x------x------x------...------x------x------x
 *	 |      |      |      |      ...      |      |      |
 *	 .      .      .      .      ...      .      .      .
 *	 .      .      .      .      ...      .      .      .
 *	 .      .      .      .      ...      .      .      .
 *	 |      |      |      |      ...      |      |      |
 *	 x------x------x------x------...------x------x------x
 *(xmax+1)                                          (2*xmax+1)
 *	 |      |      |      |      ...      |      |      |
 *	 x------x------x------x------...------x------x------x
 *	(0)    (1)    (2)            ...                  (xmax)
 *
 * Given a two-particle basis (1D or 2D), this function returns the corresponding
 * lattice indexes for the two sites
 *
 */
void getLatticeIndex(LatticeShape& lattice, Basis& basis, int &site1, int &site2) {
	switch (lattice.getDim()) {
	case 1:
		site1 = basis[0];
		site2 = basis[1];
		break;
	case 2:
		{
			int xmax = lattice.getXmax();
			site1 = basis[0] + basis[1]*(xmax + 1);
			site2 = basis[2] + basis[3]*(xmax + 1);
			break;
		}
	case 3:
	default:
		std::cout << "Dimension >= 3, not supported!" << std::endl;
	}
}



/** This function calculate the ''distance'' between two basis
 *
 * For example, in 1D the distance between (1,2) and (2,3) is 3-1=2.
 * in 2D, the distance between (1,2,3,4) and (1,2,5,6) is sqrt((5-3)^2 + (6-4)^2)
 *
 * Note if the two basis don't share a common site, it is meaniningless to
 * calculate the distance. In this case, the function return -1, indicating
 * that the function cannot handle the current case.
 *
 */
double distance(Basis& b1, Basis& b2) {
	// if the two basis are of different dimension, it is meaningless to calculate
	// the distance
	if (b1.getDim()!=b2.getDim()) {
		return -1.0;
	}

	double result;
	switch (b1.getDim()) {
	case 1: // 1D case
		if (b1[0]==b2[0]) {
			result = std::fabs(b1[1]-b2[1]);
		} else if (b1[1]==b2[1]) {
			result = std::fabs(b1[0]-b2[0]);
		} else {
			result = -1.0;
		}
		break;
	case 2:
	case 3:
		std::cout << "2D and 3D are not supported yet!" << std::endl;
		exit(-1);
	}

	return result;
}


/*************************** global variables ***********************/
/**
 * VtoG is a Matrix that is used to convert a vector V to the corresponding
 * Green's function
 *
 * for 1D case: VtoG[K, nth] gives the value of G(x1,x2) that is nth element of V
 * for 2D case: VtoG[K, nth] gives the value of G(x1, y1, x2, y2) that is nth
 *              element of V
 */
std::vector< std::vector< Basis > > VtoG;





/**
 * We group Green's functions into vectors V_K = ( .., G(i-1,j+1), G(i,j), ...)
 * where K = i + j.
 *
 * DimsOfV is used to record the size of various V_K:
 * size_of_V_K = DimsOfV[K]  (Kmin <= K <= Kmax)
 */
std::vector<int> DimsOfV;


/**
 * IndexMatrix maps a pairs of site index (site1, site2) to an number nth, which
 * has the following meaning:
 *
 * 1D case:
 *          G(x1, x2) --> 1D basis(x1, x2)
 *
 *                    --> (site1=x1, site2=x2) or it can obtained from function
 *                        getLatticeIndex(lattice, basis, site1, site2)
 *
 *                    --> IndexMatrix(site1, site2) = nth
 *
 *                    ==> G(x1, x2) is the nth item in the corresponding
 *                        vector V_{x1+x2}
 *
 * 2D case:
 *          G(x1, y1, x2, y2) --> 2D basis(x1, y1, x2, y2)
 *
 *                            --> (site1, site2) can be obtained from function
 *                                getLatticeIndex(lattice, basis, site1, site2)
 *
 *                            --> IndexMatrix(site1, site2) = nth
 *
 *                            ==> G(x1, y1, x2, y2) is the nth item in the
 *                                corresponding vector V_{x1+y1+x2+y2}
 */
IMatrix IndexMatrix;
/***********************************************************************/


/**
 * print an array of neighbors for testing purposes
 */
void printNeighbors(Neighbors& neighbors) {
	std::cout << "Print neighbors: " << std::endl;
	for (int i=0; i<neighbors.size(); ++i) {
		int dim = neighbors[i].getDim();

		if (dim==1) {
			std::cout << "x1 = " << neighbors[i][0]
			          << "\t\t x2 = " << neighbors[i][0] << std::endl;
		}

		if (dim==2) {
			// placeholder
		}

		if (dim==3) {
			// placeholder
		}

	}
}




/**
 * Given a two-particle basis, this function generates its neighbors that
 * at some distance. The distance can take both positive and negative values.
 * If distance < 0, we call the neighbors left neighbors. If distance > 0,
 * we call the neighbors right neighbors.
 */
void generateNeighbors(Basis basis, int distance, LatticeShape& lattice,
        Neighbors& neighbors) {

	int dim = lattice.getDim();
	neighbors.clear();

	switch (dim) {
		case 1: {
			// generate Neighbors for the 1D case
			// assuming x1 < x2
			int x1 = basis[0];
			int x2 = basis[1];
			int xmax = lattice.getXmax();

			// generate left neighbors
			if (distance < 0) {
				int x1_left = x1 + distance;
				int x2_left = x2 + distance;
				// left neighbor of first particle exists
				if (x1_left>=0) {
					Basis p1(min(x1_left, x2), max(x1_left, x2));
					neighbors.push_back(p1);
				}

				// left neighbor of second particle
				if (x2_left>=0 && x2_left != x1) {
					Basis p3(min(x1, x2_left),max(x1, x2_left));
					neighbors.push_back(p3);
				}
			}

			// generate right neighbors
			if ( distance > 0 ) {
				int x1_right = x1 + distance;
				int x2_right = x2 + distance;
				// right neighbor of first particle
				if (x1_right<=xmax && x1_right!=x2) {
					Basis p2(min(x1_right, x2), max(x1_right, x2));
					neighbors.push_back(p2);
				}


				// right neighbor of second particle
				if (x2_right<= xmax) {
					Basis p4(min(x1, x2_right), max(x1, x2_right));
					neighbors.push_back(p4);
				}
			}
			break;
		}
		case 2:
			// generate Neighbors for the 2D case
			std::cout<<"Dimensions not supported!"<<std::endl;
			break;
		case 3:

			break;
		default:
			std::cout<<"Dimensions not supported!"<<std::endl;

	}


}


/**
 * Generate the global variables: VtoG, DimsOfV, IndexMatrix
 *
 * This must be called before any calculation of the Green's function. When the
 * lattice changes, this function should also be called.
 */
void generateIndexMatrix(LatticeShape& lattice) {
	// reset those global variables
	VtoG.clear();
	DimsOfV.clear();
	IndexMatrix.resize(0,0);

	int dim = lattice.getDim();
	switch (dim) {
		// for the 1D case
		case 1:{
			int nth;
			int xmax = lattice.getXmax();
			//total number of sites, note the index starts from 0
			int nsite = xmax + 1;
			IndexMatrix = IMatrix(nsite, nsite);

			/**
			 * allocate space for VtoG
			 *
			 * for each K in [Kmin, Kmax] where Kmin=0+1=1,
			 * there should be a corresponding vector
			 *
			 * In order that VtoG[K] is for K, the size of VtoG
			 * is made a little larger than necessary. That is,
			 * we need to have VtoG[0], VtoG[Kmin=1], ..., VtoG[Kmax],
			 * so the size of VtoG is set to be Kmax-0+1
			 *
			 */
			int Kmin = 1;
			int Kmax = xmax + xmax - 1; //the largest value for K
			VtoG.resize(Kmax+1);
			/**
			 * similarly, we need to have DimsOfV[0], DimsOfV[Kmin=1],...,
			 * ..., DimsOfV[Kmax], so the size of DimsOfV is also Kmax+1
			 */
			DimsOfV.resize(Kmax+1);


			// set all dimension to be zero initially
			for (int i=0; i<DimsOfV.size(); ++i) {
				DimsOfV[i] = 0;
			}

			// find out the dimension of each V_K
			for (int K = Kmin; K <= Kmax; ++K) {
				nth = 0;
				for (int site1=0; site1<= K/2; ++site1) {//assuming site1 < site2
					int site2 = K - site1; // site1 + site2 = K
					if (site2 > site1 && site2<=xmax) {
						// g(site1, site2) is the nth (zero-based) item in V_{i+j} = V_K
						IndexMatrix(site1,site2) = nth;
						nth++;
					}
				}
				DimsOfV[K] = nth;
				//VtoG[K].resize(nth+1);  //index from 0 to nth

				//reserve some memory for VtoG[K]
				//to avoid reallocate and deallocate memory during inserting items
				VtoG[K].reserve(nth+1);
			}

			// store each V_K as a row with index K in VtoG
			for (int K = 1; K <= Kmax; ++K) {
				for (int site1=0; site1<= K/2; ++site1) {
					int site2 = K - site1;
					if (site2 > site1 && site2<=xmax) {
						// mapping from (K, nth) to (site1,site2)
						VtoG[K].push_back(Basis(site1, site2));
					}
				}
			}

			break;
		}
		// for the 2D case
		case 2:
		case 3:
			std::cout<<"Dimensions not supported!"<<std::endl;
			break;
	}
}







