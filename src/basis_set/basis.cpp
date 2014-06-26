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
 * position (x,y). This is a mapping from (x, y) to a positive integer nth.
 *
 * Given a basis (1D or 2D), this function return the corresponding lattice index.
 *
 */
void getLatticeIndex(LatticeShape& lattice, Basis& basis, int &site1, int &site2) {
	switch (lattice.getDim()) {
	case 1:
		site1 = basis[0];
		site2 = basis[1];
		break;
	case 2:
		int xmax = lattice.getXmax();
		site1 = basis[0] + basis[1]*(xmax + 1);
		site2 = basis[2] + basis[3]*(xmax + 1);
		break;
	case 3:
	default:
		std::cout << "Not supported!" << std::endl;
	}
}



// obtain the ''distance'' between two basis,
// for example, in 1D the distance between (1,2) and (2,3) is 3-1=2
// in 2D, the distance between (1,2,3,4) and (1,2,5,6) is sqrt((5-3)^2 + (6-4)^2)
// if the distance is negative, it means that the two basis doesn't share
// a common site. In this case, it is meaningless to calculate the distance
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
//for 1D case: VtoG[K, nth] = G(x1,x2) or basis(x1,x2)
//for 2D case: VtoG[K, nth] = G(x1, y1, x2, y2) or basis(x1,y1,x2,y2)
// each element of VtoG is a basis type (an array of size 2 or 4)
Basis **VtoG;

// to store the dimension of V = ( .., G(i-1,j+1), G(i,j), ...)
std::vector<int> DimsOfV;

//for 1D case: IndexMatrix[x1][x2]=nth  ==>  G(x1,x2) is the nth item in vector V_{x1+x2}
//for 2D case: IndexMatrix[i][j] = nth with i = x1*(xmax+1) + y1, j = x2*(xmax + 1) + y2
//              ==> G(x1, y1, x2, y2) is the nth item in vector V_{x1+y1+x2+y2}
IMatrix IndexMatrix;



void printNeighbors(Neighbors& neighbors) {
	std::cout << "Print neighbors: " << std::endl;
	for (int i=0; i<neighbors.size(); ++i) {
		int dim = neighbors[i].getDim();
		if (dim==1) {
			std::cout << "x1 = " << neighbors[i][0]
			          << "\t\t x2 = " << neighbors[i][0] << std::endl;
		}
	}
}




/**
 * Given a two-particle state, this function generate its neighbors that
 * at some distance. The distance can take both positive and negative values.
 * If distance < 0, we call the neighbors left neighbors. If distance > 0,
 * we call the neighbors right neighbors.
 */
void generateNeighbors(Basis basis, int distance, LatticeShape& lattice,
        Neighbors& neighbors) {

	int dim = lattice.getDim();
	//neighbors.clear();

	switch (dim) {
		case 1:
			// generate Neighbors for the 1D case
			int x1 = basis[0];
			int x2 = basis[1];
			int xmax = lattice.getXmax();

			// generate left neighbors
			if (distance < 0) {
				int x1_left = x1 + distance;
				int x2_left = x2 + distance;
				// left neighbor of first particle
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


// see the above documentation for VtoG, DimsOfV, and Index
void generateIndexMatrix(LatticeShape& lattice) {
	int dim = lattice.getDim();
	switch (dim) {
		// for the 1D case
		case 1:
			int nth;
			int xmax = lattice.getXmax();
			int nsite = xmax + 1; //total number of sites
			IndexMatrix = IMatrix(nsite, nsite);
			// allocate space for VtoG
			//VtoG = PairMatrix(xmax+xmax, xmax+1);
			VtoG = new Basis*[xmax+xmax];
			DimsOfV.resize(xmax+xmax); //K=xsum is from 1 to xmax+xmax-1
			// set all dimension to be zero initially
			for (int i=0; i<DimsOfV.size(); ++i) {
				DimsOfV[i] = 0;
			}

			// find out the dimension of each V
			for (int xsum = 1; xsum <= xmax + xmax-1; ++xsum) {
				nth = 0;
				for (int i=0; i<= xsum/2; ++i) {
					int j = xsum - i; // i+j = xsum
					if (j>i && j<=xmax) {
						IndexMatrix(i,j) = nth; // g(i,j) is the nth (zero-based) item in V_{i+j} = V_K
						nth++;
					}
				}
				DimsOfV[xsum] = nth;
				VtoG[xsum] = new Basis[nth+1]; //index from 0 to nth
			}

			// store each V as a row in VtoG
			for (int xsum = 1; xsum <= xmax + xmax-1; ++xsum) {
				nth = 0;
				for (int i=0; i<= xsum/2; ++i) {
					int j = xsum - i; // i+j = xsum
					if (j>i && j<=xmax) {
						Basis basis(i,j);
						VtoG[xsum][nth] = basis; //createBasis(i,j); // mapping from (K, nth) to (i,j)
						nth++;
					}
				}
			}

			break;
		// for the 2D case
		case 2:
		case 3:
			std::cout<<"Dimensions not supported!"<<std::endl;
			break;
	}
}







