/*
 * basis.cpp
 *
 *  Created on: Jun 13, 2014
 *      Author: pxiang
 */
#ifndef BASIS_H_
#define BASIS_H_

#include <utility>      // std::pair, std::make_pair
#include <string>       // std::string
#include <iostream>     // std::cout
#include <vector>
#include <stdlib.h>
#include <inttypes.h>
#include "../Utility/types.h"
#include "../Utility/misc.h"


class Basis {
public:
	Basis(){
		dim = 0;
		coordinates = NULL;
	}

	// create the basis for 1D cases
	Basis(int x1, int x2) {
		dim = 1;
		coordinates = new int[2];
		coordinates[0]=x1;
		coordinates[1]=x2;
	}

	// create the basis for 2D cases
	Basis(int x1, int y1, int x2, int y2) {
		dim = 2;
		coordinates = new int[4];
		coordinates[0]=x1;
		coordinates[1]=y1;
		coordinates[2]=x2;
		coordinates[3]=y2;
	}

	// copy constructor
	Basis(const Basis& other){
		dim = other.dim;
		int num = 2*other.dim;
		coordinates = new int[num];
	    memcpy(coordinates, other.coordinates, num*sizeof(int));
	}

	// assignment operator
	Basis& operator= (const Basis& other){
		// if the basis has already assign space to the array
		if (coordinates!=NULL) {
			delete [] coordinates;
		}

	    if (this != &other) {
	    	dim = other.dim;
	    	int num = 2*other.dim;
	    	coordinates = new int[num];
	        memcpy(coordinates, other.coordinates, num*sizeof(int));
	    }
	    return *this;
	}

	//destructor
	~Basis() {
		if (coordinates!=NULL) {
			delete [] coordinates;
		}
	}

	// [] operator
	const int& operator[](int i) const {
		return coordinates[i];
	}

	// [] operator for assignment
	int& operator[](int i) {
		return coordinates[i];
	}

	int getDim() {
		return dim;
	}

private:
	int *coordinates;
	int dim;
};

typedef std::vector< Basis > Neighbors;

/*************************** global variables ***********************/
//for 1D case: VtoG[K, nth] = G(x1,x2) or basis(x1,x2)
//for 2D case: VtoG[K, nth] = G(x1, y1, x2, y2) or basis(x1,y1,x2,y2)
// each element of VtoG is a basis type (an array of size 2 or 4)
extern Basis **VtoG;

// to store the dimension of V = ( .., G(i-1,j+1), G(i,j), ...)
extern std::vector<int> DimsOfV;

//for 1D case: Index[x1][x2]=nth  ==>  G(x1,x2) is the nth item in vector V_{x1+x2}
//for 2D case: Index[i][j] = nth with i = x1*xmax + y1, j = x2*xmax + y2
//              ==> G(x1, y1, x2, y2) is the nth item in vector V_{x1+y1+x2+y2}
extern IMatrix IndexMatrix;

double distance(Basis& b1, Basis& b2);

void printNeighbors(Neighbors& neighbors);

// a class to represent the boundary of a lattice
// you need to set up the boundary before using it
class LatticeShape {
public:

//	LatticeShape() {
//		dimension = 0;
//	}

	// use a vector of size 1, 2, 3 to initialize it
	// [xmax, ymax, zmax]
	// note xmin = ymin = zmin = 0
	LatticeShape(int dim_):dimension(0) {
		if (dim_!=1 && dim_!=2 && dim_!=3) {
			std::cout<< "Allowed sizes of the intializer vector are 1, 2, 3"<<std::endl;
			std::cout<< "The initializer size is " << std::endl;
			std::exit(-1);
		}

		dimension = dim_;
		array.resize(dimension);
	}



	int getXmax() {
		if (dimension<1) {
			std::cout << "No xmax" << std::endl;
			std::exit(-1);
		}
		return array[0];
	}



	int getYmax() {
		if (dimension<2) {
			std::cout << "No ymax" << std::endl;
			std::exit(-1);
		}
		return array[1];
	}


	int getZmax() {
		if (dimension<3) {
			std::cout << "No zmax" << std::endl;
			std::exit(-1);
		}
		return array[2];
	}

	void setXmax(int xmax) {
		if (dimension<1) {
			std::cout << "No xmax" << std::endl;
			std::exit(-1);
		}
		array[0] = xmax;
	}

	void setYmax(int ymax) {
		if (dimension<2) {
			std::cout << "No ymax" << std::endl;
			std::exit(-1);
		}
		array[1] = ymax;
	}

	void setZmax(int zmax) {
		if (dimension<3) {
			std::cout << "No zmax" << std::endl;
			std::exit(-1);
		}
		array[2] = zmax;
	}


	// to indicate whether this boundary is for 1D, 2D or 3D
	int getDim() {return dimension;}

	// if the dimension is reset, you need to set xmax, ymax, zmax again
	void setDim(int dim_) {
		dimension = dim_;
		array.clear();
		array.resize(dimension);
	}

private:
	int dimension;
	std::vector<int> array; // to store [ xmax,  [ymax,  [zmax] ]]
};

void getLatticeIndex(LatticeShape& lattice, Basis& basis, int &site1, int &site2);

void generateNeighbors(Basis basis, int distance, LatticeShape& lattice,
        Neighbors& neighbors);

// see the above documentation for VtoG, DimsOfV, and Index
void generateIndexMatrix(LatticeShape& lattice);






#endif /* BASIS_H_ */
