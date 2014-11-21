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

//#include <boost/shared_ptr.hpp> //include boost library

/**
 * Basis is used to represent the two-particle basis
 *
 * 1D:  (x1, x2)
 * 2D:  (x1, y1; x2, y2)
 * 3D:  (x1, y1, z1; x2, y2, z2)
 */
class Basis {
public:
	// default constructor
	Basis(){
		dim = 0;
		coordinates = NULL;
	}

	// create a 1D basis
	Basis(int x1, int x2) {
		dim = 1;
		coordinates = new int[2];
		coordinates[0]=x1;
		coordinates[1]=x2;
	}

	// create a 2D basis
	Basis(int x1, int y1, int x2, int y2) {
		dim = 2;
		coordinates = new int[4];
		coordinates[0]=x1;
		coordinates[1]=y1;
		coordinates[2]=x2;
		coordinates[3]=y2;
	}

	// create a 3D basis
//	Basis(int x1, int y1, int z1, int x2, int y2, int z2) {
//
//	}

	// copy constructor
	Basis(const Basis& other){
		dim = other.dim;
		int num = 2*other.dim;
		coordinates = new int[num];
//		for (int i=0; i<num; i++) {
//			coordinates[i] = other.coordinates[i];
//		}
	    memcpy(coordinates, other.coordinates, num*sizeof(int));
	}

	// assignment operator
	Basis& operator= (const Basis& other){
		// if the basis has already assigned space for the array
		if (coordinates!=NULL) {
			delete [] coordinates;
		}

	    if (this != &other) {
	    	dim = other.dim;
	    	int num = 2*other.dim;
	    	coordinates = new int[num];
//			for (int i=0; i<num; i++) {
//				coordinates[i] = other.coordinates[i];
//			}
	        memcpy(coordinates, other.coordinates, num*sizeof(int));
	    }
	    return *this;
	}

	//destructor
	~Basis() {
		if (coordinates!=NULL) {
			delete [] coordinates;
			//std::cout << "Deallocate the memory of the basis" << std::endl;
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

	// obtain the sum of all coordinates
	int getSum() {
		int result = 0;
		switch (dim) {
		case 1:
			result = coordinates[0] + coordinates[1];
			break;
		case 2:
			result = coordinates[0] + coordinates[1]
			         + coordinates[2] + coordinates[3];
			break;
		case 3:
			break;
		}
		return result;
	}

	bool operator==(Basis& other) {
		bool result=false;
		switch (dim) {
		case 1:
			result = (coordinates[0]== other.coordinates[0])
			         &&(coordinates[1]== other.coordinates[1]);
			break;
		case 2:
			result = (coordinates[0]== other.coordinates[0])
	                 &&(coordinates[1]== other.coordinates[1])
	                 &&(coordinates[2]== other.coordinates[2])
	                 &&(coordinates[3]== other.coordinates[3]);
			break;
		case 3:
			break;
		}
		return result;
	}

private:
	int *coordinates;
	int dim;
};


typedef std::vector< Basis > Neighbors;

/*************************** global variables ***********************/
extern std::vector< std::vector< Basis > > VtoG;

extern std::vector<int> DimsOfV;

extern IMatrix IndexMatrix;
/********************************************************************/

double distance(Basis& b1, Basis& b2);

void printNeighbors(Neighbors& neighbors);



/**
 * A class that describes a lattice. It includes the information about the
 * dimension and the boundaries.
 *
 * you need to set up the boundaries (xmax, ymax, zmax) before using it
 * note xmin = ymin = zmin = 0 by default
 */
class LatticeShape {
public:

	// set up the dimension of the lattice
	LatticeShape(int dim_):dimension(0) {
		if (dim_!=1 && dim_!=2 && dim_!=3) {
			std::cout<< "The input dimension is " << dim_ << std::endl;
			std::cout<< "Allowed dimensions are 1, 2, 3"<<std::endl;
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
			std::cout << "Dimension<2, no Ymax" << std::endl;
			std::exit(-1);
		}
		return array[1];
	}


	int getZmax() {
		if (dimension<3) {
			std::cout << "Dimension<3, no Zmax" << std::endl;
			std::exit(-1);
		}
		return array[2];
	}

	void setXmax(int xmax) {
		if (dimension<1) {
			std::cout << "Dimension<1, no Xmax" << std::endl;
			std::exit(-1);
		}
		array[0] = xmax;
	}

	void setYmax(int ymax) {
		if (dimension<2) {
			std::cout << "Dimension<2, no Ymax" << std::endl;
			std::exit(-1);
		}
		array[1] = ymax;
	}

	void setZmax(int zmax) {
		if (dimension<3) {
			std::cout << "Dimension<3, no Zmax" << std::endl;
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

	~LatticeShape() {
		array.clear();
	}

private:
	int dimension;
	std::vector<int> array; // to store [ xmax,  [ymax,  [zmax] ]]
};

void getLatticeIndex(LatticeShape& lattice, Basis& basis, int &site1, int &site2);

void generateNeighbors(Basis basis, int distance, LatticeShape& lattice,
        Neighbors& neighbors);

void generateIndexMatrix(LatticeShape& lattice);






#endif /* BASIS_H_ */
