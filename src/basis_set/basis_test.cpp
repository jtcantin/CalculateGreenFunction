/*
 * basis_test.cpp
 *
 *  Created on: Jun 13, 2014
 *      Author: pxiang
 */
#include "gtest/gtest.h"
#include "basis.h"


TEST(BasisTest, DISABLED_AllBasisTypes) {
	Basis b1(1,2);
	EXPECT_EQ( b1[0], 1);
	EXPECT_EQ( b1[1], 2);

	Basis b = b1;
	EXPECT_EQ( b[0], 1);
	EXPECT_EQ( b[1], 2);

	Basis b2(1,2,3,4);
	EXPECT_EQ( b2[0], 1);
	EXPECT_EQ( b2[1], 2);
	EXPECT_EQ( b2[2], 3);
	EXPECT_EQ( b2[3], 4);

	b = b2;
	EXPECT_EQ( b2[0], 1);
	EXPECT_EQ( b2[1], 2);
	EXPECT_EQ( b2[2], 3);
	EXPECT_EQ( b2[3], 4);
}


TEST(NeighborsTest, DISABLED_AssignmentAndFree) {
	Neighbors neighbors;
	Basis p(0,1);
	neighbors.push_back(p);
	EXPECT_EQ( neighbors[0][0], p[0] );
	EXPECT_EQ( neighbors[0][1], p[1] );

	p[0] = 3;
	p[1] = 5;
	neighbors.push_back(p);
	EXPECT_EQ( neighbors[1][0], p[0] );
	EXPECT_EQ( neighbors[1][1], p[1] );

	p[0] = 6;
	p[1] = 9;
	neighbors.push_back(p);
	EXPECT_EQ( neighbors[2][0], p[0] );
	EXPECT_EQ( neighbors[2][1], p[1] );

	EXPECT_EQ( neighbors.size(), 3);
}


// test the Boundary class
TEST(BoundaryTest, DISABLED_RightSize) {
	std::vector<int> maxXYZ;
	int xmax = 100;
	int ymax = 101;
	int zmax = 81;

	LatticeShape lattice1D(1);
	lattice1D.setXmax(xmax);
	EXPECT_EQ(lattice1D.getDim(),1);
	EXPECT_EQ(lattice1D.getXmax(),xmax);


	LatticeShape lattice2D(2);
	lattice2D.setXmax(xmax);
	lattice2D.setYmax(ymax);
	EXPECT_EQ(lattice2D.getDim(),2);
	EXPECT_EQ(lattice1D.getXmax(),xmax);
	EXPECT_EQ(lattice2D.getYmax(),ymax);

	LatticeShape lattice3D(3);
	lattice3D.setXmax(xmax);
	lattice3D.setYmax(ymax);
	lattice3D.setZmax(zmax);
	EXPECT_EQ(lattice3D.getDim(),3);
	EXPECT_EQ(lattice3D.getXmax(),xmax);
	EXPECT_EQ(lattice3D.getYmax(),ymax);
	EXPECT_EQ(lattice3D.getZmax(),zmax);
}



TEST(GenerateNeighborsTest, DISABLED_OneDistance1D) {
	int dim = 1;
	LatticeShape lattice1D(dim);
	lattice1D.setXmax(1);

	Basis basis(0,0);// = new Basis(0, 0); //new Basis1D;
	Neighbors neighbors;

/********** when the |distance| is 1 ****************/
	int distance = 1;
	// no neighbors
	basis[0] = 0;
	basis[1] = 1;
	generateNeighbors( basis, distance, lattice1D, neighbors);
	EXPECT_EQ(neighbors.size(), 0);
	neighbors.clear();

	// only left neighbor of first particle
	distance = -1;
	basis[0] = 1;
	basis[1] = 2;
	lattice1D.setXmax(2);
	generateNeighbors( basis, distance, lattice1D, neighbors);
	EXPECT_EQ(neighbors.size(), 1);
	EXPECT_EQ(neighbors[0][0], 0);
	EXPECT_EQ(neighbors[0][1], 2);
	neighbors.clear();

	// only right neighbor of the second particle
	distance = 1;
	basis[0] = 1;
	basis[1] = 2;
	lattice1D.setXmax(4);
	generateNeighbors( basis, distance, lattice1D, neighbors);
	EXPECT_EQ(neighbors.size(), 1);
	EXPECT_EQ(neighbors[0][0], 1);
	EXPECT_EQ(neighbors[0][1], 3);
	neighbors.clear();

	//left neighbors of the first and second particles
	distance = -1;
	basis[0] = 1;
	basis[1] = 3;
	lattice1D.setXmax(4);
	generateNeighbors( basis, distance, lattice1D, neighbors);
	EXPECT_EQ(neighbors.size(), 2);
	EXPECT_EQ(neighbors[0][0], 0);
	EXPECT_EQ(neighbors[0][1], 3);
	EXPECT_EQ(neighbors[1][0], 1);
	EXPECT_EQ(neighbors[1][1], 2);
	neighbors.clear();

/********** when the |distance| is 2 ****************/
	distance = 2;
	basis[0] = 5;
	basis[1] = 6;
	lattice1D.setXmax(10);
	generateNeighbors( basis, distance, lattice1D, neighbors);
	EXPECT_EQ(neighbors.size(), 2);
	// first right
	EXPECT_EQ(neighbors[0][0], 6);
	EXPECT_EQ(neighbors[0][1], 7);
	// second right
	EXPECT_EQ(neighbors[1][0], 5);
	EXPECT_EQ(neighbors[1][1], 8);
	neighbors.clear();

	distance = -2;
	basis[0] = 5;
	basis[1] = 6;
	lattice1D.setXmax(10);
	generateNeighbors( basis, distance, lattice1D, neighbors);
	EXPECT_EQ(neighbors.size(), 2);
	// first left
	EXPECT_EQ(neighbors[0][0], 3);
	EXPECT_EQ(neighbors[0][1], 6);
	// second left
	EXPECT_EQ(neighbors[1][0], 4);
	EXPECT_EQ(neighbors[1][1], 5);
	neighbors.clear();
}



TEST(GenerateIndexMatrixTest, DISABLED_FourSites) {
	int xmax = 3;
	LatticeShape lattice1D(1);
	lattice1D.setXmax(xmax);
	extern std::vector< std::vector< Basis > > VtoG;
	extern std::vector<int> DimsOfV;
	extern IMatrix IndexMatrix;


	generateIndexMatrix(lattice1D);
	int i = IndexMatrix(0,1);
	EXPECT_EQ(i,0);
	i = IndexMatrix(0,2);
	EXPECT_EQ(i,0);
	i = IndexMatrix(0,3);
	EXPECT_EQ(i,0);
	i = IndexMatrix(1,2);
	EXPECT_EQ(i,1);
	i = IndexMatrix(1,3);
	EXPECT_EQ(i,0);
	i = IndexMatrix(2,3);
	EXPECT_EQ(i,0);

	EXPECT_EQ(DimsOfV[0],0);
	EXPECT_EQ(DimsOfV[1],1);
	EXPECT_EQ(DimsOfV[2],1);
	EXPECT_EQ(DimsOfV[3],2);
	EXPECT_EQ(DimsOfV[4],1);
	EXPECT_EQ(DimsOfV[5],1);

	Basis basis = VtoG[1][0];
	EXPECT_EQ(basis[0],0); // K = 1
	EXPECT_EQ(basis[1],1);

	basis = VtoG[2][0];
	EXPECT_EQ(basis[0],0); // K = 2
	EXPECT_EQ(basis[1],2);


	basis = VtoG[3][0];
	EXPECT_EQ(basis[0],0); // K = 3
	EXPECT_EQ(basis[1],3);

	basis = VtoG[3][1];
	EXPECT_EQ(basis[0],1);
	EXPECT_EQ(basis[1],2);

	basis = VtoG[4][0];
	EXPECT_EQ(basis[0],1); // K = 4
	EXPECT_EQ(basis[1],3);


	basis = VtoG[5][0];
	EXPECT_EQ(basis[0],2); // K = 5
	EXPECT_EQ(basis[1],3);
}

TEST(GenerateIndexMatrixTest, DISABLED_Sites121) {
	int xmax = 121;
	LatticeShape lattice1D(1);
	lattice1D.setXmax(xmax);
	extern std::vector< std::vector< Basis > > VtoG;
	extern std::vector<int> DimsOfV;
	extern IMatrix IndexMatrix;

	generateIndexMatrix(lattice1D);
	int i = DimsOfV[121];
	EXPECT_EQ(i,61);
}

TEST(GetLatticeIndex, DISABLED_TestPointerArgument) {
	int xmax = 99;
	LatticeShape lattice1D(1);
	lattice1D.setXmax(xmax);
	LatticeShape *p = &lattice1D;
	Basis basis(23,45);
	int site1, site2;
	getLatticeIndex(*p, basis, site1, site2);
	EXPECT_EQ(site1, 23);
	EXPECT_EQ(site2, 45);

	xmax = 4;
	int ymax = 4;
	LatticeShape lattice2D(2);
	lattice2D.setXmax(xmax);
	lattice2D.setYmax(ymax);
	p = &lattice2D;
	Basis basis2D(2,3,0,3);
	basis = basis2D;
	getLatticeIndex(*p, basis, site1, site2);
	EXPECT_EQ(site1, 17);
	EXPECT_EQ(site2, 15);
}
