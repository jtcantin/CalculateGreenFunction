/*
 * binaryIO.cpp
 *
 *  Created on: Jan 2, 2014
 *      Author: pxiang
 */
#include "binaryIO.h"



// save an eigen matrix into a binary file
void saveMatrix(std::string filename, const CDMatrix& m) {
	std::ofstream f(filename.c_str(), std::ios::binary);
	// write the row_count and col_count into the file
	int rows = m.rows();
	int cols = m.cols();
	f.write((char *)&rows, sizeof(rows));
	f.write((char *)&cols, sizeof(cols));
	// write the matrix elements into the file
	f.write((char *)m.data(), sizeof(dcomplex)*rows*cols);
	f.close();
}

// load an eigen matrix from a binary file
void loadMatrix(std::string filename, CDMatrix& m) {
	int rows, cols;
	std::ifstream f(filename.c_str(), std::ios::binary);
	f.read((char *)&rows, sizeof(rows));
	f.read((char *)&cols, sizeof(cols));
	m.resize(rows, cols);
	f.read((char *)m.data(), sizeof(dcomplex)*rows*cols);
	//if (f.bad())
	//	throw std::exception("Error reading matrix");
	f.close();
}


// save an eigen matrix into a binary file
void saveMatrix(std::string filename, const DMatrix& m) {
	std::ofstream f(filename.c_str(), std::ios::binary);
	// write the row_count and col_count into the file
	int rows = m.rows();
	int cols = m.cols();
	f.write((char *)&rows, sizeof(rows));
	f.write((char *)&cols, sizeof(cols));
	// write the matrix elements into the file
	f.write((char *)m.data(), sizeof(double)*rows*cols);
	f.close();
}

// load an eigen matrix from a binary file
void loadMatrix(std::string filename, DMatrix& m) {
	int rows, cols;
	std::ifstream f(filename.c_str(), std::ios::binary);
	f.read((char *)&rows, sizeof(rows));
	f.read((char *)&cols, sizeof(cols));
	m.resize(rows, cols);
	f.read((char *)m.data(), sizeof(double)*rows*cols);
	//if (f.bad())
	//	throw std::exception("Error reading matrix");
	f.close();
}
