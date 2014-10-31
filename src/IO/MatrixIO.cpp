/*
 * MatrixIO.cpp
 *
 *  Created on: Oct 31, 2014
 *      Author: pxiang
 */

#include "MatrixIO.h"

void saveMatrix(const std::string filename, CDMatrix& gf) {
	if (filename.substr(filename.length()-3, 3)=="bin") {
		saveMatrixBin(filename, gf); // save as binary
	} else {
		// save the gf matrix into text file
		saveMatrixText(filename, gf);
	}
}


void saveMatrix(const std::string filename, DMatrix& gf) {
	if (filename.substr(filename.length()-3, 3)=="bin") {
		saveMatrixBin(filename, gf); // save as binary
	} else {
		// save the gf matrix into text file
		saveMatrixText(filename, gf);
	}
}


void loadMatrix(const std::string filename, CDMatrix& gf) {
	if (filename.substr(filename.length()-3, 3)=="bin") {
		loadMatrixBin(filename, gf); // save as binary
	} else {
		// save the gf matrix into text file
		loadMatrixText(filename, gf);
	}
}


void loadMatrix(const std::string filename, DMatrix& gf) {
	if (filename.substr(filename.length()-3, 3)=="bin") {
		loadMatrixBin(filename, gf); // save as binary
	} else {
		// save the gf matrix into text file
		loadMatrixText(filename, gf);
	}
}





