/*
 * binaryIO.h
 *
 *  Created on: Jan 2, 2014
 *      Author: pxiang
 */

#ifndef BINARYIO_H_
#define BINARYIO_H_

#include <iostream>
#include <fstream>
#include <string>
#include "../Utility/random_generator.h"
#include "../Utility/types.h"
#include <complex>


void saveMatrixBin(std::string filename, const CDMatrix& m);
void loadMatrixBin(std::string filename, CDMatrix& m);

void saveMatrixBin(std::string filename, const DMatrix& m);
void loadMatrixBin(std::string filename, DMatrix& m);

#endif /* BINARYIO_H_ */
