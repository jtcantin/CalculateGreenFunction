/*
 * textIO.h
 *
 *  Created on: Aug 20, 2014
 *      Author: pxiang
 */

#ifndef TEXTIO_H_
#define TEXTIO_H_

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <complex>
#include <limits>
#include "../Utility/types.h"
#include "../Utility/misc.h"
#include "../Utility/random_generator.h"
#include "binaryIO.h"

typedef std::numeric_limits< double > dbl;



/**
 * save a complex matrix into a text file
 *
 * row_index  col_index   real_part   imag_part
 */
void saveMatrixText(std::string filename, CDMatrix& m, bool exactValue=false);


/**
 * save a real matrix into a text file
 *
 * row_index  col_index   matrix_element
 */
void saveMatrixText(std::string filename, DMatrix& m, bool exactValue=false);


/**
 * load a complex matrix from a text file (slow compared with the binary version)
 *
 * row_index  col_index   real_part   imag_part
 */
void loadMatrixText(std::string filename, CDMatrix& m);


/**
 * load a real matrix from a text file (slow compared with the binary version)
 *
 * row_index  col_index   matrix_element
 */
void loadMatrixText(std::string filename, DMatrix& m);

#endif /* TEXTIO_H_ */
