/*
 * MatrixIO.h
 *
 *  Created on: Oct 31, 2014
 *      Author: pxiang
 */

#ifndef MATRIXIO_H_
#define MATRIXIO_H_

#include "binaryIO.h"
#include "textIO.h"

/**
 * save matrix to file
 *
 * if the filename is *.bin, save the matrix in binary format
 * otherwise save the matrix in text format
 */
void saveMatrix(const std::string filename, CDMatrix& gf);

/**
 * save matrix to file
 *
 * if the filename is *.bin, save the matrix in binary format
 * otherwise save the matrix in text format
 */
void saveMatrix(const std::string filename, DMatrix& gf);


/**
 * read a file and save its content into a matrix
 */
void loadMatrix(const std::string filename, CDMatrix& gf);

/**
 * read a file and save its content into a matrix
 */
void loadMatrix(const std::string filename, DMatrix& gf);


#endif /* MATRIXIO_H_ */
