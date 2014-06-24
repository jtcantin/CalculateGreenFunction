/*
 * mics.h
 *
 *  Created on: Dec 19, 2013
 *      Author: pxiang
 */

#ifndef MICS_H_
#define MICS_H_

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <complex>
#include "types.h"
#include "misc.h"

std::string itos(int n);
std::string stringConverter(int K, int nth);

// save two vectors of doubles to a file
void save_two_arrays(std::string filename, std::vector<double>& array1,
                     std::vector<double>& array2);

void save_two_arrays(std::string filename, std::vector<int>& array1,
                     std::vector<double>& array2);

// save two vectors of doubles to a file
void save_two_arrays(std::string filename, std::vector<double>& array1,
                     std::vector<dcomplex>& array2) ;

int max(int a, int b);

int min(int a, int b);

double stringToDouble(std::string numString);

int stringToInteger(std::string numString);

std::complex<double> stringToComplex(std::string numString);

std::string numToString(double number);

void read_double_column(std::string filename, int column, std::vector<double>& array);

std::vector<double> linspace(double start, double stop, int num=100);


#endif /* MICS_H_ */
