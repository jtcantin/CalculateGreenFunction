/*
 * types.h
 *
 *  Created on: Dec 17, 2013
 *      Author: pxiang
 */

#ifndef TYPES_H_
#define TYPES_H_

#include <math.h>
#include <cmath>
#define _USE_MATH_DEFINES

#include<iostream>

#include <complex>

//#include "mkl.h"

#define delta(x, y) ((x)==(y))?1:0




typedef std::complex<double> dcomplex;

// use eigen c++ library
#define EIGEN_USE_MKL_ALL
#include <Eigen/Dense>


typedef Eigen::MatrixXcd CDMatrix;

typedef Eigen::MatrixXd DMatrix;


typedef Eigen::MatrixXi IMatrix; // integer matrix

typedef Eigen::VectorXcd CDVector;


typedef Eigen::VectorXd DVector;

typedef Eigen::ArrayXd DArray;

typedef Eigen::ArrayXcd CDArray;



#include <list>
#include <vector>


typedef std::vector<double> Vector;










#endif /* TYPES_H_ */
