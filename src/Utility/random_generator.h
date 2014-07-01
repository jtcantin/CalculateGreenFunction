/*
 * random_generator.h
 *
 *  Created on: Dec 17, 2013
 *      Author: pxiang
 */

#ifndef RANDOM_GENERATOR_H_
#define RANDOM_GENERATOR_H_

#include "types.h"

class RandomNumberGenerator {
	private:
		unsigned seed;

	public:


// constructor
	RandomNumberGenerator(unsigned inputSeed) {
		seed = inputSeed;
		srand(seed);
	}

	// default constructor
	RandomNumberGenerator() {
		//use the current time as seed for random numbers
		seed = (unsigned)time(NULL);
		srand(seed);
	}

	void SetSeed(unsigned inputSeed) {
		seed = inputSeed;
		srand(seed);
	}

// produce a complex random number in [0, 1) + [0, 1)*I given the seed
	dcomplex randomComplex() {
		double randomReal = (double)rand()/(double)RAND_MAX;
		double randomImag = (double)rand()/(double)RAND_MAX;
		dcomplex rc = dcomplex(randomReal,randomImag);
		return rc;
	}

	// produce a complex random number in [0, 1) + [0, 1)*I given the seed
	double randomReal() {
		double random = (double)rand()/(double)RAND_MAX;
		return random;
	}


};

// generate a matrix of given dimension whose elements are random numbers in [0, 1)
class RandomNumberMatrix {
public:
	RandomNumberMatrix(int row_count, int col_count, unsigned seed) {
		rows = row_count;
		cols = col_count;
		p = NULL;
		p = new double[rows*cols];
		rng.SetSeed(seed);
		for (int i=0; i<rows*cols; i++) {
			p[i] = rng.randomReal();
		}
	}

	double operator()(const int r, const int c) const {
		return p[r*cols + c];
	}


	~RandomNumberMatrix() {
		delete [] p;
		p = NULL;
		rows = 0;
		cols = 0;
	}

private:
	double *p;
	int rows;
	int cols;
	RandomNumberGenerator rng;
};



#endif
