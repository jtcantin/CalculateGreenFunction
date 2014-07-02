/*
 * mics.cpp
 *
 *  Created on: Dec 19, 2013
 *      Author: pxiang
 */
#include "misc.h"


std::string itos(int n) {
	std::string result;
	std::ostringstream convert;   // stream used for the conversion
	convert << n;      // insert the textual representation of 'Number' in the characters in the stream
	result = convert.str(); // set 'Result' to the contents of the stream
	return result;
}

std::string stringConverter(int K, int nth) {
	return itos(K) + " " + itos(nth);
}

// save two vectors of doubles to a file
void save_two_arrays(std::string filename, std::vector<double>& array1,
                     std::vector<double>& array2) {
	std::ofstream out(filename.c_str());
  for (int i=0; i< array1.size(); ++i) {
		out << array1[i] <<"  " << array2[i] << std::endl;
  }
  out.close();
}


//
void save_two_arrays(std::string filename, std::vector<int>& array1,
                     std::vector<double>& array2) {
	std::ofstream out(filename.c_str());
  for (int i=0; i< array1.size(); ++i) {
		out << array1[i] <<"  " << array2[i] << std::endl;
  }
  out.close();
}


//
void save_two_arrays(std::string filename, std::vector<double>& array1,
                     std::vector<dcomplex>& array2) {
	std::ofstream out(filename.c_str());
  for (int i=0; i< array1.size(); ++i) {
		out << array1[i] <<"\t" << array2[i].real() << "\t" << array2[i].imag() << std::endl;
  }
  out.close();
}

int max(int a, int b) {
	return a>=b?a:b;
}

int min(int a, int b) {
	return a<=b?a:b;
}

// convert a string to double
double stringToDouble(std::string numString){
	return atoi(numString.c_str());
}

// convert a string to integer
int stringToInteger(std::string numString){
	return atof(numString.c_str());
}

// note a complex number is printed out like this:
// (0.0001703134203,0.0003920507577)
// convert a string to complex number
std::complex<double> stringToComplex(std::string numString){
	std::complex<double> complexNum;
	std::istringstream iss(numString);
	iss >> complexNum;
	return complexNum;
}


std::string numToString(double number) {
	std::stringstream ss (std::stringstream::in | std::stringstream::out);
  ss /* << std::scientific */<< number;
  std::string numString;
  ss >> numString;
  return numString;
}

// read a column of double numbers and save them into a vector
// NOTE: IT DOESN'T WORK
void read_double_column(std::string filename, int column, std::vector<double>& array) {
	std::ifstream ifs;
	ifs.open(filename.c_str());
	std::string line;
	std::string token;
	while (!ifs.eof()) {
		getline(ifs, line);
		if (line[0]!='#' && !line.empty()) { // ignore comments starting with "#"
			std::istringstream iss(line);
			for (int i=0; i<column; ++i){
				iss >> token;
			}
			array.push_back(stringToDouble(token));
		}
	}
	ifs.close();
}


// Return evenly spaced numbers over a specified interval
std::vector<double> linspace(double start, double stop, int num) {
	std::vector<double> array;
	double interval = (stop - start)/(num -1);
	double current = start;
	// add the start point
	array.push_back(current);
	for (int i=1; i<num-1; ++i) {
		current += interval;
		array.push_back(current);
	}
	// add the end point
	array.push_back(stop);
	return array;
}

/**
 * save a complex matrix into a text file
 *
 * row_index  col_index   real_part   imag_part
 */
void saveToFile(std::string filename, CDMatrix& m) {
	std::ofstream out(filename.c_str());
	for (int row=0; row< m.rows(); ++row) {
		for (int col=0; col<m.cols(); ++col) {
			dcomplex v = m(row, col);
			out << row <<"  " << col << "  "
			<< v.real() << "  "
			<< v.imag() << std::endl;
		}
	}
	out.close();
}


/**
 * save a real matrix into a text file
 *
 * row_index  col_index   matrix_element
 */
void saveToFile(std::string filename, DMatrix& m) {
	std::ofstream out(filename.c_str());
	for (int row=0; row< m.rows(); ++row) {
		for (int col=0; col<m.cols(); ++col) {
			double v = m(row, col);
			out << row <<"  " << col << "  "
			<< v << std::endl;
		}
	}
	out.close();
}
