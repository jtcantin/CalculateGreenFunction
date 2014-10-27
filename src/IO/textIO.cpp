/*
 * textIO.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: pxiang
 */


#include "textIO.h"

/**
 * save a complex matrix into a text file
 *
 * row_index  col_index   real_part   imag_part
 */
void saveMatrixText(std::string filename, CDMatrix& m, bool exactValue) {
	std::ofstream out(filename.c_str());
	// print to full precision
	if (exactValue) {
		out.precision(dbl::digits10 + 2);
	}
	// use default precision to save space
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
void saveMatrixText(std::string filename, DMatrix& m, bool exactValue) {
	std::ofstream out(filename.c_str());
	// print to full precision
	if (exactValue) {
		out.precision(dbl::digits10 + 2);
	}
	// use default precision to save space
	for (int row=0; row< m.rows(); ++row) {
		for (int col=0; col<m.cols(); ++col) {
			double v = m(row, col);
			out << row <<"  " << col << "  "
			<< v << std::endl;
		}
	}
	out.close();
}



/**
 * load a complex matrix from a text file (slow compared with the binary version)
 *
 * row_index  col_index   real_part   imag_part
 */
void loadMatrixText(std::string filename, CDMatrix& m) {
	std::ifstream ifs(filename.c_str());

	std::vector<double> row_indices;
	std::vector<double> col_indices;
	std::vector< dcomplex > values;

	std::string line;
	int row_count = 0;
	int col_count = 0;

	while(std::getline(ifs, line)) { // read one line from ifs
	    std::istringstream iss(line); // access line as a stream

	    int row_i, col_i;
	    double real_part, imag_part;

	    iss >> row_i >> col_i >> real_part >> imag_part;

	    if (row_i > row_count) {
	    	row_count = row_i;
	    }

	    if (col_i > col_count) {
	    	col_count = col_i;
	    }

	    row_indices.push_back(row_i);
	    col_indices.push_back(col_i);

	    dcomplex value(real_part, imag_part);
	    values.push_back(value);

	}

	ifs.close();

	row_count += 1;
	col_count += 1;

	m.resize(row_count, col_count);
	for (int i=0; i<row_indices.size(); ++i) {
		int row = row_indices[i];
		int col = col_indices[i];
		m(row, col) = values[i];
	}
}


/**
 * load a real matrix from a text file (slow compared with the binary version)
 *
 * row_index  col_index   matrix_element
 */
void loadMatrixText(std::string filename, DMatrix& m) {
	std::ifstream ifs(filename.c_str());

	std::vector<double> row_indices;
	std::vector<double> col_indices;
	std::vector<double> values;

	std::string line;
	int row_count = 0;
	int col_count = 0;

	while(std::getline(ifs, line)) { // read one line from ifs
	    std::istringstream iss(line); // access line as a stream

	    int row_i, col_i;
	    double value;

	    iss >> row_i >> col_i >> value;

	    if (row_i > row_count) {
	    	row_count = row_i;
	    }

	    if (col_i > col_count) {
	    	col_count = col_i;
	    }

	    row_indices.push_back(row_i);
	    col_indices.push_back(col_i);
	    values.push_back(value);

	}

	ifs.close();

	row_count += 1;
	col_count += 1;

	m.resize(row_count, col_count);
	for (int i=0; i<row_indices.size(); ++i) {
		int row = row_indices[i];
		int col = col_indices[i];
		m(row, col) = values[i];
	}
}
