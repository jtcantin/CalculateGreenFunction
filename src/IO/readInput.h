/*
 * readInput.h
 *
 *  Created on: Aug 18, 2014
 *      Author: pxiang
 */

#ifndef READINPUT_H_
#define READINPUT_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cctype>

#include <stdexcept>

class BadConversion : public std::runtime_error {
public:
  BadConversion(std::string const& s)
    : std::runtime_error(s)
    { }
};

// convert a variable into a string
#define xstr(s) str(s)
#define str(s) #s

/**
 * convert from string to boolean
 */
bool to_bool(std::string str);

class InputVariable {
public:
	// default constructor
	InputVariable(std::string type, std::string name, void * pValue) {
		type_ = type;
		name_ = name;
		pValue_ = pValue;
	}

	// go over each line of the file and find the input for the current variable
	void read(std::string filename) {
		std::ifstream iFile(filename.c_str());
		std::string line;

		/* While there is still a line. */
		while(getline(iFile, line)) {
			std::istringstream iss(line);
			std::string typeString;
			std::string nameString;
			std::string valueString;

			iss >> typeString >> nameString >> valueString;
			if (typeString==type_ && nameString==name_) {
				// convert valueString to the right type
				if (typeString=="int") {
					int * p = (int *) pValue_;
					*p = atoi(valueString.c_str());
				} else if (typeString=="double") {
					double * p = (double *) pValue_;
					*p = atof(valueString.c_str());
				} else if (typeString=="bool") {
					bool * p = (bool *) pValue_;
					*p = to_bool(valueString);
				} else if (typeString=="unsigned") {
					unsigned * p = (unsigned *) pValue_;
					*p = (unsigned) atoi(valueString.c_str());
				}
				break; // found the input value and jump out of the while loop
			}
		}

		iFile.close();
	}

private:
	std::string type_;
	std::string name_;
	void * pValue_;

};

#endif /* READINPUT_H_ */
