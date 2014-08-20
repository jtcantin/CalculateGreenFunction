/*
 * readInput.cpp
 *
 *  Created on: Aug 18, 2014
 *      Author: pxiang
 */

#include "readInput.h"

/**
 * convert from string to boolean
 */
bool to_bool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}





//std::string stringify(double x) {
////	std::ostringstream strs;
////	strs << x;
////	std::string str = strs.str();
//	return std::to_string(x);
//}
