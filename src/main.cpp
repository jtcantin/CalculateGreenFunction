/*
 * main.cpp
 *
 *  Created on: Jun 11, 2014
 *      Author: pxiang
 */

#include "IO/binaryIO.h"
#include "Utility/misc.h"
#include "Utility/random_generator.h"
#include "Utility/types.h"
#include "gtest/gtest.h"

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	int i= RUN_ALL_TESTS();
	return 0;
}


