/*
 * readInput_test.cpp
 *
 *  Created on: Aug 18, 2014
 *      Author: pxiang
 */

#include "gtest/gtest.h"
#include "readInput.h"
#include "../recursiveCalculation/recursiveCalculation.h"

TEST(StringificationTest, IntDouble) {
	int var1;
	double var2;
	std::string var1string = xstr(var1);
	std::string var2string = xstr(var2);
	EXPECT_TRUE(var1string=="var1");
	EXPECT_TRUE(var2string=="var2");
}


TEST(StringificationTest, Struct) {
	struct {
		int var1;
		double var2;
	} aStruct;

	std::string var1string = xstr(aStruct.var1);
	std::string var2string = xstr(aStruct.var2);
	EXPECT_TRUE(var1string=="aStruct.var1");
	EXPECT_TRUE(var2string=="aStruct.var2");
}


TEST(InputVariable, Unsigned) {
	InteractionData interactionData = {1.0,1.0,1.0,true,
			false,false,10,88};
	InputVariable inputVar("unsigned", xstr(interactionData.seed),
			                &(interactionData.seed));
	std::ofstream myfile;
	myfile.open ("seed.txt");
	myfile << "some random string"<<std::endl;
	myfile << "unsigned interactionData.seed 100"<<std::endl;
	myfile << "more random string"<<std::endl;
	myfile.close();
	inputVar.read("seed.txt");
	EXPECT_TRUE(interactionData.seed == 100);
}
