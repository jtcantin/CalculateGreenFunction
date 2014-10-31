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
			false,false,10,88,true,true};
	InputVariable inputVar("unsigned", xstr(interactionData.seed),
			                &(interactionData.seed));
	std::ofstream myfile;
	myfile.open ("seed.txt");
	myfile << "some random string"<<std::endl;
	myfile << "unsigned interactionData.seed 100"<<std::endl;
	myfile << "more random string"<<std::endl;
	myfile.close();
	inputVar.read("seed.txt");
	system("rm seed.txt");
	EXPECT_TRUE(interactionData.seed == 100);
}


TEST(ReadInputTest, AllParameters) {
	std::string inputFile = "tmp.txt";
	// create the input file
	std::ofstream myfile;
	myfile.open (inputFile.c_str());
	myfile << "some random string: qqn;^(7qytwlyh"<<std::endl;
	myfile << "unsigned seed 100"<<std::endl;
	myfile << "more random string: 2g;wyqp&926y"<<std::endl;
	myfile << "int xmax 101"<<std::endl;
	myfile << "double onsiteE 1.0"<<std::endl;
	myfile << "double hop 2.0"<<std::endl;
	myfile << "double dyn 3.0"<<std::endl;
	myfile << "more random string: agn;ahy;$264r3"<<std::endl;
	myfile << "bool randomOnsite true"<<std::endl;
	myfile << "bool randomHop false"<<std::endl;
	myfile << "bool randomDyn true"<<std::endl;
	myfile << "int initialSeparation 5"<<std::endl;
	myfile << "bool longRangeHop true"<<std::endl;
	myfile << "bool longRangeDyn true"<<std::endl;
	myfile << "int maxDistance 10"<<std::endl;
	myfile << "double E 0.0" << std::endl;
	myfile << "double eta 0.01" << std::endl;
	myfile.close();

	int xmax;
	int initialSeparation;
	double onsiteE, hop, dyn;
	bool randomOnsite, randomHop, randomDyn;
	int maxDistance;
	unsigned seed;
	bool longRangeHop, longRangeDyn;
	double E, eta;

	// read from the input file
	READ_INPUT(inputFile, "int", xmax);
	READ_INPUT(inputFile, "double", onsiteE);
	READ_INPUT(inputFile, "double", hop);
	READ_INPUT(inputFile, "double", dyn);
	READ_INPUT(inputFile, "bool", randomOnsite);
	READ_INPUT(inputFile, "bool", randomHop);
	READ_INPUT(inputFile, "bool", randomDyn);
	READ_INPUT(inputFile, "int", maxDistance);
	READ_INPUT(inputFile, "unsigned", seed);
	READ_INPUT(inputFile, "bool", longRangeHop);
	READ_INPUT(inputFile, "bool", longRangeDyn);
	READ_INPUT(inputFile, "int", initialSeparation);
	READ_INPUT(inputFile, "double", E);
	READ_INPUT(inputFile, "double", eta);

	InteractionData interactionData = { onsiteE, hop, dyn,
			                            randomOnsite, randomHop, randomDyn,
			                            maxDistance, seed,
	                                    longRangeHop, longRangeDyn};

	EXPECT_EQ(xmax, 101);
	EXPECT_EQ(interactionData.seed, 100);
	EXPECT_DOUBLE_EQ(interactionData.onsiteE, 1.0);
	EXPECT_DOUBLE_EQ(interactionData.hop, 2.0);
	EXPECT_DOUBLE_EQ(interactionData.dyn, 3.0);
	EXPECT_EQ(interactionData.randomOnSite, true);
	EXPECT_EQ(interactionData.randomHop, false);
	EXPECT_EQ(interactionData.randomDyn, true);
	EXPECT_EQ(interactionData.maxDistance, 10);
	EXPECT_EQ(interactionData.longRangeHop, true);
	EXPECT_EQ(interactionData.longRangeDyn, true);
	EXPECT_EQ(initialSeparation, 5);
	EXPECT_DOUBLE_EQ(E, 0.0);
	EXPECT_DOUBLE_EQ(eta, 0.01);
}
