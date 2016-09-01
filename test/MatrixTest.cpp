/**
 * @file MatrixTest.cpp
 * @brief Tests for the Matrix and Vector classes with Google Test.
 *
 *
 * This file contains tests for the Matrix/Vector classes. It is an adaptation of
 * SMU Mathematics' Dr. Daniel R. Reynolds' tests in the
 * original "matrix_test.cpp" file. This file mainly reorganizes those tests
 * into units and handles API differences in this version of the code.
 *
 * Original repository: [drreynolds/matrix](https://bitbucket.org/drreynolds/matrix)
 *
 * @author Paul Herz
 * @date 2016-08-27
 */

#include <iostream>

#include "gtest/gtest.h"
#include "Matrix.h"
#include "Vector.h"

class MatrixTest: public ::testing::Test {};

class VectorTest: public ::testing::Test {};

//TEST_F(MatrixTest, Constructor) {
//	EXPECT_EQ(1, 1);
//}

TEST_F(VectorTest, DefaultConstructorCreatesEmptyVector) {
	auto v = PH::Vector();
	
}

int main(int argc, char** argv) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


