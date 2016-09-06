/**
 * @file MatrixTest.cpp
 * @brief Tests for Matrix/Vector with the Catch testing library.
 *
 *
 * This file contains tests for the Matrix/Vector classes.
 * it is an adaptation of SMU Mathematics' Dr. Daniel R. Reynolds' tests in the
 * original "matrix_test.cpp" file. This file mainly reorganizes those tests
 * into units and handles API differences in this version of the code.
 *
 * Original repository: [drreynolds/matrix](https://bitbucket.org/drreynolds/matrix)
 *
 * @author Paul Herz
 * @date 2016-08-27
 */

#define CATCH_CONFIG_MAIN
#include "catch.h"

#include "Matrix.h"
#include "Vector.h"


TEST_CASE("Dimensions struct should initialize values properly", "[dimensions]") {
	auto d = PH::Dimensions(1,1);
	REQUIRE(d.rows == 1);
	REQUIRE(d.columns == 1);
}


TEST_CASE("Dimensions should be comparable", "[dimensions]") {
	auto d1 = PH::Dimensions(0,1);
	auto d2 = d1;
	auto d3 = PH::Dimensions(1,0);
	
	REQUIRE(d1 == d2);
	REQUIRE_FALSE(d1 != d2);
	
	REQUIRE(d1 != d3);
	REQUIRE_FALSE(d1 == d3);
}


TEST_CASE("Matrix default constructor should create a (0,0) matrix", "[matrix]") {
	auto m = PH::Matrix();
	
	REQUIRE(m.rows() == 0);
	REQUIRE(m.columns() == 0);
}

TEST_CASE("Matrix(r,c) constructor should fill a (r,c) matrix with zeroes", "[matrix]") {
	auto m = PH::Matrix(12, 10);
	
	REQUIRE(m(11, 9) == 0.0);
}

TEST_CASE("Matrix should be constructable by nested std::initializer_list objects", "[matrix]") {
	auto m = PH::Matrix({
		{1, 2, 3},
		{4, 5, 6},
		{7, 8, 9}
	});
	
	for(Index r = 0; r < m.rows(); ++r) {
		for(Index c = 0; c < m.columns(); ++c) {
			REQUIRE(m(r, c) == 1 + r * m.columns() + c);
		}
	}
	
	PH::Matrix m2 = {
		{1, 2, 3},
		{4, 5, 6},
		{7, 8, 9}
	};
	
	for(Index r = 0; r < m2.rows(); ++r) {
		for(Index c = 0; c < m2.columns(); ++c) {
			REQUIRE(m2(r, c) == 1 + r * m2.columns() + c);
		}
	}
	
	REQUIRE_THROWS(PH::Matrix m2({
		{1, 2, 3},
		{4, 5},
		{7}
	}));
}


TEST_CASE("fromArray should build a 2D Matrix from a 1D array", "[matrix]") {
	std::vector<double> data = {1,2,3,4,5,6,7,8,9};
	auto m = PH::Matrix::fromArray<3,3>(data);
	
	for(Index r = 0; r < m.rows(); ++r) {
		for(Index c = 0; c < m.columns(); ++c) {
			CAPTURE("\n" + m.str());
			CAPTURE(m(r, c));
			CAPTURE(r);
			CAPTURE(c);
			REQUIRE(m(r, c) == 1 + r * m.columns() + c);
		}
	}
	
	// empty data
	std::vector<double> emptyData = {};
	REQUIRE( (PH::Matrix::fromArray<0,0>(emptyData) == PH::Matrix()) );
	
	// too small
	std::vector<double> smallData = {1, 2};
	REQUIRE_THROWS( (PH::Matrix::fromArray<1, 3>(smallData)) );
	
	// too big
	std::vector<double> bigData = {1, 2, 3, 4};
	REQUIRE_THROWS( (PH::Matrix::fromArray<1, 3>(bigData)) );
}


TEST_CASE("Matrix should be constructable by 2D std::vector", "[matrix]") {
	Raw2DArray goodData = {{1,2,3},{4,5,6},{7,8,9}};
	Raw2DArray badData = {{1,2,3},{4}};
	
	// functional constructor
	auto m1 = PH::Matrix(goodData);
	
	for(Index r = 0; r < m1.rows(); ++r) {
		for(Index c = 0; c < m1.columns(); ++c) {
			REQUIRE(m1(r, c) == 1 + r * m1.columns() + c);
		}
	}
	
	// assignment constructor
	PH::Matrix m2 = goodData;
	
	for(Index r = 0; r < m2.rows(); ++r) {
		for(Index c = 0; c < m2.columns(); ++c) {
			REQUIRE(m2(r, c) == 1 + r * m2.columns() + c);
		}
	}
	
	REQUIRE_THROWS(auto a = PH::Matrix(badData));
}


TEST_CASE("eye(size) should create an identity matrix", "[matrix]") {
	// size zero
	REQUIRE(PH::Matrix::eye(0) == PH::Matrix());
	
	// some size
	auto eye = PH::Matrix::eye(10);
	REQUIRE(eye.isSquare());
	
	// test for definition of identity matrix
	for(Index r = 0; r < eye.rows(); ++r) {
		for(Index c = 0; c < eye.columns(); ++c) {
			if(r == c) {
				REQUIRE(eye(r,c) == 1);
			} else {
				REQUIRE(eye(r,c) == 0);
			}
		}
	}
}


TEST_CASE("random(r,c) should create a matrix of uniform distribution 0...1", "[matrix]") {
	// size zero
	REQUIRE(PH::Matrix::random(0,0) == PH::Matrix());
	
	// some size
	auto random = PH::Matrix::random(10, 10);
	
	// test for not all being same, 0 ≤ x ≤ 1
	auto allSame = true;
	auto prevValue = 0.0;
	random.mapElements([&](MathNumber& element, Index r, Index c) {
		// make sure not all same
		if(allSame && r != 0 || c != 0) {
			if(prevValue != element) allSame = false;
		}
		
		// check range
		REQUIRE(element >= 0.0);
		REQUIRE(element <= 1.0);
	});
	
	REQUIRE_FALSE(allSame);
}


TEST_CASE("linSpace(a,b,r,c) should create a linear span along a matrix", "[matrix]") {
	auto l1 = PH::Matrix::linSpace(1, 9, 3, 3);
	CAPTURE("\n"+l1.str());
	
	l1.mapElements([&](MathNumber& element, Index r, Index c) {
		auto linearIndex = 1 + r * l1.columns() + c;
		REQUIRE(element == linearIndex);
	});
	
	// upper bound below lower bound
	REQUIRE_THROWS(PH::Matrix::linSpace(9, 1, 3, 3));
}


TEST_CASE("logSpace(a,b,r,c) should create a logarithmic span along a matrix", "[matrix]") {
	auto l2 = PH::Matrix::logSpace(1, 9, 3, 3);
	CAPTURE("\n"+l2.str());
	
	l2.mapElements([&](MathNumber& element, Index r, Index c) {
		auto linearIndex = 1 + r * l2.columns() + c;
		REQUIRE(element == std::pow(10.0, linearIndex));
	});
	
	// upper bound below lower bound
	REQUIRE_THROWS(PH::Matrix::logSpace(9, 1, 3, 3));
}


