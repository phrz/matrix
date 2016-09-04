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
			CAPTURE(r);
			CAPTURE(c);
			REQUIRE(m(r, c) == 1 + r * m.columns() + c);
		}
	}
	
	// TODO: test that size mismatch (too big, too small) throws
}

