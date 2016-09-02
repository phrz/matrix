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

TEST_CASE("one equals one", "[matrix]") {
	REQUIRE(1 == 1);
}
