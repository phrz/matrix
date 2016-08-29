//
//  PHMath.h
//  matrix
//
//  Created by Paul Herz on 8/28/16.
//  Copyright © 2016 Paul Herz. All rights reserved.
//

#ifndef PHMath_h
#define PHMath_h

// for accessing pi via M_PI
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <cmath>
#else
#include <cmath>
#endif

#include <vector>

using  MathNumber = double;
using  Raw1DArray = std::vector<MathNumber>;
using  Raw2DArray = std::vector<Raw1DArray>;
using       Index = std::size_t;

#include "Matrix.h"
#include "Vector.h"

namespace PH {
	
	// sqrt(sum_ij Cij^2)
	double Norm(const Matrix& A);
	// max_i sum_j |Cij|
	double InfNorm(const Matrix& A);
	
	// max_j sum_i |Cij|
	double OneNorm(const Matrix& A);
	
	// C = A+B
	Matrix operator+(const Matrix& A, const Matrix &B);
	
	// C = A-B
	Matrix operator-(const Matrix& A, const Matrix &B);
	
	// C = A*B
	Matrix operator*(const Matrix& A, const Matrix &B);
	
	// C = A*b
	Matrix operator*(const Matrix& A, const double b);
	
	// C = a*B
	Matrix operator*(const double a, const Matrix& B);
	
	// linear span
	Matrix linSpace(double a, double b, Index m, Index n);
	
	 // logarithmic span
	Matrix logSpace(double a, double b, Index m, Index n);
	
	// Cij random in [0,1]
	Matrix randomMatrixOfSize(Index m, Index n);
	
	// Cij = delta_i(j)
	Matrix Eye(Index n);
	
	// creates from input file
	Matrix MatrixRead(const char *infile);
	
	// C = A^{-1}
	Matrix Inverse(const Matrix& A);
	
	// sum_ij (Cij * Aij)
	MathNumber dotProduct(const Matrix& A, const Matrix& B);
	
	// inner product between two vectors
	MathNumber dotProduct(const Vector& v1, const Vector& v2);
	
	/// sqrt(sum_i vi^2)    (vector 2-norm)
	MathNumber norm(const Vector& v);
	
	/// max_i |vi|   (vector inf-norm)
	MathNumber infNorm(const Vector& v);
	
	/// sum_i |vi|   (vector 1-norm)
	MathNumber oneNorm(const Vector& v);
	
	// creates a vector of n linearly spaced values from a through b
	Vector linSpace(double a, double b, Index n);
	
	// creates a vector of n logarithmically spaced values from 10^a through 10^b
	Vector logSpace(double a, double b, Index n);
	
	// creates a vector of n uniformly-distributed random values
	Vector randomVectorOfSize(Index n);
	
	
	//--- linear algebra routines ---
	
	// backward substitution on the linear system U*X = B, returning X as a new Matrix
	//    U and B remain unchanged in this operation
	Matrix backSubstitution(const Matrix& U, const Matrix& B);
	
	// backward substitution on U*x = b, returning x as a new Vector
	//    U and b remain unchanged in this operation
	Vector backSubstitution(const Matrix& U, const Vector& b);
	
	// forward substitution on the linear system L*X = B, returning X as a new Matrix
	// L and B remain unchanged in this operation
	Matrix forwardSubstitution(const Matrix& L, const Matrix& B);
	
	// forward substitution on L*x = b, returning x as a new vector<double>
	// L and b remain unchanged in this operation
	Vector ForwardSubstitution(const Matrix& L, const Vector& b);
	
	// solves a linear system A*X = B, returning X as a new Matrix
	// A and B are modified in this operation; X holds the result
	Matrix linearSolve(Matrix& A, Matrix& B);
	
	// solves a linear system A*x = b, filling in the input vector<double> x
	// A and b are modified in this operation; x holds the result
	Vector linearSolve(Matrix& A, Vector& b);
	
	//--- supplementary matrix-vector arithmetic routines ---
	
	// standard matrix-vector product -> new vector (function form)
	Vector matrixVectorProduct(const Matrix& A, const Vector& v);
	
	// standard matrix-vector product -> new vector
	Vector operator*(const Matrix& A, const Vector& v);
	
}

#endif /* PHMath_h */
