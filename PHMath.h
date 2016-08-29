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
#include "optional.h"

using  MathNumber = double;
using  Raw1DArray = std::vector<MathNumber>;
using  Raw2DArray = std::vector<Raw1DArray>;
using       Index = std::size_t;

template<class T>
using Optional = std::experimental::optional<T>;

#include "Matrix.h"
#include "Vector.h"

namespace PH {
	
	//--- supplementary matrix arithmetic routines ---
	
	double Dot(const Matrix& A, const Matrix& B);    // sum_ij (Cij * Aij)
	double Norm(const Matrix& A);                    // sqrt(sum_ij Cij^2)
	double InfNorm(const Matrix& A);                 // max_i sum_j |Cij|
	double OneNorm(const Matrix& A);                 // max_j sum_i |Cij|
	
	//--- new matrix creation routines (C is the output, A and B are the operands) ---
	
	Matrix operator+(const Matrix& A, const Matrix &B);        // C = A+B
	Matrix operator-(const Matrix& A, const Matrix &B);        // C = A-B
	Matrix operator*(const Matrix& A, const Matrix &B);        // C = A*B
	Matrix operator*(const Matrix& A, const double b);         // C = A*b
	Matrix operator*(const double a, const Matrix& B);         // C = a*B
	Matrix linSpace(double a, double b, Index m, Index n);   // linear span
	Matrix logSpace(double a, double b, Index m, Index n);   // logarithmic span
	Matrix randomMatrixOfSize(Index m, Index n);                         // Cij random in [0,1]
	Matrix Eye(Index n);                                      // Cij = delta_i(j)
	Matrix MatrixRead(const char *infile);                     // creates from input file
	Matrix Inverse(const Matrix& A);                           // C = A^{-1}
	
	
	//--- linear algebra routines ---
	
	// backward substitution on the linear system U*X = B, returning X as a new Matrix
	//    U and B remain unchanged in this operation
	Optional<Matrix> backSubstitution(const Matrix& U, const Matrix& B);
	
	// backward substitution on U*x = b, returning x as a new Vector
	//    U and b remain unchanged in this operation
	Optional<Vector> backSubstitution(const Matrix& U, const Vector& b);
	
	// forward substitution on the linear system L*X = B, returning X as a new Matrix
	// L and B remain unchanged in this operation
	Optional<Matrix> forwardSubstitution(const Matrix& L, const Matrix& B);
	
	// forward substitution on L*x = b, returning x as a new vector<double>
	// L and b remain unchanged in this operation
	Optional<Vector> ForwardSubstitution(const Matrix& L, const Vector& b);
	
	// solves a linear system A*X = B, returning X as a new Matrix
	// A and B are modified in this operation; X holds the result
	Optional<Matrix> linearSolve(Matrix& A, Matrix& B);
	
	// solves a linear system A*x = b, filling in the input vector<double> x
	// A and b are modified in this operation; x holds the result
	Optional<Vector> linearSolve(Matrix& A, Vector& b);
	
}

#endif /* PHMath_h */
