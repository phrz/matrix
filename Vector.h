//
//  Vector.h
//  matrix
//
//  Created by Paul Herz on 8/28/16.
//  Copyright © 2016 Paul Herz. All rights reserved.
//

#ifndef Vector_h
#define Vector_h

#include "PHMath.h"
#include "Serializable.h"

namespace PH {
	
	class Vector: public Serializable<Vector> {
	protected:
		Raw1DArray _data;
	public:
		
		// Serializable implementation
		void serialize(std::ostream& output) const;
		Vector& deserialize(std::istream& input);
		
		// inner product between two vectors
		double dot(const MathVector& v1, const MathVector& v2);
		
		// norms
		double norm(const MathVector& v);      // sqrt(sum_i vi^2) (vector 2-norm)
		double infNorm(const MathVector& v);   // max_i |vi|       (vector inf-norm)
		double oneNorm(const MathVector& v);   // sum_i |vi|       (vector 1-norm)
		
		// creates a vector of n linearly spaced values from a through b
		MathVector linSpace(double a, double b, Index n);
		
		// creates a vector of n logarithmically spaced values from 10^a through 10^b
		MathVector logSpace(double a, double b, Index n);
		
		// creates a vector of n uniformly-distributed random values
		MathVector randomVectorOfSize(Index n);
		
		// output routines
		std::ostream& operator<<(std::ostream& os, const MathVector& v);
		
		// extract/insert routines for portions of vectors
		MathVector vecExtract(MathVector& x,       // y = x(is:ie)
							  long int is, long int ie);
		int vecInsert(MathVector& x, long int is,           // x(is:ie) = y
					  long int ie, MathVector& y);
		
		// arithmetic operators for MathVector
		MathVector& operator+=(MathVector& v, const double c);
		MathVector& operator+=(MathVector& v, const MathVector& w);
		
		MathVector& operator-=(MathVector& v, const double c);
		MathVector& operator-=(MathVector& v, const MathVector& w);
		
		MathVector& operator*=(MathVector& v, const double c);
		MathVector& operator*=(MathVector& v, const MathVector& w);
		
		MathVector& operator/=(MathVector& v, const double c);
		MathVector& operator/=(MathVector& v, const MathVector& w);
		
		MathVector& operator^=(MathVector& v, const double c);
		MathVector& operator^=(MathVector& v, const MathVector& w);
		
		MathVector operator+(const MathVector& v, const double c);
		MathVector operator+(const double c, const MathVector& v);
		MathVector operator+(const MathVector& v, const MathVector& w);
		MathVector operator-(const MathVector& v, const double c);
		MathVector operator-(const double c, const MathVector& v);
		MathVector operator-(const MathVector& v, const MathVector& w);
		MathVector operator*(const MathVector& v, const double c);
		MathVector operator*(const double c, const MathVector& v);
		MathVector operator*(const MathVector& v, const MathVector& w);
		MathVector operator/(const MathVector& v, const double c);
		MathVector operator/(const double c, const MathVector& v);
		MathVector operator/(const MathVector& v, const MathVector& w);
		MathVector operator^(const MathVector& v, const double c);
		MathVector operator^(const double c, const MathVector& v);
		MathVector operator^(const MathVector& v, const MathVector& w);
		
	}; // class Vector
	
	//--- supplementary matrix-vector arithmetic routines ---
	
	// standard matrix-vector product -> new vector (function form)
	Optional<Vector> matrixVectorProduct(const Matrix& A, const MathVector& v);
	
	// standard matrix-vector product -> new vector
	Vector operator*(const Matrix& A, const MathVector& v);
	
} // namespace PH

#endif /* Vector_h */
