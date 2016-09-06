//
//  Vector.h
//  matrix
//
//  Created by Paul Herz on 8/28/16.
//  Copyright © 2016 Paul Herz. All rights reserved.
//

#ifndef Vector_h
#define Vector_h

#include "Types.h"
#include "Serializable.h"

namespace PH {
	
	class Vector: public Serializable<Vector> {
	protected:
		Raw1DArray _data;
	public:
		
		// empty constructor
		Vector();
		
		// size constructor, zero-fill
		Vector(Index n);
		
		// size and custom fill constructor
		Vector(Index n, MathNumber fill);
		
		Vector(Raw1DArray data);
		
		Vector(std::initializer_list<double> l);
		
		static Vector random(const Index size);
		
		// creates a vector of n linearly spaced values from a through b
		static Vector linSpace(double a, double b, Index n);
		
		// creates a vector of n logarithmically spaced values from 10^a through 10^b
		static Vector logSpace(double a, double b, Index n);
		
		// Serializable implementation
		void serialize(std::ostream& output) const;
		Vector deserialize(std::istream& input);
		
		void mapElements(std::function<void(MathNumber&,Index)> callback) {
			for(auto it = _data.begin(); it != _data.end(); ++it) {
				Index i = it - _data.begin();
				callback(*it, i);
			}
		}
		
		MathNumber& operator[](const Index i) {
			return _data[i];
		}
		
		MathNumber operator[](const Index i) const {
			return _data[i];
		}
		
		Index size() const {
			return _data.size();
		}
		
		// extract/insert routines for portions of vectors
		// y = x(is:ie)
		Vector range(Index begin, Index end);
		
		// v(begin:end) = source
		void insert(Index begin, Index end, Vector& source);
		
		// arithmetic operators for MathVector
		Vector& operator+=(const MathNumber addend);
		Vector& operator+=(const Vector& addend);
		
		Vector& operator-=(const MathNumber subtrahend);
		Vector& operator-=(const Vector& subtrahend);
		
		Vector& operator*=(const MathNumber c);
		Vector& operator*=(const Vector& w);
		
		Vector& operator/=(const MathNumber divisor);
		Vector& operator/=(const Vector& divisors);
		
		Vector& operator^=(const MathNumber exponent);
		Vector& operator^=(const Vector& exponents);
		
		MathNumber static dot(const Vector& v1, const Vector& v2);
		
		static MathNumber norm(const Vector& v);
		static MathNumber infNorm(const Vector& v);
		static MathNumber oneNorm(const Vector& v);
		
		std::string str(int precision = fullPrecision) const;
		
	}; // class Vector
	
	Vector operator+(const Vector& v, const MathNumber c);
	Vector operator+(const MathNumber c, const Vector& v);
	Vector operator+(const Vector& v, const Vector& w);
	
	Vector operator-(const Vector& v, const double c);
	Vector operator-(const MathNumber c, const Vector& v);
	Vector operator-(const Vector& v, const Vector& w);
	
	Vector operator*(const Vector& v, const double c);
	Vector operator*(const MathNumber c, const Vector& v);
	Vector operator*(const Vector& v, const Vector& w);
	Vector operator-(const Vector& v);
	
	Vector operator/(const Vector& v, const double c);
	Vector operator/(const MathNumber c, const Vector& v);
	Vector operator/(const Vector& v, const Vector& w);
	
	Vector operator^(const Vector& v, const double c);
	Vector operator^(const MathNumber c, const Vector& v);
	Vector operator^(const Vector& v, const Vector& w);
	
	// output routines
	std::ostream& operator<<(std::ostream& os, const Vector& v);
	
} // namespace PH

#endif /* Vector_h */
