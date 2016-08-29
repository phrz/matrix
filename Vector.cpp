//
//  Vector.cpp
//  matrix
//
//  Created by Paul Herz on 8/28/16.
//  Copyright © 2016 Paul Herz. All rights reserved.
//

#include "Vector.h"
#include "NotImplementedException.h"

namespace PH {
	
	Vector::Vector(Index n) {
		_data = Raw1DArray(n);
	}
	
	Vector::Vector(Index n, MathNumber fill) {
		_data = Raw1DArray(n,fill);
	}
	
	// Serializable interface
	void Vector::serialize(std::ostream& output) const {
		
		// [CITE] http://stackoverflow.com/a/6693088/3592716
		// quick way to concatenate STL container with commas
		std::ostringstream ss;
		std::copy(_data.begin(), _data.end() - 1, std::ostream_iterator<int>(ss, ", "));
	}
	
	
	static Vector& Vector::deserialize(std::istream& input) {
		throw new NotImplementedException();
	}
	
	
	// extract routine for portions of vectors, output = this(begin:end)
	Vector Vector::range(Index begin, Index end) {
		
		// update start and end if any are negative
		start = (start < 0) ? start+this->size() : start;
		end = (start < 0) ? end+this->size() : end;
		
		if (source.size() != (end - begin + 1)) {
			// Checking that source vector matches provided subvector range
			throw new std::invalid_argument("insert: size mismatch, supplied vector has " + std::to_string(source.size()) + " entries, but requested subvector has " + std::to_string(end - begin + 1) + " entries.");
		} else if (begin < 0 || end < 0 || end >= this.size()) {
			// Checking that range indices are legal
			// (counting numbers, with end in bounds)
			throw new std::invalid_argument("insert: requested subvector does not exist: (" + std::to_string(begin) + ":" + std::to_string(end) + ") on a vector with " + std::to_string(this->size()) + " elements.");
		} else if (end < begin) {
			// Checking that lower index is below upper index
			// (no backwards ranges)
			throw new std::invalid_argument("insert: requested subvector does not exist, upper index is below lower index: (" + std::to_string(begin) + ":" + std::to_string(end) + ")");
		}
		
		// create new vector of desired size
		MathVector result(end - start + 1);
		
		// copy requested data
		for (Index i = start; i <= end; i++) {
			result[i-start] = _data[i];
		}
		
		return result;
	}
	
	
	// insert routine for portions of vectors, this(begin:end) = source
	void Vector::insert(Index begin,
				  Index end, Vector& source) {
		
		// update start and end if any are negative
		begin = (begin < 0) ? begin + this->size() : begin;
		end = (end < 0) ? end + this->size() : end;
		
		if (source.size() != (end - begin + 1)) {
			// Checking that source vector matches provided subvector range
			throw new std::invalid_argument("insert: size mismatch, supplied vector has " + std::to_string(source.size()) + " entries, but requested subvector has " + std::to_string(end - begin + 1) + " entries.");
		} else if (begin < 0 || end < 0 || end >= this.size()) {
			// Checking that range indices are legal
			// (counting numbers, with end in bounds)
			throw new std::invalid_argument("insert: requested subvector does not exist: (" + std::to_string(begin) + ":" + std::to_string(end) + ") on a vector with " + std::to_string(this->size()) + " elements.");
		} else if (end < begin) {
			// Checking that lower index is below upper index
			// (no backwards ranges)
			throw new std::invalid_argument("insert: requested subvector does not exist, upper index is below lower index: (" + std::to_string(begin) + ":" + std::to_string(end) + ")");
		}
		
		// perform operation
		for (Index i = 0; i < source.size(); i++) {
			_data[i + start] = source[i];
		}
	}
	
	
	// In place arithmetic operators (members)
	#pragma mark In-Place Member Operators
	
	/// ADDITION
	
	/// In-place, constant (vector-by-number) addition.
	MathVector& Vector::operator+=(const MathNumber addend) {
		this->mapElements([addend](MathNumber& element, Index i) {
			element += addend;
		});
		
		return *this;
	}
	
	/// In-place, elementwise (vector-by-vector) addition
	Vector& Vector::operator+=(const Vector& addends) {
		if (this->size() != addends.size()) {
			throw new std::invalid_argument("vector-by-vector addition: incompatible vector sizes.");
		}
		
		this->mapElements([addends](MathNumber& element, Index i) {
			element += addends[i];
		});
		
		return *this;
	}
	
	/// SUBTRACTION
	
	/// In-place, constant (vector-by-number) subtraction
	Vector& Vector::operator-=(const MathNumber subtrahend) {
		this->mapElements([subtrahend](MathNumber& element, Index i) {
			element -= subtrahend;
		});
		
		return *this;
	}
	
	/// In-place, elementwise (vector-by-vector) subtraction
	Vector& Vector::operator-=(const Vector& subtrahends) {
		if (this->size() != subtrahends.size()) {
			throw new std::invalid_argument("vector-by-vector subtraction: incompatible vector sizes.");
		}
		
		this->mapElements([subtrahends](MathNumber& element, Index i) {
			element -= subtrahends[i];
		});
		
		return *this;
	}
	
	/// MULTIPLICATION
	
	/// In-place, constant (vector-by-number) multiplication
	Vector& Vector::operator*=(const MathNumber multiplier) {
		this->mapElements([multiplier](MathNumber& element, Index i) {
			element *= multiplier;
		});
		
		return *this;
	}
	
	/// In-place, elementwise (vector-by-vector) multiplication
	Vector& Vector::operator*=(const MathVector& multipliers) {
		if (this->size() != multipliers.size()) {
			throw new std::invalid_argument("vector-by-vector multiplication: incompatible vector sizes.");
		}

		this->mapElements([multipliers](MathNumber& element, Index i) {
			element *= multipliers[i];
		});
		
		return *this;
	}
	
	/// DIVISION
	
	/// In-place, constant (vector-by-number) division
	Vector& Vector::operator/=(const MathNumber divisor) {
		this->mapElements([divisor](MathNumber& element, Index i) {
			element /= divisor;
		});
		
		return *this;
	}
	
	/// In-place, elementwise (vector-by-vector) division
	Vector& Vector::operator/=(const Vector& divisors) {
		if (this->size() != divisors.size()) {
			throw new std::invalid_argument("vector-by-vector division: incompatible vector sizes.");
		}
		
		this->mapElements([divisors](MathNumber& element, Index i) {
			element /= divisors[i];
		});
		
		return *this;
	}
	
	/// EXPONENTIATION
	
	/// In-place, constant (vector-by-number) exponentiation
	Vector& Vector::operator^=(const MathNumber exponent) {
		this->mapElements([exponent](MathNumber& element, Index i) {
			element = std::pow(element, exponent);
		});
		
		return *this;
	}
	
	/// In-place, elementwise (vector-by-vector) exponentiation
	MathVector& Vector::operator^=(const Vector& exponents) {
		if (this->size() != exponents.size()) {
			throw new std::invalid_argument("vector-by-vector exponentiation: incompatible vector sizes.");
		}
		
		this->mapElements([exponents](MathNumber& element, Index i){
			element = std::pow(element, exponents[i]);
		});
		
		return *this;
	}
	
	#pragma mark Functional Operators
	
	// Result = Vector + Constant
	Vector operator+(const Vector& vector, const MathNumber constant) {
		Vector result = vector;
		result += constant;
		return result;
	}
	
	// Result = Constant + Vector
	Vector operator+(const MathNumber constant, const Vector& vector) {
		return vector + constant; // commutative
	}
	
	// Result = Vector + Vector
	Vector operator+(const Vector& vector1, const Vector& vector2) {
		Vector result = vector1;
		result += vector2;
		return result;
	}
	
	// Result = Vector - Constant
	Vector operator-(const Vector& vector, const MathNumber constant) {
		Vector result = vector;
		result -= constant;
		return result;
	}
	
	// Result = Constant - Vector
	Vector operator-(const MathNumber constant, const Vector& vector) {
		Vector result = vector;
		result.mapElements([constant](MathNumber& element, Index i){
			element = constant - element;
		});
		return result;
	}
	
	// Result = Vector - Vector
	Vector operator-(const Vector& vector1, const Vector& vector2) {
		Vector result = vector1;
		result -= vector2;
		return result;
	}
	
	// Result = Vector * Constant
	Vector operator*(const Vector& vector, const MathNumber constant) {
		Vector result = vector;
		result *= constant;
		return result;
	}
	
	// Result = Constant * Vector
	Vector operator*(const MathNumber constant, const Vector& vector) {
		return vector * constant; // commutative
	}
	
	// Result = Vector * Vector
	Vector operator*(const Vector& vector1, const Vector& vector2) {
		Vector result = vector1;
		vector1 *= vector2;
		return result;
	}
	
	// Result = Vector / Constant
	Vector operator/(const Vector& vector, const MathNumber constant) {
		Vector result = vector1;
		result /= constant;
		return result;
	}
	
	// Result = Constant / Vector
	Vector operator/(const MathNumber constant, const Vector& vector) {
		Vector result = vector;
		result.mapElements([constant](MathNumber& element, Index i){
			element = constant / element;
		});
		return result;
	}
	
	// Result = Vector / Vector
	Vector operator/(const Vector& vector1, const Vector& vector2) {
		MathVector result = vector1;
		result /= vector2;
		return result;
	}
	
	// Result = Vector ^ Constant
	Vector operator^(const Vector& vector, const MathNumber constant) {
		Vector result = vector
		result ^= constant;
		return result;
	}
	
	// Result = Constant ^ Vector
	Vector operator^(const MathNumber constant, const MathVector& vector) {
		Vector result = vector;
		result.mapElements([constant](MathNumber& element, Index i){
			element = std::pow(constant, element);
		});
		return result;
	}
	
	// Result = Vector ^ Vector
	Vector operator^(const Vector& vector1, const Vector& vector2) {
		Vector result = vector1;
		result ^= vector2;
		return result;
	}
	
	// streaming output routine
	ostream& operator<<(ostream& os, const MathVector& v) {
		for (Index i=0; i<v.size(); i++)
			os << "  " << v[i];
		os << "\n";
		return os;
	}
	
} // namespace PH
