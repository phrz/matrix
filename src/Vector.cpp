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
	
	Vector::Vector() {
		_data = Raw1DArray(0);
	}
	
	
	Vector::Vector(Index n) {
		_data = Raw1DArray(n);
	}
	
	
	Vector::Vector(Index n, MathNumber fill) {
		_data = Raw1DArray(n,fill);
	}
	
	
	Vector::Vector(Raw1DArray data) {
		_data = data;
	}
	
	
	Vector::Vector(std::initializer_list<double> l): _data(l) {}
	
	
	Vector Vector::random(Index size) {
		std::random_device randomDevice;
		std::mt19937 generator(randomDevice());
		std::uniform_real_distribution<> dist(0, 1);
		
		Vector result(size);
		result.mapElements([&dist, &generator](MathNumber& element, Index i){
			element = dist(generator);
		});
		
		return result;
	}
	
	
	// create a new vector of linearly spaced data
	Vector Vector::linSpace(MathNumber a, MathNumber b, Index n) {
		if (n<2) {
			throw std::invalid_argument("linSpace expects a vector size argument > 1.");
		}
		
		Vector v(n);
		MathNumber h = (b-a)/(n-1);
		
		for (Index i=0; i<n; i++) {
			// i * h + a
			v[i] = std::fma(i, h, a);
		}
		
		return v;
	}
	
	// create a new vector of logarithmically spaced data
	Vector Vector::logSpace(MathNumber a, MathNumber b, Index n) {
		if (n<2) {
			throw std::invalid_argument("logSpace expects a vector size argument > 1.");
		}
		
		Vector v(n);
		MathNumber h = (b-a)/(n-1);
		
		for (Index i=0; i<n; i++)
			v[i] = pow(10.0, a + i*h);
		
		return v;
	}
	
	
	// Serializable interface
	void Vector::serialize(std::ostream& output) const {
		output << str();
	}
	
	
	Vector Vector::deserialize(std::istream& input) {
		
		auto newData = Raw1DArray();
		std::string line = "";
		
		while (getline(input, line)) {
			
			std::istringstream iss(line);
			double element = 0.0;
			
			while (iss >> element) {
				newData.push_back(element);
			}
		}
		
		return Vector(newData);
	}
	
	
	// extract routine for portions of vectors, output = this(begin:end)
	Vector Vector::range(Index begin, Index end) {
		
		if (this->size() < (end - begin + 1)) {
			// Checking that source vector matches provided subvector range
			throw std::invalid_argument("insert: size mismatch, this vector has " + std::to_string(this->size()) + " entries, but requested subvector has " + std::to_string(end - begin + 1) + " entries.");
		} else if (end >= this->size()) {
			// Checking that range indices are legal
			// (counting numbers, with end in bounds)
			throw std::invalid_argument("insert: requested subvector does not exist: (" + std::to_string(begin) + ":" + std::to_string(end) + ") on a vector with " + std::to_string(this->size()) + " elements.");
		} else if (end < begin) {
			// Checking that lower index is below upper index
			// (no backwards ranges)
			throw std::invalid_argument("insert: requested subvector does not exist, upper index is below lower index: (" + std::to_string(begin) + ":" + std::to_string(end) + ")");
		}
		
		// create new vector of desired size
		Vector result(end - begin + 1);
		
		// copy requested data
		for (Index i = begin; i <= end; i++) {
			result[i-begin] = _data[i];
		}
		
		return result;
	}
	
	
	// insert routine for portions of vectors, this(begin:end) = source
	void Vector::insert(Index begin,
				  Index end, Vector& source) {
		
		if (source.size() != (end - begin + 1)) {
			// Checking that source vector matches provided subvector range
			throw std::invalid_argument("insert: size mismatch, supplied vector has " + std::to_string(source.size()) + " entries, but requested subvector has " + std::to_string(end - begin + 1) + " entries.");
		} else if (end >= this->size()) {
			// Checking that range indices are legal
			// (counting numbers, with end in bounds)
			throw std::invalid_argument("insert: requested subvector does not exist: (" + std::to_string(begin) + ":" + std::to_string(end) + ") on a vector with " + std::to_string(this->size()) + " elements.");
		} else if (end < begin) {
			// Checking that lower index is below upper index
			// (no backwards ranges)
			throw std::invalid_argument("insert: requested subvector does not exist, upper index is below lower index: (" + std::to_string(begin) + ":" + std::to_string(end) + ")");
		}
		
		// perform operation
		for (Index i = 0; i < source.size(); i++) {
			_data[i + begin] = source[i];
		}
	}
	
	
	// In place arithmetic operators (members)
	#pragma mark In-Place Member Operators
	
	/// ADDITION
	
	/// In-place, constant (vector-by-number) addition.
	Vector& Vector::operator+=(const MathNumber addend) {
		this->mapElements([&addend](MathNumber& element, Index i) {
			element += addend;
		});
		
		return *this;
	}
	
	/// In-place, elementwise (vector-by-vector) addition
	Vector& Vector::operator+=(const Vector& addends) {
		if (this->size() != addends.size()) {
			throw std::invalid_argument("vector-by-vector addition: incompatible vector sizes.");
		}
		
		this->mapElements([&addends](MathNumber& element, Index i) {
			element += addends[i];
		});
		
		return *this;
	}
	
	/// SUBTRACTION
	
	/// In-place, constant (vector-by-number) subtraction
	Vector& Vector::operator-=(const MathNumber subtrahend) {
		this->mapElements([&subtrahend](MathNumber& element, Index i) {
			element -= subtrahend;
		});
		
		return *this;
	}
	
	/// In-place, elementwise (vector-by-vector) subtraction
	Vector& Vector::operator-=(const Vector& subtrahends) {
		if (this->size() != subtrahends.size()) {
			throw std::invalid_argument("vector-by-vector subtraction: incompatible vector sizes.");
		}
		
		this->mapElements([&subtrahends](MathNumber& element, Index i) {
			element -= subtrahends[i];
		});
		
		return *this;
	}
	
	/// MULTIPLICATION
	
	/// In-place, constant (vector-by-number) multiplication
	Vector& Vector::operator*=(const MathNumber multiplier) {
		this->mapElements([&multiplier](MathNumber& element, Index i) {
			element *= multiplier;
		});
		
		return *this;
	}
	
	/// In-place, elementwise (vector-by-vector) multiplication
	Vector& Vector::operator*=(const Vector& multipliers) {
		if (this->size() != multipliers.size()) {
			throw std::invalid_argument("vector-by-vector multiplication: incompatible vector sizes.");
		}

		this->mapElements([&multipliers](MathNumber& element, Index i) {
			element *= multipliers[i];
		});
		
		return *this;
	}
	
	/// DIVISION
	
	/// In-place, constant (vector-by-number) division
	Vector& Vector::operator/=(const MathNumber divisor) {
		this->mapElements([&divisor](MathNumber& element, Index i) {
			element /= divisor;
		});
		
		return *this;
	}
	
	/// In-place, elementwise (vector-by-vector) division
	Vector& Vector::operator/=(const Vector& divisors) {
		if (this->size() != divisors.size()) {
			throw std::invalid_argument("vector-by-vector division: incompatible vector sizes.");
		}
		
		this->mapElements([&divisors](MathNumber& element, Index i) {
			element /= divisors[i];
		});
		
		return *this;
	}
	
	/// EXPONENTIATION
	
	/// In-place, constant (vector-by-number) exponentiation
	Vector& Vector::operator^=(const MathNumber exponent) {
		this->mapElements([&exponent](MathNumber& element, Index i) {
			element = std::pow(element, exponent);
		});
		
		return *this;
	}
	
	/// In-place, elementwise (vector-by-vector) exponentiation
	Vector& Vector::operator^=(const Vector& exponents) {
		if (this->size() != exponents.size()) {
			throw std::invalid_argument("vector-by-vector exponentiation: incompatible vector sizes.");
		}
		
		this->mapElements([&exponents](MathNumber& element, Index i){
			element = std::pow(element, exponents[i]);
		});
		
		return *this;
	}
	
	// Vector Dot Product
	// inner product between two vectors
	MathNumber Vector::dot(const Vector& v1, const Vector& v2) {
		if (v1.size() != v2.size()) {
			throw std::invalid_argument("dot product: incompatible vector sizes (must be same)");
		}
		
		MathNumber result = 0.0;
		for(Index i=0; i<v2.size(); i++) {
			result += v1[i] * v2[i];
		}
		
		return result;
	}
	
	// vector 2-norm
	MathNumber Vector::norm(const Vector& v) {
		MathNumber sum = 0.0;
		for(Index i = 0; i < v.size(); i++) {
			// sum = (v[i]*v[i]) + sum
			sum = std::fma(v[i], v[i], sum);
		}
		return sqrt(sum);
	}
	
	
	// vector infinity norm
	MathNumber Vector::infNorm(const Vector& v) {
		MathNumber mx = 0.0;
		for (Index i = 0; i < v.size(); i++) {
			mx = std::max(mx,std::abs(v[i]));
		}
		return mx;
	}
	
	
	// vector one-norm
	MathNumber Vector::oneNorm(const Vector& v) {
		MathNumber sum = 0.0;
		for (Index i = 0; i < v.size(); i++) {
			sum += std::abs(v[i]);
		}
		return sum;
	}
	
	// build a string representation and return it
	std::string Vector::str(int precision) const {
		auto ss = std::stringstream();
		
		for (Index i = 0; i < size(); ++i) {
			ss << "    " << std::setprecision(precision) << (*this)[i];
		}
		
		return ss.str();
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
		result.mapElements([&constant](MathNumber& element, Index i){
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
		result *= vector2;
		return result;
	}
	
	// Result = Vector * -1
	Vector operator-(const Vector& v) {
		return v * -1;
	}
	
	// Result = Vector / Constant
	Vector operator/(const Vector& vector, const MathNumber constant) {
		Vector result = vector;
		result /= constant;
		return result;
	}
	
	// Result = Constant / Vector
	Vector operator/(const MathNumber constant, const Vector& vector) {
		Vector result = vector;
		result.mapElements([&constant](MathNumber& element, Index i){
			element = constant / element;
		});
		return result;
	}
	
	// Result = Vector / Vector
	Vector operator/(const Vector& vector1, const Vector& vector2) {
		Vector result = vector1;
		result /= vector2;
		return result;
	}
	
	// Result = Vector ^ Constant
	Vector operator^(const Vector& vector, const MathNumber constant) {
		Vector result = vector;
		result ^= constant;
		return result;
	}
	
	// Result = Constant ^ Vector
	Vector operator^(const MathNumber constant, const Vector& vector) {
		Vector result = vector;
		result.mapElements([&constant](MathNumber& element, Index i){
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
	std::ostream& operator<<(std::ostream& os, const Vector& v) {
		os << v.str(displayPrecision);
		return os;
	}
	
} // namespace PH
