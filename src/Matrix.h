/**
 * @file Matrix.h
 * @brief defines a Matrix class and linear algebra functions.
 *
 * This is an adaptation of Dr. Daniel Reynold's
 * original work, a C++ Matrix class, with modifications
 * for broader compatibility, API consistency,
 * and documentation.
 *
 * The repository for the original work can be found here:
 * [drreynolds/matrix](https://bitbucket.org/drreynolds/matrix) (bitbucket)
 *
 * This file defines the Matrix class, as well a variety of linear
 * algebra functions defined based on matrices and vectors.  Here,
 * "vectors" are considered as MathVector objects.
 *
 * @author Daniel R. Reynolds (drreynolds)
 * @date 2015-06-19
 *
 * @author Paul Herz
 * @date 2016-08-26
 */

#ifndef Matrix_h
#define Matrix_h

// Inclusions
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>

#include "Types.h"
#include "Serializable.h"
#include "Vector.h"
#include "NotImplementedException.h"

namespace PH {
	
	struct Dimensions {
		Index rows;
		Index columns;
		
		Dimensions(Index r, Index c): rows(r), columns(c) {}
		bool operator==(const Dimensions& rhs) const {
			return (this->rows==rhs.rows)
			&& (this->columns==rhs.columns);
		}
		bool operator!=(const Dimensions& rhs) const {
			return !(*this==rhs);
		}
	};

	/**
	 * @class Matrix
	 * @brief an arithmetic matrix class.
	 *
	 * This class wraps a Raw2DArray, alias for std::vector<MathVector>
	 *
	 */
	class Matrix: public Serializable<Matrix> {
	protected:
		
		/// Row count
		Index _rowCount = 0;
		
		/// Column count
		Index _columnCount = 0;
		
	public:
		
		/**
		 * @brief raw matrix data
		 * stored by column.
		 * public for the sake of operator manipulation.
		 */
		Raw2DArray _data;

		/** @defgroup Constructors
			@{
		 */
		#pragma mark Constructors
		
		/// Default constructor for an empty matrix
		Matrix();
		
		/// Size constructor (zero-fills matrix)
		Matrix(Index r, Index c);
		
		Matrix(std::initializer_list<std::initializer_list<double>> l);
		Matrix& operator=(std::initializer_list<std::initializer_list<double>> l);
		
		/// Array wrap constructor using raw array
		template<Index r, Index c>
		static Matrix fromArray(Raw1DArray source);
		
		Matrix(Raw2DArray source);
		Matrix& operator=(const Raw2DArray& source);
		
		static Matrix eye(const Index size);
		static Matrix random(const Index rows, Index columns);
		
		// linear span
		static Matrix linSpace(double a, double b, Index m, Index n);
		
		// logarithmic span
		static Matrix logSpace(double a, double b, Index m, Index n);
		
		/// Copy constructor
		Matrix(const Matrix& A);
		
		/// Copy assignment operator
		Matrix& operator=(const Matrix& A);
		
		/// Move constructor
		Matrix(Matrix&& A);
		
		/// Move assignment operator
		Matrix& operator=(Matrix&& A);
		/// @}
		// end of Constructors
		
		
		/** @defgroup Serializable Implementation
			@{
		 */
		#pragma mark Serializable
		
		void serialize(std::ostream& output) const;
		Matrix deserialize(std::istream& input);
		/// @}
		// end of Serializable implementation
		
		
		void resize(Index r, Index c) {
			_rowCount = r;
			_columnCount = c;
			
			_data.resize(_columnCount);
			for (Index column = 0; column < _columnCount; column++) {
				_data[column].resize(_rowCount);
			}
		}
		
		/** @defgroup Iterators
			@{
		 */
		#pragma mark Iterators
		
		void mapColumns(std::function<void(Raw1DArray&,Index)> callback) {
			for(auto it = _data.begin(); it != _data.end(); ++it) {
				Index i = it - _data.begin();
				callback(*it, i);
			}
		}
		
		void mapElements(std::function<void(MathNumber&,Index,Index)> callback) {
			mapColumns([&callback](Raw1DArray& columnArray, Index c){
				for(auto it = columnArray.begin(); it != columnArray.end(); ++it) {
					Index r = it - columnArray.begin();
					callback(*it, r, c);
				}
			});
		}
		
		/// @}
		// end of Iterators
		
		
		/** @defgroup Dimension Accessors
			@{
		 */
		#pragma mark Dimension Accessors
		
		bool const isSquare() const { return _rowCount == _columnCount; }
		Index const size() const { return _rowCount * _columnCount; }
		Index const& rows() const { return _rowCount; }
		Index const& columns() const { return _columnCount; }
		Dimensions const dimensions() const {
			return Dimensions(_rowCount, _columnCount);
		}
		/// @}
		// end of Dimension Accessors

		
		/** @defgroup Row or Column Accessors
			@{
		 */
		#pragma mark Row or Column Accessors
		
		Raw1DArray column(Index i);
		Raw1DArray row(Index i);
		Raw1DArray& operator[](Index i) {
			// Removed
			throw new NotImplementedException();
		}
		/// @}
		// end of Row or Column Accessors
		
		
		/** @defgroup Cell accessors
			MATLAB/Fortran-style accessors 
			e.g. `M(row,column)` and `M(linearIndex)`
			@{
		 */
		#pragma mark Cell accessors
		MathNumber& operator()(Index row, Index column);
		const MathNumber& operator()(Index row, Index column) const;
		MathNumber& operator()(Index linearIndex);
		const MathNumber& operator()(Index linearIndex) const;
		/// @}
		// end of Cell accessors
		

		/** @defgroup Output routines
			@{
		 */
		#pragma mark Output routines
		
		std::string str(int precision = fullPrecision) const;
		
		friend std::ostream& operator<<(std::ostream& os, const Matrix& A);
		/// @}
		// end of Output routines
		
		
		/** @defgroup In-place operations
			@{
		 */
		#pragma mark In-place operations
		
		/// C = a*A + b*B
		Matrix& linearSumInPlace(MathNumber matrix1Constant, MathNumber b, const Matrix& B);
		
		Matrix& operator+=(const MathNumber constant);
		Matrix& operator+=(const Matrix& matrix);
		
		Matrix& operator-=(const MathNumber constant);
		Matrix& operator-=(const Matrix& matrix);
		
		Matrix& elementwiseMultiply(const Matrix& matrix);
		Matrix& operator*=(const MathNumber constant);
		
		Matrix& elementwiseDivide(const Matrix& matrix);
		Matrix& operator/=(const MathNumber constant);
		
		/// extract/insert routines for portions of vectors
		/// y = x(is,js:ie,je)
		Matrix range(Index beginRow, Index beginColumn, Index endRow, Index endColumn);
		
		/// This(beginRow:beginColumn,endRow:endColumn) = Source
		void insert(const Matrix& source, Index beginRow, Index beginColumn, Index endRow, Index endColumn);
		
		/// Matrix fill operator
		Matrix& operator=(const MathNumber constant);
		
		/// C = C.^p
		Matrix& elementwisePower(const MathNumber constant);
		Matrix& elementwisePower(const Matrix& matrix);
		
		/// Cij = |Cij|
		Matrix abs();
		Matrix& absInPlace();
		
		/// C = C^T
		Matrix& transposeInPlace(); // requires square matrix
		Matrix transpose();
		
		/// C = C^{-1}
		Matrix inverse() const;
		/// @}
		// end of In-place operations
		
		Matrix submatrix(Index beginRow, Index beginColumn,
		 Index endRow, Index endColumn);
		
		// backward substitution on the linear system U*X = B, returning X as a new Matrix
		//    U and B remain unchanged in this operation
		static Matrix backSubstitution(const Matrix& U, const Matrix& B);
		
		// backward substitution on U*x = b, returning x as a new Vector
		//    U and b remain unchanged in this operation
		static Vector backSubstitution(const Matrix& U, const Vector& b);
		
		// forward substitution on the linear system L*X = B, returning X as a new Matrix
		// L and B remain unchanged in this operation
		static Matrix forwardSubstitution(const Matrix& L, const Matrix& B);
		
		// forward substitution on L*x = b, returning x as a new vector<double>
		// L and b remain unchanged in this operation
		static Vector forwardSubstitution(const Matrix& L, const Vector& b);
		
		// solves a linear system A*X = B, returning X as a new Matrix
		// A and B are modified in this operation; X holds the result
		static Matrix linearSolve(Matrix& A, Matrix& B);
		
		// solves a linear system A*x = b, filling in the input vector<double> x
		// A and b are modified in this operation; x holds the result
		static Vector linearSolve(Matrix& A, Vector& b);
		
		static MathNumber dot(const Matrix& A, const Matrix& B);
		
		// Scalar output operators on matrices
		MathNumber min() const;
		MathNumber max() const;
		
		static MathNumber norm(const Matrix& A);
		static MathNumber infNorm(const Matrix& A);
		static MathNumber oneNorm(const Matrix& A);
		
	}; // class Matrix
	
	# pragma mark Functional Operators
	
	Matrix operator*(const Matrix& matrix1, const Matrix& matrix2);
	Vector operator*(const Matrix& matrix, const Vector& vector);
	
	bool operator==(const Matrix& lhs, const Matrix& rhs);
	
	Matrix operator-(const Matrix& matrix, const MathNumber constant);
	Matrix operator-(const MathNumber constant, const Matrix& matrix);
	Matrix operator-(const Matrix& matrix1, const Matrix& matrix2);
	
	bool withinTolerance(const MathNumber a, const MathNumber b, const double precision = 1e-4);

} // namespace PH

#endif // Matrix_h
