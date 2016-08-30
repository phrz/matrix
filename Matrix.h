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

#include "PHMath.h"
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
	class Matrix {
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
		
		/// Array wrap constructor using raw array
		template<int r, int c>
		static Matrix& fromArray(Raw1DArray source);
		
		Matrix(Raw2DArray source);
		Matrix& operator=(const Raw2DArray& source);
		
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
			mapColumns([callback](Raw1DArray& columnArray, Index c){
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
		
		Raw1DArray& column(Index i);
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
		double& operator()(Index row, Index column);
		double operator()(Index row, Index column) const;
		double& operator()(Index linearIndex);
		double operator()(Index linearIndex) const;
		/// @}
		// end of Cell accessors
		

		/** @defgroup Output routines
			@{
		 */
		#pragma mark Output routines
		/// write to stdout
		int Write() const;
		
		/// write to file
		int Write(const char *outfile) const;
		
		/// streaming output (cout << M)
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
		
		/// C = C.*A
		Matrix& elementwiseMultiply(const Matrix& other);
		
		/// C = a*C
		Matrix operator*=(const MathNumber a);
		
		/// C = C./A
		Matrix& elementwiseDivide(const Matrix& A);
		
		/// C(is:ie,js:je) = A
		int Insert(const Matrix& A, long int is, long int ie,
		 long int js, long int je);
		
		/// C = a
		int Constant(double a);
		
		/// C = C.^p
		int Power(double p);
		
		/// Cij = |Cij|
		int Abs();
		
		/// C = C^T
		int Trans();
		
		/// C = C^{-1}
		int Inverse();
		/// @}
		// end of In-place operations
		
		// derived matrix creation operations (C is the output, A calls the operation)
		Matrix T();                                              // C = A^T
		Matrix Extract(long int is, long int ie,                 // C = A(is:ie,js:je)
		 long int js, long int je);

		// Scalar output operators on matrices
		double Min() const;                             // min_ij Cij
		double Max() const;                             // min_ij Cij
		
	}; // class Matrix
	
	# pragma mark Functional Operators
	
	Matrix operator*(const Matrix& matrix1, const Matrix& matrix2);
	
	bool operator==(const Matrix& lhs, const Matrix& rhs);

} // namespace PH

#endif // Matrix_h
