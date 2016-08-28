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

#ifndef MATRIX_DEFINED__
#define MATRIX_DEFINED__

// Inclusions
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <experimental/optional>

// for accessing pi via M_PI
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <cmath>
#else
#include <cmath>
#endif

namespace PH {

	using  MathNumber = double;
	using  Raw1DArray = std::vector<MathNumber>;
	using  Raw2DArray = std::vector<Raw1DArray>;
	using       Index = std::size_t;

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
		
		template<typename... ParameterTypes>
		void mapColumns(void (*callback)(ParameterTypes...)) {
			
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
		/// @}
		// end of Dimension Accessors

		
		/** @defgroup Row or Column Accessors
			@{
		 */
		#pragma mark Row or Column Accessors
		
		MathVector& column(Index i);
		MathVector row(Index i);
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
			C is the matrix calling the routine, i.e. C = A*X is C.Product(A,X)
			zero is success, 1 is failure.
			@{
		 */
		#pragma mark In-place operations
		
		/// C = A*X
		int Product(const Matrix& A, const Matrix& X);
		
		/// C = a*A + b*B
		int LinearSum(MathNumber a, const Matrix& A, MathNumber b, const Matrix& B);
		
		/// C = C+A
		int Add(const Matrix& A) { return LinearSum(1.0,*this,1.0,A); };
		
		/// C = C+a
		int Add(MathNumber a);
		
		/// C = C-A
		int Subtract(const Matrix& A) { return LinearSum(1.0,*this,-1.0,A); };
		
		/// C = C-a
		int Subtract(MathNumber a) { return Add(-a); };
		
		/// C = C.*A
		int Multiply(const Matrix& A);
		
		/// C = a*C
		int Multiply(MathNumber a);
		
		/// C = C./A
		int Divide(const Matrix& A);
		
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
		
		
		/** @defgroup In-place operation aliases
			C is the matrix calling the routine, i.e. C = A*X is C.Product(A,X)
			zero is success, 1 is failure.
			@{
		 */
		#pragma mark In-place operation aliases
		
		/// C = C+A
		int operator+=(const Matrix& A) {
			return Add(A);
		}
		
		/// C = C+a (add a to all entries)
		int operator+=(double a)        { return Add(a); }
		
		/// C = C-A
		int operator-=(const Matrix& A) { return Subtract(A); }
		
		/// C = C-a (add -a to all entries)
		int operator-=(double a)        { return Subtract(a); };
		
		/// C = C.*A (componentwise multiply)
		int operator*=(const Matrix& A) { return Multiply(A); };
		
		/// C = a*C
		int operator*=(double a)        { return Multiply(a); };
		
		/// C = C./A (componentwise divide)
		int operator/=(const Matrix& A) { return Divide(A); };
		
		/// C = (1/a)*C
		int operator/=(double a)        { return Multiply(1.0/a); };
		
		/// C = C.^p (componentwise power)
		int operator^=(double p)        { return Power(p); };
		
		/// C = a (all entries equal a)
		int operator=(double a)         { return Constant(a); };
		/// @}
		// end of In-place operation aliases
		
		
		// derived matrix creation operations (C is the output, A calls the operation)
		Matrix T();                                              // C = A^T
		Matrix Extract(long int is, long int ie,                 // C = A(is:ie,js:je)
		 long int js, long int je);

		// Scalar output operators on matrices
		double Min() const;                             // min_ij Cij
		double Max() const;                             // min_ij Cij
		bool operator==(const Matrix& A) const;         // check for Matrix equality
	};



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
	Matrix Random(Index m, Index n);                         // Cij random in [0,1]
	Matrix Eye(Index n);                                      // Cij = delta_i(j)
	Matrix MatrixRead(const char *infile);                     // creates from input file
	Matrix Inverse(const Matrix& A);                           // C = A^{-1}



	//--- supplementary matrix-vector arithmetic routines ---

	// standard matrix-vector product -> new vector (function form)
	MathVector matrixVectorProduct(const Matrix& A, const MathVector& v);

	// standard matrix-vector product -> new vector
	MathVector operator*(const Matrix& A, const MathVector& v);


	//--- linear algebra routines ---

	// backward substitution on the linear system U*X = B, filling in an existing Matrix X
	//    U and B remain unchanged in this operation; X holds the result
	//    B and X may have multiple columns
	int backSubstitution(const Matrix& U, Matrix& X, const Matrix& B);

	// backward substitution on the linear system U*X = B, returning X as a new Matrix
	//    U and B remain unchanged in this operation
	Matrix backSubstitution(const Matrix& U, const Matrix& B);

	// backward substitution on U*x = b, filling in an existing vector<double> x
	//    U and b remain unchanged in this operation; x holds the result
	int backSubstitution(const Matrix& U, MathVector& x, const MathVector& b);

	// backward substitution on U*x = b, returning x as a new vector<double>
	//    U and b remain unchanged in this operation
	MathVector backSubstitution(const Matrix& U, const MathVector& b);

	// forward substitution on the linear system L*X = B, returning X as a new Matrix
	// L and B remain unchanged in this operation
	std::optional<Matrix> forwardSubstitution(const Matrix& L, const Matrix& B);

	// forward substitution on L*x = b, returning x as a new vector<double>
	// L and b remain unchanged in this operation
	std::optional<Vector> ForwardSubstitution(const Matrix& L, const Vector& b)

	// solves a linear system A*X = B, returning X as a new Matrix
	// A and B are modified in this operation; X holds the result
	std::optional<Matrix> linearSolve(Matrix& A, Matrix& B)

	// solves a linear system A*x = b, filling in the input vector<double> x
	// A and b are modified in this operation; x holds the result
	std::optional<Vector> linearSolve(Matrix& A, Vector& b)

	// solves a linear system A*x = b, returning x as a new vector<double>
	// A and b are modified in this operation; x holds the result
	MathVector linearSolve(Matrix& A, MathVector& b);

} // namespace PH

#endif // MATRIX_DEFINED__
