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
 * "vectors" are considered as Vector objects.
 *
 * @author Daniel R. Reynolds (drreynolds)
 * @date 2015-06-19
 *
 * @author Paul Herz
 * @date 2016-08-26
 */

#ifndef Matrix_h
#define Matrix_h

#include "Types.h"
#include "Serializable.h"
#include "Vector.h"
#include "NotImplementedException.h"

namespace PH {
	
	
	/**
	 * @brief a struct that wraps a 2D Matrix's size in rows and columns.
	 *
	 * @property rows
	 * @property columns
	 */
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
	 * This class abstracts a container of numeric data, and provides
	 * a variety of factory functions (`eye`, `random`, etc.) and
	 * linear algebra routines.
	 *
	 */
	class Matrix: public Serializable<Matrix> {
	protected:
		
		/// The number of rows.
		Index _rowCount = 0;
		
		/// The number of columns.
		Index _columnCount = 0;
		
	public:
		
		/**
		 * @brief the underlying data structure for class's numeric data.
		 *
		 * This property is a Raw2DArray, or `std::vector<std::vector<double>>`.
		 * It stores the matrix in column-major form, i.e. each inner vector is
		 * a column. This is contrary to how matrices are represented outwardly,
		 * but initializers, accessors, and output routines all take this into
		 * account, so the API consistenly represents location in row-major form.
		 */
		Raw2DArray _data;
		
		
#pragma mark - Constructors
		
		/**
		 * @brief the default constructor, creates an empty (0,0) matrix.
		 */
		Matrix();
		
		
		/**
		 * @brief creates a matrix of a given size, filled with zeroes.
		 *
		 * @param r the number of rows
		 * @param c the number of columns
		 */
		Matrix(Index r, Index c);

		
#pragma mark Initializer list constructors
		
		/**
		 * @brief initializer list constructor, aliases the assignment form.
		 *
		 * Takes a nested (2D) initializer list (e.g. `{{1,2},{3,4},{5,6}}`
		 * as a parameter and builds a matrix if the dimensions are consistent.
		 *
		 * @param l a nested (2D) initializer list of numbers, row-major.
		 *
		 * @throws std::invalid_argument rows in the initializer list 
		 * must have identical length.
		 *
		 * @see Matrix::operator=(std::initializer_list<std::initializer_list<double>> l)
		 */
		Matrix(std::initializer_list<std::initializer_list<double>> l);
		
		
		/**
		 * @brief initializer list assignment constructor
		 *
		 * Handles a nested (2D) initializer list (e.g. `{{1,2},{3,4},{5,6}}`
		 * as a rvalue and builds a matrix if the dimensions are consistent.
		 *
		 * @param l a nested (2D) initializer list of numbers, row-major.
		 *
		 * @throws std::invalid_argument rows in the initializer list
		 * must have identical length.
		 */
		Matrix& operator=(std::initializer_list<std::initializer_list<double>> l);
		
		
#pragma mark Copy constructors
		
		/**
		 * @brief matrix copy constructor, aliases copy assignment operator.
		 *
		 * sets the matrix size to the source matrix's size and copies its
		 * data structure over.
		 *
		 * @param source a source matrix from which to copy.
		 *
		 * @see Matrix::operator=(const Matrix& source)
		 */
		Matrix(const Matrix& source);
		
		
		/**
		 * @brief matrix copy assignment operator
		 *
		 * sets the matrix size to the source matrix's size and copies its
		 * data structure over.
		 *
		 * @param source another matrix from which to nondestructively 
		 * copy data.
		 *
		 * @return a reference to the newly copied matrix (`*this`)
		 */
		Matrix& operator=(const Matrix& source);
		
		
#pragma mark Move constructors
		
		/**
		 * @brief matrix move constructor, aliases move assignment operator.
		 *
		 * sets the matrix size to the source matrix's size and copies its
		 * data structure over.
		 *
		 * @param source a source matrix from which to move the data.
		 *
		 * @see Matrix::operator=(const Matrix&& source)
		 */
		Matrix(Matrix&& source);
		
		
		/**
		 * @brief matrix move assignment operator
		 *
		 * sets the matrix size to the source matrix's size and copies its
		 * data structure over.
		 *
		 * @param source a source matrix from which to move the data.
		 */
		Matrix& operator=(Matrix&& source);
		
		
#pragma mark Type conversion constructors
		
		/**
		 * @brief array wrap constructor
		 *
		 * Given dimensions and a one-dimensional array of numbers, this
		 * constructor attempts to linearly fill the 2D matrix from the 1D
		 * source data. Fails if the size of the matrix is not equivalent
		 * to the length of the 1D array given.
		 *
		 * @tparam r row size of the created matrix.
		 * @tparam c column size of the created matrix.
		 *
		 * @param source a one-dimensional numeric array of source data.
		 *
		 * @throws std::invalid_argument the size of the matrix must equal the
		 * length of the source data.
		 *
		 * @return a Matrix of size (r,c) filled with the given source data.
		 */
		template<Index r, Index c>
		static Matrix fromArray(Raw1DArray source) {
			if (r * c != source.size()) {
				throw std::invalid_argument("Matrix::fromArray: the given array could not properly fit into the requested matrix dimensions.");
			}
			
			Matrix result = Matrix(r,c);
			Index columnCount = result.columns();
			
			result.mapElements([&source, &columnCount]
			(MathNumber& element, Index row, Index column) {
				element = source[row * columnCount + column];
			});
			
			return result;
		}
		
		
		/**
		 * @brief constructor, aliases row-major 2D array assignment operator.
		 *
		 * Given a `Raw2DArray` (`std::vector<std::vector<double>>`)
		 * representing a row-major matrix (human-readable form) as a parameter,
		 * create a Matrix containing this data. Fails if row length is not
		 * consistent.
		 *
		 * @param source a nested (2D) array in row-major form
		 * (i.e. each inner array is a row)
		 *
		 * @throws std::invalid_argument rows in the 2D array must have
		 * identical length.
		 *
		 * @see Matrix::operator=(const Raw2DArray& source)
		 */
		Matrix(const Raw2DArray& source);
		
		
		/**
		 * @brief row-major 2D array assignment operator
		 *
		 * Given a `Raw2DArray` (`std::vector<std::vector<double>>`)
		 * representing a row-major matrix (human-readable form) as an rvalue,
		 * create a Matrix containing this data. Fails if row length is not
		 * consistent.
		 *
		 * @param source a nested (2D) array in row-major form
		 * (i.e. each inner array is a row)
		 *
		 * @throws std::invalid_argument rows in the 2D array must have
		 * identical length.
		 *
		 * @return a reference to the newly constructed matrix (`*this`)
		 */
		Matrix& operator=(const Raw2DArray& source);
		
		
#pragma mark - Convenience factories
		
		/**
		 * @brief creates an identity matrix of a given size (square)
		 *
		 * an identity matrix is a square matrix of size (n,n) that is all
		 * zero-value except for its diagonals ((i,i) for all 0...n-1),
		 * which are ones.
		 *
		 * @param size the size of the identity matrix, i.e. both its row and
		 * column size.
		 *
		 * @return the generated identity matrix
		 */
		static Matrix eye(const Index size);
		
		
		/**
		 * @brief creates a randomly-filled matrix of a given size.
		 *
		 * creates a matrix of a given size and sets all elements to a random
		 * value in a uniform distribution 0.0...1.0. This is much like the
		 * random matrix generator in MATLAB in that it is not merely
		 * pseudorandom but on a uniform distribution.
		 *
		 * @param rows the desired row count for the random matrix.
		 * @param columns the desired column count for the random matrix.
		 *
		 * @return the generated uniform distribution matrix
		 */
		static Matrix random(const Index rows, Index columns);
		
		
		/**
		 * @brief creates a linear span from a...b along a matrix.
		 *
		 * Given numbers a and b to form a range, this creates a linear span
		 * (constant spacing through the range) along the linear indices
		 * of a matrix.
		 *
		 * Layman's terms: first creates a array of size rows*columns
		 * and fills it with values a...b, incrementing the value
		 * by a constant that will get it to `b` by the end of the array.
		 * Then it linearly fills the created matrix (size (rows,columns))
		 * with these values, going left-to-right and top-to-bottom.
		 *
		 * @param a the lower bound of the linear span.
		 * @param b the upper bound of the linear span.
		 * @param rows the desired row count for the linear span matrix.
		 * @param columns the desired row count for the linear span matrix.
		 *
		 * @return the generated linear span matrix.
		 */
		static Matrix linSpace(double a, double b, Index rows, Index columns);
		
		
		/**
		 * @brief creates a logarithmic span from a...b along a matrix.
		 *
		 * Given numbers a and b to form a range, this creates a logarithmic span
		 * (logarithmic spacing through the range) along the linear indices
		 * of a matrix.
		 *
		 * Layman's terms: first creates a array of size rows*columns
		 * and fills it with values 10^a...10^b, incrementing the exponent
		 * by a constant that will get it to `b` by the end of the array.
		 * Then it linearly fills the created matrix (size (rows,columns))
		 * with these values, going left-to-right and top-to-bottom.
		 *
		 * @param a the exponent of the lower bound of the logarithmic span.
		 * @param b the exponent of the upper bound of the logarithmic span.
		 * @param rows the desired row count for the logarithmic span matrix.
		 * @param columns the desired row count for the logarithmic span matrix.
		 *
		 * @return the generated logarithmic span matrix.
		 */
		static Matrix logSpace(double a, double b, Index rows, Index columns);
		
		
#pragma mark - Serializable implementation
		
		/**
		 * @brief represents the matrix in full precision as a stream.
		 *
		 * This method is required by Serializable in order to serialize
		 * the contents of the Matrix to full precision. This function is
		 * rarely called directly, instead saveTo(path) is used for File I/O.
		 *
		 * Uses the str() method for standardized string representation.
		 *
		 * @param output the output stream to pipe out the serialized data to.
		 *
		 * @see Serializable::saveTo
		 * @see Matrix::str
		 */
		void serialize(std::ostream& output) const;
		
		
		/**
		 * @brief converts a serialized Matrix string back to a Matrix object.
		 *
		 * This method is required by Serializable in order to deserialize
		 * a serialized matrix string back to a Matrix object. This function is
		 * rarely called directly, instead loadFrom(path) is used for File I/O.
		 *
		 * Deserializes based on the standard representation in the str() function.
		 *
		 * @param input the input stream from which to receive serialized data.
		 *
		 * @see Serializable::loadFrom
		 * @see Matrix::str
		 */
		Matrix deserialize(std::istream& input);

		
#pragma mark - Iterators
		
		/**
		 * @brief iterates over columns with a lambda function.
		 *
		 * @param callback a function repeatedly called for each column, with
		 * two parameters: a `Raw1DArray` reference (the column data), and
		 * an `Index` number (the column number)
		 *
		 */
		void mapColumns(std::function<void(Raw1DArray&,Index)> callback) {
			for(auto it = _data.begin(); it != _data.end(); ++it) {
				Index i = it - _data.begin();
				callback(*it, i);
			}
		}
		
		
		/**
		 * @brief iterates over all elements with a lambda function.
		 *
		 * @param callback a function repeatedly called for each element, with
		 * three parameters: a `MathNumber` reference (the element data), and
		 * two Index parameters, the row and column indices respectively.
		 *
		 */
		void mapElements(std::function<void(MathNumber&,Index,Index)> callback) {
			mapColumns([&callback](Raw1DArray& columnArray, Index c){
				for(auto it = columnArray.begin(); it != columnArray.end(); ++it) {
					Index r = it - columnArray.begin();
					callback(*it, r, c);
				}
			});
		}
		
		
#pragma mark - Dimension accessors and mutators
		
		/**
		 * @brief resets row and column count properties and resizes the data structure.
		 *
		 * @param r the new row count
		 * @param c the new column count
		 *
		 */
		void resize(Index r, Index c) {
			_rowCount = r;
			_columnCount = c;
			
			_data.resize(_columnCount);
			for (Index column = 0; column < _columnCount; column++) {
				_data[column].resize(_rowCount);
			}
		}
		
		
		/**
		 * @brief quick square matrix checker
		 *
		 * @return whether the matrix is square or not.
		 */
		bool const isSquare() const { return _rowCount == _columnCount; }
		
		
		/**
		 * @brief quick size calculator
		 *
		 * @return the element count of the matrix, effectively rows * columns.
		 */
		Index const size() const { return _rowCount * _columnCount; }
		
		
		/**
		 * @brief row count accessor
		 *
		 * @return the number of rows in the matrix.
		 */
		Index const& rows() const { return _rowCount; }
		
		
		/**
		 * @brief column count accessor
		 *
		 * @return the number of column in the matrix.
		 */
		Index const& columns() const { return _columnCount; }
		
		
		/**
		 * @brief dimensions struct generator
		 *
		 * @return a Dimensions struct containing the row and column counts.
		 */
		Dimensions const dimensions() const {
			return Dimensions(_rowCount, _columnCount);
		}
		

#pragma mark - Data accessors
		
		/**
		 * @brief matrix column copier
		 *
		 * @return a copy of the nth column.
		 */
		Raw1DArray copyColumn(Index i);
		
		
		/**
		 * @brief matrix row copier
		 *
		 * @return a copy of the nth row.
		 */
		Raw1DArray copyRow(Index i);
		
		
		/// Deprecated column accessor, removed because of
		/// ambiguous syntax (confused with linear index accessor)
		Raw1DArray& operator[](Index i) {
			// Removed
			throw NotImplementedException();
		}
		
		
		/**
		 * @brief 2D position element accessor (mutable)
		 *
		 * @return a mutable reference to the element at (row,column).
		 */
		MathNumber& operator()(Index row, Index column);
		
		
		/**
		 * @brief 2D position element accessor
		 *
		 * @return an immutable reference to the element at (row,column).
		 */
		const MathNumber& operator()(Index row, Index column) const;
		
		
		/**
		 * @brief 2D position element accessor (mutable)
		 *
		 * @return a mutable reference to the element at a given linear index.
		 */
		MathNumber& operator()(Index linearIndex);
		
		
		/**
		 * @brief 2D position element accessor (immutable)
		 *
		 * @return an immutable reference to the element at a given linear index.
		 */
		const MathNumber& operator()(Index linearIndex) const;
		
		
#pragma mark - Output routines
		
		/**
		 * @brief string representation generator
		 *
		 * @param precision the precision to which to display the numbers.
		 * in `operator<<`, this defaults to `displayPrecision`. Elsewhere,
		 * most importantly in serialization to file, this defaults to 
		 * `fullPrecision`.
		 *
		 * @return a string representation of the matrix.
		 */
		std::string str(int precision = fullPrecision) const;
		
		
		/**
		 * @brief stream insertion operator, aliases str()
		 *
		 * calls `str(precision)`, specifying a precision of `displayPrecision`
		 *
		 * @param os the output stream (left hand side)
		 * @param A the matrix to insert in the stream (right hand side)
		 *
		 * @see Matrix::str
		 */
		friend std::ostream& operator<<(std::ostream& os, const Matrix& A);
		
		
#pragma mark - In-place transformations
		
		/// in-place elementwise sum of this matrix (times a constant) and
		/// another matrix (times a constant)
		Matrix& linearSumInPlace(MathNumber matrix1Constant,
								 MathNumber matrix2Constant,
								 const Matrix& matrix2);
		
		
#pragma mark Addition
		
		/// in-place sum of this matrix and a constant.
		Matrix& operator+=(const MathNumber constant);
		
		
		/// in-place elementwise sum of this matrix and another.
		Matrix& operator+=(const Matrix& matrix);
		
		
#pragma mark Subtraction
		
		/// in-place difference of this matrix and a constant.
		Matrix& operator-=(const MathNumber constant);
		
		
		/// in-place elementwise difference of this matrix and another.
		Matrix& operator-=(const Matrix& matrix);
		
		
#pragma mark Multiplication
		
		/// in-place elementwise product of this matrix and another.
		///
		/// explicitly named to avoid confusion with dot product or
		/// non-elementwise multiplication (actual matrix multiplication)
		/// which is not an in-place operation.
		Matrix& elementwiseMultiply(const Matrix& matrix);
		
		
		/// in-place product of this matrix and a constant.
		Matrix& operator*=(const MathNumber constant);
		
		
#pragma mark Division
		
		/// in-place elementwise quotient of this matrix and another.
		///
		/// explicitly named to avoid confusion with "right division" or
		/// "left division", entirely different operations and typically
		/// performed or at least written differently (right division A/B is
		/// written AB^-1, which can be performed as `A*B.inverse()`
		Matrix& elementwiseDivide(const Matrix& matrix);
		
		
		/// in-place quotient of this matrix and a constant.
		Matrix& operator/=(const MathNumber constant);

		
#pragma mark Power
		
		/// in-place exponentiation by a constant.
		///
		/// explicitly named to avoid confusion with actual matrix exponentials,
		/// a different operation altogether, which can only be performed on
		/// square matrices.
		Matrix& elementwisePower(const MathNumber constant);
		
		
		/// in-place elementwise exponentiation by another matrix.
		///
		/// explicitly named to avoid confusion with different exponential
		/// operations involving matrices.
		Matrix& elementwisePower(const Matrix& matrix);
		
		
#pragma mark - Submatrix operations
		
		/**
		 * @brief copies a submatrix from this matrix.
		 *
		 * given four parameters, the first pair being a 2D index in this
		 * matrix, and the second pair being another 2D index in this matrix,
		 * create a new matrix ("submatrix") and fill it with the values between
		 * the two 2D indices, effectively "copying" a rectangular section of
		 * this matrix.
		 *
		 * @param beginRow the row dimension of the starting 2D index
		 * @param beginColumn the column dimension of the starting 2D index
		 * @param endRow the row dimension of the ending 2D index
		 * @param endColumn the column dimension of the ending 2D index
		 *
		 * @throws std::invalid_argument one or more 2D index components exceed
		 * the actual dimensions of this matrix.
		 * @throws std::invalid_argument the intended upper 2D index occurs
		 * before thelower 2D index - backwards ranges are not supported.
		 * @throws std::invalid_argument cannot build a submatrix of size (0,0)
		 *
		 * @return the generated submatrix copy
		 */
		Matrix range(Index beginRow, Index beginColumn, Index endRow, Index endColumn);
		
		
		/**
		 * @brief copies a submatrix from this matrix.
		 *
		 * given five parameters, the first being a "source" matrix,
		 * the following pair being a 2D index in *this*
		 * matrix, and the last pair being another 2D index in *this* matrix,
		 * fill in the rectangular area of this matrix defined by the two indices
		 * (top-left and bottom-right) with all the values of the source matrix,
		 * effectively "pasting" it into this matrix. Fails if the source matrix
		 * is not identical in size to the rectangular area of this matrix.
		 *
		 * @param beginRow the row dimension of the starting 2D index
		 * @param beginColumn the column dimension of the starting 2D index
		 * @param endRow the row dimension of the ending 2D index
		 * @param endColumn the column dimension of the ending 2D index
		 *
		 * @throws std::invalid_argument the source matrix is not identical
		 * in size to the rectangular area defined for this matrix.
		 * @throws std::invalid_argument one or more 2D index components exceed
		 * the actual dimensions of this matrix.
		 * @throws std::invalid_argument the intended upper 2D index occurs
		 * before thelower 2D index - backwards ranges are not supported.
		 *
		 */
		void insert(const Matrix& source, Index beginRow, Index beginColumn, Index endRow, Index endColumn);

		
#pragma mark - Other operations
		
		/// Matrix fill operator
		///
		/// fills all elements with a constant.
		Matrix& operator=(const MathNumber constant);
		
		
		/// Replaces all elements with their absolute values.
		Matrix& absInPlace();
		
		
		/// Creates a new matrix where all elements are absolute values.
		Matrix abs();
		
		
		/// Transposes a square matrix in place.
		Matrix& transposeInPlace();
		
		
		/// Creates a new matrix that is this matrix's transpose.
		Matrix transpose();
		
		
		/// Creates a new matrix that is this matrix's inverse.
		///
		/// Fails if this matrix is square or singular.
		Matrix inverse() const;
		
		
		/// backward substitution on the linear system U*X = B, returning X as a new Matrix
		/// U and B remain unchanged in this operation
		static Matrix backSubstitution(const Matrix& U, const Matrix& B);
		
		
		/// backward substitution on U*x = b, returning x as a new Vector
		/// U and b remain unchanged in this operation
		static Vector backSubstitution(const Matrix& U, const Vector& b);
		
		
		/// forward substitution on the linear system L*X = B, returning X as a new Matrix
		/// L and B remain unchanged in this operation
		static Matrix forwardSubstitution(const Matrix& L, const Matrix& B);
		
		
		/// forward substitution on L*x = b, returning x as a new vector<double>
		/// L and b remain unchanged in this operation
		static Vector forwardSubstitution(const Matrix& L, const Vector& b);
		
		
		/// solves a linear system A*X = B, returning X as a new Matrix
		/// A and B are modified in this operation; X holds the result
		static Matrix linearSolve(Matrix& A, Matrix& B);
		
		
		/// solves a linear system A*x = b, filling in the input vector<double> x
		/// A and b are modified in this operation; x holds the result
		static Vector linearSolve(Matrix& A, Vector& b);
		
		
		/// Matrix dot product
		static MathNumber dot(const Matrix& A, const Matrix& B);
		
		
#pragma mark Scalar output methods
		
		/// gets minimum value element
		MathNumber min() const;
		
		
		/// gets maximum value element
		MathNumber max() const;
		
		
		/// matrix Frobenius norm (column/row vector 2-norm)
		static MathNumber norm(const Matrix& A);
		
		
		/// matrix infinity norm (column vector infinity norm, row vector one norm)
		static MathNumber infNorm(const Matrix& A);
		
		
		/// matrix one norm (column vector one norm, row vector infinity norm)
		static MathNumber oneNorm(const Matrix& A);
		
		
	}; // class Matrix
	
# pragma mark Infix operators
	
	Matrix operator*(const Matrix& matrix1, const Matrix& matrix2);
	Vector operator*(const Matrix& matrix, const Vector& vector);
	
	bool operator==(const Matrix& lhs, const Matrix& rhs);
	
	Matrix operator-(const Matrix& matrix, const MathNumber constant);
	Matrix operator-(const MathNumber constant, const Matrix& matrix);
	Matrix operator-(const Matrix& matrix1, const Matrix& matrix2);
	
	bool withinTolerance(const MathNumber a, const MathNumber b, const double precision = 1e-4);

} // namespace PH

#endif // Matrix_h
