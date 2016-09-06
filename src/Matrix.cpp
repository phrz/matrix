#include "Matrix.h"

namespace PH {
	
	
#pragma mark - Constructors
	
	Matrix::Matrix() {
		this->resize(0,0);
	}
	
	
	Matrix::Matrix(Index r, Index c) {
		this->resize(r,c);
	}
	
	
#pragma mark Initializer list constructors
	
	Matrix::Matrix(std::initializer_list<std::initializer_list<double>> l) {
		(*this) = l;
	}
	
	
	Matrix& Matrix::operator=(std::initializer_list<std::initializer_list<double>> l) {
		// row major array to column major data structure
		Index listRowCount = l.size();
		Index listColumnCount = l.begin()->size();
		this->resize(listRowCount, listColumnCount);
		
		// iterate rows
		Index row = 0;
		for(auto listRow: l) {
			if(listRow.size() != listColumnCount) {
				throw new std::invalid_argument("Matrix initializer list constructor: rows in the initializer list must have identical length.");
			}
			
			// iterate columns within row
			Index column = 0;
			for(auto listRowElement: listRow) {
				(*this)(row, column) = listRowElement;
				++column;
			}
			
			++row;
		} // row iterator
		
		return *this;
	}
	
	
#pragma mark Copy constructors
	
	Matrix::Matrix(const Matrix& source) {
		(*this) = source;
	}
	
	
	Matrix& Matrix::operator=(const Matrix& source) {
		this->resize(source.rows(), source.columns());
		this->_data = source._data;
		
		return *this;
	}
	
	
#pragma mark Move constructors
	
	Matrix::Matrix(Matrix&& source) {
		(*this) = std::move(source);
	}
	
	
	Matrix& Matrix::operator=(Matrix&& source) {
		if(this != &source) {
			this->resize(source.rows(), source.columns());
			
			_data = std::move(source._data);
			
			source.resize(0,0);
			source._data = Raw2DArray();
		}
		
		return *this;
	}
	
	
#pragma mark Type conversion constructors
	
	Matrix::Matrix(const Raw2DArray& source) {
		(*this) = source;
	}
	
	
	Matrix& Matrix::operator=(const Raw2DArray& source) {
		// row major matrix constructor
		this->resize(source.size(), source[0].size());
		
		for (Index row = 0; row < _rowCount; ++row) {
			if (source[row].size() != _columnCount) {
				throw new std::invalid_argument("Matrix Raw2DArray constructor: rows in the 2D array must have identical length. Make sure the array provided is row-major.");
			}
			
			for (Index column = 0; column < _columnCount; ++column) {
				(*this)(row, column) = source[row][column];
			}
		}
		
		return *this;
	}
	
	
#pragma mark - Convenience factories
	
	Matrix Matrix::eye(const Index size) {
		Matrix I(size, size);
		for (Index i = 0; i < size; i++) {
			I(i,i) = 1.0;
		}
		return I;
	}
	
	
	Matrix Matrix::random(Index rows, Index columns) {
		std::random_device randomDevice;
		std::mt19937 generator(randomDevice());
		std::uniform_real_distribution<> dist(0, 1);
		
		Matrix result(rows, columns);
		result.mapElements([&dist, &generator](MathNumber& element, Index r, Index c){
			element = dist(generator);
		});
		
		return result;
	}


	Matrix Matrix::linSpace(MathNumber a, MathNumber b, Index rows, Index columns) {
		Matrix result(rows, columns);
		
		if(b < a) {
			throw new std::invalid_argument("linSpace: b (upper bound) less than a (lower bound)");
		}
		
		MathNumber range = b - a;
		MathNumber delta = range / (rows * columns - 1);
		
		result.mapElements([&](MathNumber& element, Index r, Index c) {
			element = a + (r * columns + c) * delta;
		});
		
		return result;
	}
	
	
	Matrix Matrix::logSpace(MathNumber a, MathNumber b, Index rows, Index columns) {
		Matrix result(rows, columns);
		
		if(b < a) {
			throw new std::invalid_argument("linSpace: b (upper bound) less than a (lower bound)");
		}
		
		MathNumber range = b - a;
		MathNumber delta = range / (rows * columns - 1);
		
		result.mapElements([&](MathNumber& element, Index r, Index c) {
			element = std::pow(10.0, a + (r * columns + c) * delta);
		});
		
		return result;
	}
	
	
#pragma mark - Serializable implementation
	
	void Matrix::serialize(std::ostream& output) const {
		output << str();
	}
	
	
	Matrix Matrix::deserialize(std::istream& input) {
		
		// determine matrix size
		Index newRowCount = 0, newColumnCount = 0, currentRow = 0, currentColumn = 0;
		auto newData = Raw2DArray();
		std::string line = "";
		
		while (getline(input, line)) {
			std::istringstream iss(line);
			double element = 0.0;
			currentColumn = 0;
			
			newData.push_back(Raw1DArray());
			auto* row = &(newData.back());
			
			while (iss >> element) {
				row->push_back(element);
				element = 0.0;
				++currentColumn;
			}
			
			if((currentColumn != 0) && (newRowCount == 0)) {
				// first row, set column number
				newColumnCount = currentColumn;
			}
			
			if ((currentColumn != 0) && (currentColumn != newColumnCount)) {
				throw new std::runtime_error("matrix deserialize: not all rows had same column count.");
			}
			
			if (currentColumn != 0) {
				++currentRow;
			} else {
				// pop the empty row array we added
				newData.pop_back();
			}
		}
		
		// last row number is total row count
		newRowCount = currentRow;
		
		// create matrix of desired size
		Matrix result = Matrix(newData);
		
		return result;
	}


#pragma mark - Data accessors
	
	Raw1DArray Matrix::copyColumn(Index i) {
	  return _data[i];
	}
	
	
	Raw1DArray Matrix::copyRow(Index row) {
		Raw1DArray rowArray = Raw1DArray(_columnCount);
		
		for (Index column = 0; column < _columnCount; column++) {
			rowArray[column] = _data[column][row];
		}
		
		return rowArray;
	}

	
	MathNumber& Matrix::operator()(Index row, Index column) {
		return _data[column][row];
	}
	
	
	const MathNumber& Matrix::operator()(Index row, Index column) const {
		return _data[column][row];
	}
	
	
	MathNumber& Matrix::operator()(Index linearIndex) {
		return _data[linearIndex / _rowCount][linearIndex % _rowCount];
	}
	
	
	const MathNumber& Matrix::operator()(Index linearIndex) const {
		return _data[linearIndex / _rowCount][linearIndex % _rowCount];
	}
	

#pragma mark - Output routines
	
	std::string Matrix::str(int precision) const {
		auto ss = std::stringstream();
		
		for (Index row = 0; row < _rowCount; row++) {
			for (Index column = 0; column < _columnCount; column++) {
				ss << "    " << std::setprecision(precision) << (*this)(row, column);
			}
			ss << std::endl;
		}
		
		return ss.str();
	}

	
	std::ostream& operator<<(std::ostream& os, const Matrix& A) {
		os << A.str(displayPrecision);
		return os;
	}
	
	
# pragma mark - In-place transformations
	
	Matrix& Matrix::linearSumInPlace(MathNumber matrix1Constant, MathNumber matrix2Constant, const Matrix& matrix2) {
		// check that array sizes match
		if (this->dimensions() != matrix2.dimensions()) {
			throw new std::invalid_argument("matrix-by-matrix linear sum (a*A+b*B): incompatible matrix dimensions (must be same)");
		}
		
		mapElements([&matrix1Constant, &matrix2Constant, &matrix2]
					(MathNumber& matrix1Element, Index row, Index column)
		{
			matrix1Element *= matrix1Constant;
			// fused multiply-add operation ($0 * $1) + $2
			// usually handled in one operation on a modern processor.
			matrix1Element = std::fma(matrix2Constant, matrix2(row, column), matrix1Element);
		});
		
		return *this;
	}


#pragma mark Addition
	
	Matrix& Matrix::operator+=(const MathNumber constant) {
		mapElements([&constant](MathNumber& element, Index row, Index column) {
			element += constant;
		});
		
		return *this;
	}
	
	
	Matrix& Matrix::operator+=(const Matrix& matrix) {
		linearSumInPlace(1.0, 1.0, matrix);
		
		return *this;
	}
	
	
#pragma mark Subtraction
	
	Matrix& Matrix::operator-=(const MathNumber constant) {
		(*this) += -constant;
		return *this;
	}
	
	
	Matrix& Matrix::operator-=(const Matrix& matrix) {
		linearSumInPlace(1.0, -1.0, matrix);
		return *this;
	}
	

#pragma mark Multiplication
	
	Matrix& Matrix::elementwiseMultiply(const Matrix& other) {
		// check that array sizes match
		if (this->dimensions() != other.dimensions()) {
			throw new std::invalid_argument("elementwise matrix multiplication (A.*B): incompatible matrix dimensions (must be same)");
		}

		// perform operation
		mapElements([&other](MathNumber& element, Index row, Index column) {
			element *= other(row, column);
		});

		return *this;
	}
	
	
	Matrix& Matrix::operator*=(const MathNumber constant) {
		// perform operation
		mapElements([&constant](MathNumber& element, Index row, Index column) {
			element *= constant;
		});

		return *this;
	}
	
	
#pragma mark Division
	
	Matrix& Matrix::elementwiseDivide(const Matrix& matrix) {
		// check that array sizes match
		if (this->dimensions() != matrix.dimensions()) {
			throw new std::invalid_argument("elementwise matrix division (A./B): incompatible matrix dimensions (must be same)");
		}
		
		mapElements([&matrix](MathNumber& element, Index row, Index column) {
			MathNumber otherElement = matrix(row,column);
			if (otherElement == 0.0) {
				throw new std::runtime_error("elementwise matrix division (A./B): division by zero");
			}
			element /= otherElement;
		});
		
		return *this;
	}
	
	
	Matrix& Matrix::operator/=(const MathNumber constant) {
		// perform operation
		mapElements([&constant](MathNumber& element, Index row, Index column) {
			element /= constant;
		});
		
		return *this;
	}
	

#pragma mark Power
	
	// C = C.^p
	Matrix& Matrix::elementwisePower(const MathNumber constant) {
		mapElements([&constant](MathNumber& element, Index r, Index c) {
			element = std::pow(element, constant);
		});
		
		return *this;
	}
	
	
	Matrix& Matrix::elementwisePower(const Matrix& matrix) {
		mapElements([&matrix](MathNumber& element, Index r, Index c) {
			element = std::pow(element, matrix(r,c));
		});
		
		return *this;
	}
	
	
#pragma mark - Submatrix operations
	
	Matrix Matrix::range(Index beginRow, Index beginColumn, Index endRow, Index endColumn) {
		
		Index newRows = endRow - beginRow + 1;
		Index newColumns = endColumn - beginColumn + 1;
		
		if (beginRow >= rows() || endRow >= rows() || beginColumn >= columns() || endColumn >= columns()) {
			
			throw new std::invalid_argument("range: requested submatrix does not exist: Begin(" + std::to_string(beginRow) + ", " + std::to_string(beginColumn) + "), End(" + std::to_string(endRow) + ", " + std::to_string(endColumn) + ") on a matrix of size ("+std::to_string(rows())+", "+std::to_string(columns())+")");
			
		} else if (endRow < beginRow || endColumn < beginColumn) {
			
			throw new std::invalid_argument("range: requested submatrix does not exist, upper index is below lower index: Begin(" + std::to_string(beginRow) + ", " + std::to_string(beginColumn) + "), End(" + std::to_string(endRow) + ", " + std::to_string(endColumn) + ")");
			
		} else if (newRows == 0 || newColumns == 0) {
			
			throw new std::invalid_argument("range: cannot build a submatrix with a zero-size dimension.");
			
		}
		
		Matrix submatrix(newRows, newColumns);
		
		submatrix.mapElements([&](MathNumber& newElement, Index r, Index c) {
			// map this matrix's (r+beginRow, c+beginColumn) => (r, c) in the result
			// for all (r, c) in this matrix.
			newElement = (*this)(r+beginRow, c+beginColumn);
		});
		
		return submatrix;
	}
	
	
	void Matrix::insert(const Matrix& source, Index beginRow, Index beginColumn, Index endRow, Index endColumn) {
		
		// VALIDATION
		
		if (source.rows() != (endRow-beginRow+1) || source.columns() != (endColumn-beginColumn+1)) {
			// the source matrix is not the same size as the submatrix that it is to replace
			throw new std::invalid_argument("insert: size mismatch, supplied matrix is ("+std::to_string(source.rows())+", "+std::to_string(source.columns())+"), but requested submatrix is ("+std::to_string(endRow-beginRow+1)+", "+std::to_string(endColumn-beginColumn+1)+").");
			
		} else if (beginRow >= rows() || endRow >= rows() || beginColumn >= columns() || endColumn >= columns()) {
			
			throw new std::invalid_argument("insert: requested submatrix does not exist: Begin(" + std::to_string(beginRow) + ", " + std::to_string(beginColumn) + "), End(" + std::to_string(endRow) + ", " + std::to_string(endColumn) + ") on a matrix of size ("+std::to_string(source.rows())+", "+std::to_string(source.columns())+")");
			
		} else if (endRow < beginRow || endColumn < beginColumn) {
			
			throw new std::invalid_argument("insert: requested submatrix does not exist, upper index is below lower index: Begin(" + std::to_string(beginRow) + ", " + std::to_string(beginColumn) + "), End(" + std::to_string(endRow) + ", " + std::to_string(endColumn) + ")");
		}

		// perform operation
		for (Index column = 0; column < source.columns(); ++column) {
			for (Index row = 0; row < source.rows(); ++row) {
				(*this)(row+beginRow, column+beginColumn) = source(row, column);
			}
		}
	}
	
	
#pragma mark - Other operations
	
	Matrix& Matrix::operator=(const MathNumber constant) {
		mapElements([&constant](MathNumber& element, Index r, Index c) {
			element = constant;
		});
		return *this;
	}
	
	
	Matrix& Matrix::absInPlace() {
		mapElements([](MathNumber& element, Index r, Index c) {
			element = std::abs(element);
		});
		return *this;
	}
	
	
	Matrix Matrix::abs() {
		Matrix result = *this;
		result.absInPlace();
		return result;
	}
	
	
	Matrix& Matrix::transposeInPlace() {
		if(not this->isSquare()) {
			throw new std::logic_error("Cannot transpose a non-square matrix in place");
		}
		
		for (Index row = 0; row < rows(); ++row) {
			for (Index column = 0; column < row; ++column) {
				std::swap( (*this)(row,column), (*this)(column, row) );
			}
		}
		
		return *this;
	}

	
	Matrix Matrix::transpose() {
		
		// new rows = columns
		// new columns = rows
		Matrix result(columns(), rows());
		
		// copy transposed data over
		mapElements([&result](MathNumber& oldElement, Index r, Index c){
			result(c, r) = oldElement;
		});
		
		return result;
	}

	
	Matrix Matrix::inverse() const {
		// check that matrix sizes match
		if (not this->isSquare()) {
			throw new std::logic_error("Cannot invert a non-square matrix.");
		}

		// create two temporary matrices for operation
		Matrix A = *this;
		Matrix B = Matrix::eye(rows());
		
		// solve the linear system A*X=B, filling in X, the inverse
		auto X = Matrix::linearSolve(A, B);
		
		return X;
	}
	
	
# pragma mark linear algebra routines
	
	// solves a linear system A*X = B, returning Matrix X
	Matrix Matrix::linearSolve(Matrix& A, Matrix& B) {
		
		// check that matrix sizes match
		if (A.rows() != B.rows() or not A.isSquare()) {
			throw new std::invalid_argument("linearSolve: size mismatch between matrices.");
		}
		
		// create temporary variables
		long int i, j, k, p;
		MathNumber Amax;
		
		// determine magnitude of entries in A (for singularity check later)
		Amax = Matrix::infNorm(A);
		
		// perform Gaussian elimination to convert A,B to an upper-triangular system
		for (k=0; k<A.rows()-1; k++) {   // loop over diagonals
			
			// find the pivot row p
			p=k;
			
			for (i=k; i<A.rows(); i++) {
				if (std::abs(A._data[k][i]) > std::abs(A._data[k][p])) {
					p=i;
				}
			}
			
			// swap rows in A
			for (j=k; j<A.rows(); j++) {
				std::swap(A._data[j][p], A._data[j][k]);
			}
			
			// swap rows in B
			for (j=0; j<B.columns(); j++) {
				std::swap(B._data[j][p], B._data[j][k]);
			}
			
			// check for nonzero matrix diagonal
			if (std::abs(A._data[k][k]) < STOL*Amax) {
				throw new std::runtime_error("linearSolve: numerically singular matrix.");
			}
			
			// perform elimination using row k
			// store multipliers in column below pivot
			for (i=k+1; i<A.rows(); i++) {
				A._data[k][i] /= A._data[k][k];
			}
			
			for (j=k+1; j<A.rows(); j++) {      // loop over columns of A, to right of pivot
				for (i=k+1; i<A.rows(); i++) {   // update rows in column
					A._data[j][i] -= A._data[k][i]*A._data[j][k];
				}
			}
			
			for (j=0; j<B.columns(); j++) {
				// update entries in B
				for (i=k+1; i<A.rows(); i++) {
					B._data[j][i] -= A._data[k][i]*B._data[j][k];
				}
			}
		} // end Gaussian elimination
		
		// check for singularity at end (only need to check final diagonal entry)
		if (std::abs(A._data[A.rows()-1][A.rows()-1]) < STOL*Amax) {
			throw new std::runtime_error("linearSolve: numerically singular matrix.");
		}
		
		// perform Backward Substitution on result
		Matrix X = backSubstitution(A, B);
		return X;
	}
	
	// solves a linear system A*x = b, returning Vector x
	//    A and b are modified in this operation; x holds the result
	Vector Matrix::linearSolve(Matrix& A, Vector& b) {
		
		// check that matrix sizes match
		if (A.rows() != b.size() or not A.isSquare()) {
			throw new std::invalid_argument("linearSolve: size mismatch between matrix and vector.");
		}
		
		// create temporary variables
		Index i, j, k, p;
		MathNumber Amax;
		
		// determine magnitude of entries in A (for singularity
		// check later)
		Amax = infNorm(A);
		
		// perform Gaussian elimination to convert A,B to an
		// upper-triangular system
		for (k = 0; k < A.rows() - 1; k++) {	// loop over diagonals
			
			// find the pivot row p
			p = k;
			
			for (i = k; i < A.rows(); i++) {
				if (std::abs(A._data[k][i]) > std::abs(A._data[k][p])) {
					p = i;
				}
			}
			
			// swap rows in A
			for (j = k; j < A.rows(); j++) {
				std::swap(A._data[j][p], A._data[j][k]);
			}
			
			// swap rows in b
			std::swap(b[p], b[k]);
			
			// check for nonzero matrix diagonal
			if (std::abs(A._data[k][k]) < STOL * Amax) {
				throw new std::runtime_error("linearSolve: numerically singular matrix.");
			}
			
			// perform elimination using row k
			for (i = k + 1; i < A.rows(); i++) {
				// Store multipliers in column below pivot
				A._data[k][i] /= A._data[k][k];
			}
			
			// Loop over columns of A, to right of pivot
			for (j = k + 1; j < A.rows(); j++) {
				// update rows in column
				for (i = k + 1; i < A.rows(); i++) {
					A._data[j][i] -= A._data[k][i] * A._data[j][k];
				}
			}
			// update entries in b
			for (i = k + 1; i < A.rows(); i++) {
				b[i] -= A._data[k][i] * b[k];
			}
		}
		
		// check for singularity at end (only need to check final
		// diagonal entry)
		if (std::abs(A._data[A.rows() - 1][A.rows() - 1]) < STOL * Amax) {
			throw new std::runtime_error("linearSolve: numerically singular matrix.");
		}
		
		// perform Backward Substitution on result
		Vector x = backSubstitution(A, b);
		return x;
		
	} // linearSolve
	
	
	// backward substitution on the linear system U*X = B, filling in a Matrix X
	Matrix Matrix::backSubstitution(const Matrix& U, const Matrix& B) {
		
		auto X = Matrix(U.rows(), B.columns());
		
		// check that matrix sizes match
		if (U.rows() != B.rows() || not U.isSquare()) {
			throw new std::invalid_argument("backSubstitution: size mismatch between matrix and vector.");
		}
		
		// copy B into X
		X = B;
		
		// analyze matrix for magnitude
		MathNumber Umax = Matrix::infNorm(U);
		
		// perform column-oriented Backward Substitution algorithm
		for (Index j = U.rows() - 1; j-- > 0;) {
			// check for nonzero matrix diagonal
			if (std::abs(U._data[j][j]) < STOL*Umax) {
				throw new std::runtime_error("backSubstitution: numerically singular matrix.");
			}
			
			// solve for this row of solution
			for (Index k=0; k < X.columns(); k++) {
				X(j,k) /= U(j,j);
			}
			
			// update all remaining rhs
			for (Index k=0; k<X.columns(); k++) {
				for (Index i=0; i<j; i++) {
					X._data[k][i] -= U._data[j][i]*X._data[k][j];
				}
			}
		} // end Back sub.
		
		// return success
		return X;
	}
	
	// backward substitution on U*x = b, returning Vector x
	Vector Matrix::backSubstitution(const Matrix& U, const Vector& b) {
		
		auto x = Vector(U.rows());
		
		// check that matrix sizes match
		if (U.rows() != b.size() || U.rows() != U.columns()) {
			throw new std::invalid_argument("backSubstitution: size mismatch between matrix and vector.");
		}
		
		// copy b into x
		x = b;
		
		// analyze matrix for magnitude
		MathNumber Umax = Matrix::infNorm(U);
		
		// perform column-oriented Backward Substitution algorithm
		for (long int j=U.rows()-1; j>=0; j--) {
			// check for nonzero matrix diagonal
			if (std::abs(U._data[j][j]) < STOL * Umax) {
				throw new std::runtime_error("backSubstitution: numerically singular matrix.");
			}
			
			// solve for this row of solution
			x[j] /= U._data[j][j];
			
			// update all remaining rhs
			for (long int i=0; i<j; i++) {
				x[i] -= U._data[j][i]*x[j];
			}
		}
		
		// return success
		return x;
	}
	
	
	// forward substitution on the linear system L*X = B, filling in Matrix X
	//    L and B remain unchanged in this operation; X holds the result
	//    B and X may have multiple columns
	Matrix Matrix::forwardSubstitution(const Matrix& L, const Matrix& B) {
		
		auto X = Matrix(L.rows(), B.columns());
		
		// check that matrix sizes match
		if (L.rows() != B.rows() or not L.isSquare()) {
			throw new std::invalid_argument("forwardSubstitution: size mismatch between matrix and vector.");
		}
		
		// copy B into X
		X = B;
		
		// analyze matrix magnitude
		MathNumber Lmax = Matrix::infNorm(L);
		
		// perform column-oriented Forwards Substitution algorithm
		for (long int j=0; j<L.rows(); j++) {
			// check for nonzero matrix diagonal
			if (std::abs(L._data[j][j]) < STOL*Lmax) {
				throw new std::runtime_error("backSubstitution: numerically singular matrix.");
			}
			
			// solve for this row of solution
			for (long int k=0; k<X.columns(); k++) {
				X._data[k][j] /= L._data[j][j];
			}
			
			// update all remaining rhs
			for (long int k=0; k<X.columns(); k++) {
				for (long int i=j+1; i<L.rows(); i++) {
					X._data[k][i] -= L._data[j][i]*X._data[k][j];
				}
			}
		} // end Column-oriented forward sub.
		
		// return success
		return X;
	}
	
	
	// forward substitution on L*x = b, filling in a resulting Vector x
	//    L and b remain unchanged in this operation; x holds the result
	Vector Matrix::forwardSubstitution(const Matrix& L, const Vector& b) {
		
		auto x = Vector(L.rows());
		
		// check that matrix sizes match
		if (L.rows() != b.size() || not L.isSquare()) {
			throw new std::invalid_argument("forwardSubstitution: size mismatch between matrix and vector.");
		}
		
		x = b;
		
		// analyze matrix for magnitude
		MathNumber Lmax = Matrix::infNorm(L);
		
		// perform column-oriented Forwards Substitution algorithm
		for (Index j=0; j<L.rows(); j++) {
			
			// check for nonzero matrix diagonal
			if (std::abs(L._data[j][j]) < STOL*Lmax) {
				throw new std::runtime_error("forwardSubstitution: numerically singular matrix.");
			}
			
			// solve for this row of solution
			x[j] /= L._data[j][j];
			
			// update all remaining rhs
			for (long int i=j+1; i<L.rows(); i++) {
				x[i] -= L._data[j][i]*x[j];
			}
		}
		
		// return success
		return x;
	}
	
	
	// Matrix Dot Product
	MathNumber dotProduct(const Matrix& A, const Matrix& B) {
		
		if (A.dimensions() != B.dimensions()) {
			throw new std::invalid_argument("dot product: incompatible matrix dimensions (must be same)");
		}
		
		MathNumber sum = 0.0;
		
		for(Index column = 0; column < A.columns(); ++column) {
			for(Index row = 0; row < A.rows(); ++row) {
				// sum = A(r,c) * B(r,c) + sum
				sum = std::fma(A(row,column), B(row, column), sum);
			}
		}
		
		return sum;
	}

	
# pragma mark Scalar output methods

	// minimum entry in the matrix
	MathNumber Matrix::min() const {
		MathNumber mn = _data[0][0];
		for (Index column = 0; column < columns(); ++column) {
			for (Index row = 0; row < rows(); ++row) {
				mn = std::min(mn, (*this)(row,column));
			}
		}
		return mn;
	}

	
	// maximum entry in the matrix
	MathNumber Matrix::max() const {
		MathNumber mx = _data[0][0];
		for (Index column = 0; column < columns(); ++column) {
			for (Index row = 0; row < rows(); ++row) {
				mx = std::max(mx, (*this)(row,column));
			}
		}
		return mx;
	}
	
	
	// matrix Frobenius norm (column/row vector 2-norm)
	MathNumber Matrix::norm(const Matrix& A) {
		MathNumber sum=0.0;
		for (Index j=0; j<A.columns(); j++)
			for (Index i=0; i<A.rows(); i++)
				sum += A(i,j)*A(i,j);
		return sqrt(sum);
	}
	
	
	// matrix infinity norm (column vector infinity norm, row vector one norm)
	MathNumber Matrix::infNorm(const Matrix& A) {
		MathNumber mx=0.0;
		for (Index i=0; i<A.rows(); i++) {
			MathNumber sum=0.0;
			for (Index j=0; j<A.columns(); j++)
				sum += std::abs(A(i,j));
			mx = std::max(mx,sum);
		}
		return mx;
	}
	
	
	// matrix one norm (column vector one norm, row vector infinity norm)
	MathNumber Matrix::oneNorm(const Matrix& A) {
		MathNumber mx=0.0;
		for (Index j=0; j<A.columns(); j++) {
			MathNumber sum=0.0;
			for (Index i=0; i<A.rows(); i++)
				sum += std::abs(A(i,j));
			mx = std::max(mx,sum);
		}
		return mx;
	}
	
	
# pragma mark Infix operators
	
	// Result = Matrix * Matrix
	Matrix operator*(const Matrix& matrix1, const Matrix& matrix2) {
		
		// if A has (X rows, Y cols) and B has (Y rows, Z cols),
		// they can be multiplied, producing a matrix of size
		// (X rows, Z cols)
		
		if (matrix1.columns() != matrix2.rows()) {
			throw new std::invalid_argument("matrix-by-matrix multiplication: incompatible matrix dimensions (inner dimensions must be same, A*B requires A's columns to be the same as B's rows)");
		}
		
		// perform operation
		Matrix result = Matrix(matrix1.rows(), matrix2.columns());
		
		for (Index k = 0; k< matrix2.columns(); ++k) {
			for (Index j = 0; j < matrix1.columns(); ++j) {
				for (Index i = 0; i < matrix1.rows(); ++i) {
					result(i,k) += matrix1(i,j)*matrix2(j,k);
				}
			}
		}
		
		return result;
	}
	
	
	// Result = Matrix * Vector
	Vector operator*(const Matrix& matrix, const Vector& vector) {
		
		Vector result(matrix.rows());
		
		if (matrix.columns() != vector.size()) {
			throw new std::invalid_argument("matrix-vector product: incompatible matrix/vector sizes in A*v");
		}
		
		for(Index column = 0; column < matrix.columns(); ++column) {
			for(Index row = 0; row < matrix.rows(); ++row) {
				result += matrix(row,column) * vector[column];
			}
		}
		
		return result;
	}
	
	
	bool operator==(const Matrix& lhs, const Matrix& rhs) {
		
		if (lhs.dimensions() != rhs.dimensions()) {
			return false;
		}
		
		// assume equality, disprove with single unequal element
		
		for (Index column = 0; column < lhs.columns(); column++) {
			for (Index row = 0; row < lhs.rows(); row++) {
				MathNumber left = lhs(row,column);
				MathNumber right = rhs(row,column);
				bool isEqual = withinTolerance(left,right);
				if(not isEqual) { return false; }
			}
		}
		
		return true;
	}
	
	
	// Result = Matrix - Constant
	Matrix operator-(const Matrix& matrix, const MathNumber constant) {
		Matrix result = matrix;
		result -= constant;
		return result;
	}
	
	// Result = Constant - Matrix
	Matrix operator-(const MathNumber constant, const Matrix& matrix) {
		Matrix result = matrix;
		result.mapElements([&constant](MathNumber& element, Index r, Index c){
			element = constant - element;
		});
		return result;
	}
	
	// Result = Matrix - Matrix
	Matrix operator-(const Matrix& matrix1, const Matrix& matrix2) {
		Matrix result = matrix1;
		result -= matrix2;
		return result;
	}
	
	
	bool withinTolerance(const MathNumber a, const MathNumber b, const double precision) { // precision default = 1e-6
		int aExp, bExp;
		double aFrac, bFrac;
		
		aFrac = std::frexp(a, &aExp);
		bFrac = std::frexp(b, &bExp);
		
		// are the significands within the difference tolerance,
		// and are exponents the same?
		return (std::abs(aFrac - bFrac) < precision) && (aExp == bExp);
	}
	
	
} // namespace PH
