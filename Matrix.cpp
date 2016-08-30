
#include "Matrix.h"

using namespace std;

// singularity tolerance
#define STOL 1.e-15

namespace PH {

	// general constructor (initializes values to 0.0)
	Matrix::Matrix(Index r, Index c) {
		this->resize(r,c);
	}
	
	
	template<int r, int c>
	Matrix Matrix::fromArray(Raw1DArray source) {
		if (r * c != source.size()) {
			throw new std::invalid_argument("Matrix::fromArray: the given array could not properly fit into the requested matrix dimensions.");
		}
		
		Matrix result = Matrix(r,c);
		
		result.mapElements([source](MathNumber& element, Index row, Index column) {
			element = source[row * column + column];
		});
		
		return result;
	}

	// constructor that copies input data (2D vector)
	Matrix::Matrix(Raw2DArray source) {
		_columnCount = source.size();
		_rowCount = source[0].size();
		
		for (Index column = 0; column < _columnCount; column++) {
			if (source[column].size() != _rowCount) {
				throw new std::invalid_argument("Matrix Raw2DArray constructor: rows in the 2D array must have identical length.");
			}
			
			for (Index row = 0; row < _rowCount; row++) {
				(*this)(row, column) = source[column][row];
			}
		}
	}
	
	
	// create a new n by n identity matrix
	Matrix Matrix::eye(const Index size) {
		Matrix I(size, size);
		for (Index i = 0; i < size; i++) {
			I(i,i) = 1.0;
		}
		return I;
	}
	
	
	// create a new matrix of linearly spaced data
	Matrix Matrix::linSpace(MathNumber a, MathNumber b, Index m, Index n) {
		Matrix C(m,n);
		MathNumber h = (b-a)/(m*n-1);
		Index idx=0;
		for (Index j=0; j<n; j++)
			for (Index i=0; i<m; i++)
				C._data[j][i] = a + (idx++)*h;
		return C;
	}
	
	// create a new column-vector matrix of logarithmically spaced data
	Matrix Matrix::logSpace(MathNumber a, MathNumber b, Index m, Index n) {
		Matrix C(m,n);
		MathNumber h = (b-a)/(m*n-1);
		Index idx=0;
		for (Index j=0; j<n; j++)
			for (Index i=0; i<m; i++)
				C._data[j][i] = pow(10.0, a + (idx++)*h);
		return C;
	}
	
	
	// Raw2DArray -> conversion by assignment constructor
	Matrix& Matrix::operator=(const Raw2DArray& source) {
		
	}

	// copy constructor
	Matrix::Matrix(const Matrix& A) {
		this->resize(A.rows(), A.columns());

		for (Index column = 0; column < _columnCount; column++) {
			_data[column] = A._data[column];
		}
	}

	// C = A
	Matrix& Matrix::operator=(const Matrix& A) {
		this->resize(A.rows(), A.columns());

		for (Index column = 0; column < _columnCount; column++) {
			_data[column] = A._data[column];
		}
		
		return *this;
	}

	Matrix::Matrix(Matrix&& A) {
		*this = std::move(A);
	}

	/// Move assignment operator
	Matrix& Matrix::operator=(Matrix&& A) {
		if(this != &A) {
			_rowCount = _columnCount = 0;
			_data.resize(0);
			
			_rowCount = A.rows();
			_columnCount = A.columns();
			
			_data = std::move(A._data);
			A._rowCount = 0;
			A._columnCount = 0;
			A._data = Raw2DArray();
		}
		return *this;
	}
	
	
	// Serializable interface
	void Matrix::serialize(std::ostream& output) const {
		output << str();
	}
	
	Matrix& Matrix::deserialize(std::istream& input) {
		// determine matrix size
		Index newRowCount = 0, newColumnCount = 0;
		std::string line;
		
		while (getline(input, line)) {
			std::istringstream iss(line);
			float columnsOnRow;
			Index n = 0;
			while (iss >> columnsOnRow) {
				n++;
			}
			
			if((n > 0) && (newRowCount == 0)) {
				// first row, set _columnCount
				newColumnCount = n;
			}
			
			if ((n > 0) && (n != _columnCount)) {
				throw new std::runtime_error("matrix deserialize: not all rows had same column count.");
			}
			
			if (n > 0) {
				++newRowCount;
			}
		}
		
		input.seekg(0, ios::beg);
		
		// create matrix of desired size
		Matrix result = Matrix(newRowCount, newColumnCount);
		
		// load matrix based on data from file
		for (Index i = 0; i < _rows; i++) {
			
			getline(input, line);
			std::istringstream iss(line);
			
			for (Index column = 0; column < newColumnCount; ++column) {
				iss >> result(row, column);
			}
		}
		
		return result;

	}


	// column accessor routines
	MathVector& Matrix::column(Index i) {
	  return _data[i];
	}
	
	// row accessor (copy) routine
	MathVector Matrix::row(Index row) {
		MathVector rowVector = MathVector(_columnCount);
		
		for (Index column = 0; column < _columnCount; column++) {
			rowVector[column] = _data[column][row];
		}
		
		return rowVector;
	}

	// Matlab/Fortran Matrix accessors (row, column)
	MathNumber& Matrix::operator()(Index row, Index column) {
		return _data[column][row];
	}
	MathNumber Matrix::operator()(Index row, Index column) const {
		return _data[column][row];
	}
	MathNumber& Matrix::operator()(Index linearIndex) {
		return _data[linearIndex / _rows][linearIndex % _rows];
	}
	MathNumber Matrix::operator()(Index linearIndex) const {
		return _data[linearIndex / _rows][linearIndex % _rows];
	}
	
	// build a string representation and return it
	std::string Matrix::str() const {
		std::stringstream ss();
		
		for (Index row = 0; row < _rows; row++) {
			for (Index column = 0; column < _columnCount; column++) {
				ss << _data[column][row];
			}
			ss << "\n";
		}
		
		return ss.str();
	}

	// streaming output routine
	ostream& operator<<(ostream& os, const Matrix& A) {
		os << this->str();
	}
	
	# pragma mark In-Place Member Operations

	// A = (a*A) + (b*B)
	Matrix& Matrix::linearSumInPlace(MathNumber matrix1Constant, MathNumber matrix2Constant, const Matrix& matrix2) {
		// check that array sizes match
		if (this.dimensions() == matrix2.dimensions()) {
			throw new std::invalid_argument("matrix-by-matrix linear sum (a*A+b*B): incompatible matrix dimensions (must be same)");
		}
		
		mapElements([](MathNumber& matrix1Element, Index row, Index column) {
			matrix1Element *= matrix1Constant;
			// fused multiply-add operation ($0 * $1) + $2
			// usually handled in one operation on a modern processor.
			matrix1Element = std::fma(matrix2Constant, matrix2(row, column), matrix1Element);
		});
		
		return *this;
	}

	
	// ADDITION
	
	Matrix& Matrix::operator+=(const MathNumber constant) {
		mapElements([constant](MathNumber& element, Index row, Index column) {
			element += constant;
		});
		
		return *this;
	}
	
	Matrix& Matrix::operator+=(const Matrix& matrix) {
		linearSumInPlace(1.0, 1.0, matrix);
		
		return *this;
	}
	
	
	// SUBTRACTION
	
	Matrix& Matrix::operator-=(const MathNumber constant) {
		this += -constant;
		return *this;
	}
	
	Matrix& Matrix::operator-=(const Matrix& matrix) {
		linearSumInPlace(1.0, -1.0, matrix);
		
		return *this;
	}
	
	
	// MULTIPLICATION
	
	Matrix& Matrix::elementwiseMultiply(const Matrix& other) {
		// check that array sizes match
		if (this.dimensions() != other.dimensions()) {
			throw new std::invalid_argument("elementwise matrix multiplication (A.*B): incompatible matrix dimensions (must be same)");
		}

		// perform operation
		mapElements([other](MathNumber& element, Index row, Index column) {
			element *= other(row, column);
		});

		return *this;
	}
	
	Matrix& Matrix::operator*=(const MathNumber constant) {
		// perform operation
		mapElements([constant](MathNumber& element, Index row, Index column) {
			element *= constant;
		});

		return *this;
	}
	
	
	// DIVISION
	
	Matrix& Matrix::elementwiseDivide(const Matrix& matrix) {
		// check that array sizes match
		if (this.dimensions() != matrix.dimensions()) {
			throw new std::invalid_argument("elementwise matrix division (A./B): incompatible matrix dimensions (must be same)");
		}
		
		mapElements([matrix](MathNumber& element, Index row, Index column) {
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
		mapElements([constant](MathNumber& element, Index row, Index column) {
			element /= constant;
		});
		
		return *this;
	}
	
	
	
	//     is,ie,js,je negative  =>  offset from end of dimension (-1 == end)
	int Matrix::insert(const Matrix& source, Index beginRow, Index beginColumn, Index endRow, Index endColumn) {
		
		// update indices if negative
		
		beginRow = (beginRow < 0) ? beginRow + _rows : beginRow;
		endRow = (endRow < 0) ? endRow + _rows : endRow;
		beginColumn = (beginColumn < 0) ? beginColumn + _columnCount : beginColumn;
		endColumn = (endColumn < 0) ? endColumn + _columnCount : endColumn;

		// VALIDATION
		
		if (source.rows() != (endRow-beginRow+1) || source.columns() != (endColumn-beginColumn+1)) {
			// the source matrix is not the same size as the submatrix that it is to replace
			throw new std::invalid_argument("insert: size mismatch, supplied matrix is ("+std::to_string(source.rows())+", "+std::to_string(source.columns())+"), but requested submatrix is ("+std::to_string(endRow-beginRow+1)+", "+std::to_string(endColumn-beginColumn+1)+").");
			
		} else if (beginRow < 0 || beginRow >= rows() || endRow < 0 || endRow >= rows() || beginColumn < 0 || beginColumn >= columns() || endColumn < 0 || endColumn >= columns()) {
			
			throw new std::invalid_argument("insert: requested submatrix does not exist: Begin(" + std::to_string(beginRow) + ":" + std::to_string(beginColumn) + "), End(" + std::to_string(endRow) + ":" + std::to_string(endColumn) + ") on a matrix of size ("+std::to_string(source.rows())+", "+std::to_string(source.columns())+")");
			
		} else if (endRow < beginRow || endColumn < beginColumn) {
			
			throw new std::invalid_argument("insert: requested submatrix does not exist, upper index is below lower index: Begin(" + std::to_string(beginRow) + ":" + std::to_string(beginColumn) + "), End(" + std::to_string(endRow) + ":" + std::to_string(endColumn) + ")");
		}

		// perform operation
		for (Index column = 0; j < A.columns(); ++column) {
			for (Index row = 0; row < A.rows(); ++row) {
				(*this)(row+beginRow, column+beginColumn) = source(row, column);
			}
		}
	}
	
	
	// constant fill operator
	Matrix& Matrix::operator=(const MathNumber constant) {
		mapElements([constant](MathNumber& element, Index r, Index c) {
			element = constant;
		});
	}

	
	// C = C.^p
	Matrix& Matrix::elementwisePower(const MathNumber constant) {
		mapElements([constant](MathNumber& element, Index r, Index c) {
			element = std::pow(element, constant);
		});
		
		return *this;
	}
	
	
	Matrix& Matrix::elementwisePower(const Matrix& matrix) {
		mapElements([other](MathNumber& element, Index r, Index c) {
			element = std::pow(element, other(r,c));
		});
		
		return *this;
	}
	
	
	Matrix Matrix::abs() {
		Matrix result = this;
		result.absInPlace();
		return result;
	}
	
	
	Matrix& Matrix::absInPlace() {
		mapElements([](MathNumber& element, Index r, Index c) {
			element = std::abs(element);
		});
		return *this;
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

	
	// Cij = Cji
	Matrix Matrix::transpose() {
		
		// new rows = columns
		// new columns = rows
		Matrix result(columns(), rows());
		
		// copy transposed data over
		mapElements([result](MathNumber& oldElement, Index r, Index c){
			result(c, r) = oldElement;
		});
		
		return result;
	}

	// computes the inverse of a nonsingular matrix 
	Matrix Matrix::inverse() const {
		// check that matrix sizes match
		if (this.isSquare()) {
			throw new std::logic_error("Cannot invert a non-square matrix.");
		}

		// create two temporary matrices for operation
		Matrix A = *this;
		Matrix X = Matrix(rows(), columns());
		Matrix B = Matrix::eye(rows());
		
		// solve the linear system A*X=B, filling in X, the inverse
		Matrix::linearSolve(A, X, B);
		
		// return success
		return 0;
	}


	///// Derived matrix creation operations /////

	// C = A^T
	Matrix Matrix::T() {
	  Matrix C(*this);
	  C.Trans();
	  return C;
	}

	// submatrix extraction routine (creates a matrix from a portion of an existing Mat)
	//     is,ie,js,je negative  =>  offset from end of dimension (-1 == end)
	Matrix Matrix::Extract(long int is, long int ie, long int js, long int je) {

	  // update is,ie,js,je if any are negative
	  is = (is < 0) ? is+_rows : is;
	  ie = (ie < 0) ? ie+_rows : ie;
	  js = (js < 0) ? js+_columnCount : js;
	  je = (je < 0) ? je+_columnCount : je;

	  // check that requested submatrix exists
	  if (is < 0 || is >= _rows) {
		cerr << "Matrix::Extract error, requested submatrix does not exist\n";
		cerr << "  illegal is = " << is << " (matrix has " << _rows << " rows)\n";
	  }
	  if (ie < 0 || ie >= _rows) {
		cerr << "Matrix::Extract error, requested submatrix does not exist\n";
		cerr << "  illegal ie = " << ie << " (matrix has " << _rows << " rows)\n";
	  }
	  if (js < 0 || js >= _columnCount) {
		cerr << "Matrix::Extract error, requested submatrix does not exist\n";
		cerr << "  illegal js = " << js << " (matrix has " << _columnCount << " columns)\n";
	  }
	  if (je < 0 || je >= _columnCount) {
		cerr << "Matrix::Extract error, requested submatrix does not exist\n";
		cerr << "  illegal je = " << je << " (matrix has " << _columnCount << " columns)\n";
	  }
	  if (ie < is || je < js) {
		cerr << "Matrix::Extract error, requested submatrix does not exist\n";
		cerr << "  upper index is below lower index: is = " << is << ", ie = " 
		 << ie << ", js = " << js << ", je = " << je << endl;
	  }

	  // create new matrix of desired size
	  Matrix C(ie-is+1, je-js+1);

	  // copy requested data
	  for (Index j=js; j<=je; j++) 
		for (Index i=is; i<=ie; i++) 
		  C._data[j-js][i-is] = _data[j][i];

	  // return object
	  return C;
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
		Amax = InfNorm(A);
		
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
	static Matrix Matrix::backSubstitution(const Matrix& U, const Matrix& B) {
		
		auto X = Matrix(U.rows(), B.columns());
		
		// check that matrix sizes match
		if (U.rows() != B.rows() || not U.isSquare()) {
			throw new std::invalid_argument("backSubstitution: size mismatch between matrix and vector.");
		}
		
		// copy B into X
		X = B;
		
		// analyze matrix for magnitude
		MathNumber Umax = InfNorm(U);
		
		// perform column-oriented Backward Substitution algorithm
		for (Index j=U.rows()-1; j>=0; j--) {
			// check for nonzero matrix diagonal
			if (std::abs(U._data[j][j]) < STOL*Umax) {
				cerr << "BackSubstitution error: numerically singular matrix!\n";
				return emptyResult;
			}
			
			// solve for this row of solution
			for (Index k=0; k<X.columns(); k++) {
				X._data[k][j] /= U._data[j][j];
			}
			
			// update all remaining rhs
			for (Index k=0; k<X.columns(); k++) {
				for (Index i=0; i<j; i++) {
					X._data[k][i] -= U._data[j][i]*X._data[k][j];
				}
			}
		} // end Back sub.
		
		// return success
		return std::optional<Matrix>(X);
	}
	
	// backward substitution on U*x = b, returning Vector x
	Vector backSubstitution(const Matrix& U, const Vector& b) {
		
		auto x = Vector(U.rows());
		
		// check that matrix sizes match
		if (U.rows() != b.size() || U.rows() != U.columns()) {
			throw new std::invalid_argument("backSubstitution: size mismatch between matrix and vector.");
		}
		
		// copy b into x
		x = b;
		
		// analyze matrix for magnitude
		MathNumber Umax = InfNorm(U);
		
		// perform column-oriented Backward Substitution algorithm
		for (long int j=U.rows()-1; j>=0; j--) {
			// check for nonzero matrix diagonal
			if (std::abs(U._data[j][j]) < STOL * Umax) {
				cerr << "BackSubstitution error: numerically singular matrix!\n";
				return emptyResult;
			}
			
			// solve for this row of solution
			x[j] /= U._data[j][j];
			
			// update all remaining rhs
			for (long int i=0; i<j; i++) {
				x[i] -= U._data[j][i]*x[j];
			}
		}
		
		// return success
		return std::optional<Vector>(x);
	}
	
	
	// forward substitution on the linear system L*X = B, filling in Matrix X
	//    L and B remain unchanged in this operation; X holds the result
	//    B and X may have multiple columns
	Matrix forwardSubstitution(const Matrix& L, const Matrix& B) const {
		
		auto X = Matrix(L.rows(), B.columns());
		
		// check that matrix sizes match
		if (L.rows() != B.rows() or not L.isSquare()) {
			throw new std::invalid_argument("forwardSubstitution: size mismatch between matrix and vector.");
		}
		
		// copy B into X
		X = B;
		
		// analyze matrix magnitude
		MathNumber Lmax = InfNorm(L);
		
		// perform column-oriented Forwards Substitution algorithm
		for (long int j=0; j<L.rows(); j++) {
			// check for nonzero matrix diagonal
			if (std::abs(L._data[j][j]) < STOL*Lmax) {
				cerr << "ForwardSubstitution error: singular matrix!\n";
				return emptyResult;
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
		return std::optional<Matrix>(X);
	}
	
	
	// forward substitution on L*x = b, filling in a resulting Vector x
	//    L and b remain unchanged in this operation; x holds the result
	Vector forwardSubstitution(const Matrix& L, const Vector& b) const {
		
		auto x = Vector(L.rows());
		
		// check that matrix sizes match
		if (L.rows() != b.size() || not L.isSquare()) {
			throw new std::invalid_argument("forwardSubstitution: size mismatch between matrix and vector.");
		}
		
		x = b;
		
		// analyze matrix for magnitude
		MathNumber Lmax = infNorm(L);
		
		// perform column-oriented Forwards Substitution algorithm
		for (long int j=0; j<L.rows(); j++) {
			
			// check for nonzero matrix diagonal
			if (std::abs(L._data[j][j]) < STOL*Lmax) {
				cerr << "ForwardSubstitution error: singular matrix!\n";
				return emptyResult;
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

	
# pragma mark Scalar Output methods

	// minimum entry in the matrix
	MathNumber Matrix::min() const {
		MathNumber mn = _data[0][0];
		for (Index column = 0; column < columns(); ++column) {
			for (Index row = 0; i < rows(); i++) {
				mn = std::min(mn, (*this)(row,column));
			}
		}
		return mn;
	}

	// maximum entry in the matrix
	MathNumber Matrix::max() const {
		MathNumber mx = _data[0][0];
		for (Index column = 0; column < columns(); ++column) {
			for (Index row = 0; i < rows(); i++) {
				mx = std::max(mx, (*this)(row,column));
			}
		}
		return mx;
	}
	
	// matrix Frobenius norm (column/row vector 2-norm)
	static MathNumber Matrix::norm(const Matrix& A) {
		MathNumber sum=0.0;
		for (Index j=0; j<A.columns(); j++)
			for (Index i=0; i<A.rows(); i++)
				sum += A(i,j)*A(i,j);
		return sqrt(sum);
	}
	
	// matrix infinity norm (column vector infinity norm, row vector one norm)
	static MathNumber Matrix::infNorm(const Matrix& A) {
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
	static MathNumber Matrix::oneNorm(const Matrix& A) {
		MathNumber mx=0.0;
		for (Index j=0; j<A.columns(); j++) {
			MathNumber sum=0.0;
			for (Index i=0; i<A.rows(); i++)
				sum += std::abs(A(i,j));
			mx = std::max(mx,sum);
		}
		return mx;
	}
	
	# pragma mark Functional Operators
	
	/// Result = Matrix * Matrix
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
					(*this)(i,k) += A(i,j)*X(j,k);
				}
			}
		}
		
		return result;
	}
	
	bool operator==(const Matrix& lhs, const Matrix& rhs) {
		
		if (lhs.dimensions() != rhs.dimensions()) {
			return false;
		}
		
		// assume equality, disprove with single unequal element
		bool isEqual = true;
		
		lhs.mapElements([](){
		
		});
		
		for (Index j=0; j<_columnCount; j++)
			for (Index i=0; i<_rows; i++)
		  equal &= (A._data[j][i] == _data[j][i]);
		
		return isEqual;
	}
	
} // namespace PH
