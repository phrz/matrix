
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
	static Matrix& Matrix::fromArray(Raw1DArray source) {
		if (r * c != source.size()) {
			throw new std::invalid_argument("Matrix::fromArray: the given array could not properly fit into the requested matrix dimensions.");
		}
		
		this->resize(r,c);
		
		Index sourceIndex = 0;
		
		for (Index column = 0; column < _columnCount; column++) {
			for (Index row = 0; row < _rows; row++) {
				_data[column][row] = source[sourceIndex++];
			}
		}
	}

	// constructor that copies input data (2D vector)
	Matrix::Matrix(Raw2DArray source) {
		_columnCount = source.size();
		_rows = source[0].size();
		
		for (Index column = 0; column < _columnCount; column++) {
			if (source[column].size() != _rows) {
				throw new std::invalid_argument("Matrix Raw2DArray constructor: rows in the 2D array must have identical length.");
			}
			
			for (Index row = 0; row < _rows; row++) {
				_data[column][row] = source[column][row];
			}
		}
	}
	
	// Raw2DArray -> conversion by assignment constructor
	Matrix& Matrix::operator=(const Raw2DArray& source) {
		
	}

	// copy constructor
	Matrix::Matrix(const Matrix& A) {
		_rows = A.rows();
		_columnCount = A.columns();
		_data.resize(_columnCount);

		for (Index column = 0; column < _columnCount; column++) {
			_data[column].resize(_rows);
			_data[column] = A._data[column];
		}
	}

	// C = A
	Matrix& Matrix::operator=(const Matrix& A) {
		_rows = A.rows();
		_columnCount = A.columns();
		_data.resize(_columnCount);

		for (Index column = 0; column < _columnCount; column++) {
			_data[column].resize(_rows);
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
			_rows = _columnCount = 0;
			_data.resize(0);
			
			_rows = A.rows();
			_columnCount = A.columns();
			
			_data = std::move(A._data);
			A._rows = 0;
			A._columnCount = 0;
			A._data.clear();
		}
		return *this;
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

	// write myself to a file
	int Matrix::Write(const char *outfile) const {

	  // return with failure if 'outfile' is empty
	  if (strlen(outfile) < 1) {
		cerr << "Matrix::Write error, empty outfile\n";
		return 1;
	  }

	  // open output file
	  FILE *fptr = NULL;
	  fptr = fopen(outfile, "w");
	  if (fptr == NULL) {
		cerr << "Matrix::Write error, unable to open " << outfile << " for writing\n";
		return 1;
	  }

	  // print data to file
	  for (Index row = 0; row < _rows; row++) {
		for (Index column = 0; column < _columnCount; column++)
		  fprintf(fptr, "  %.16g", _data[column][row]);
		fprintf(fptr, "\n");
	  }

	  // close output file and return
	  fclose(fptr);
	  return 0;
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
	
	Matrix& Matrix::elementwiseDivide(const Matrix& other) {

	  // check that array sizes match
	  if (A._rows != _rows || A._columnCount != _columnCount) {
		cerr << "Matrix::Divide error, matrix size mismatch\n";
		cerr << "  Matrix 1 is " << _rows << " x " << _columnCount 
		 << ",  Matrix 2 is " << A._rows << " x " << A._columnCount << endl;
		return 1;
	  }
	  
	  // perform operation where A is nonzero, otherwise return with error
	  for (Index j=0; j<_columnCount; j++)
		for (Index i=0; i<_rows; i++) {
		  if (A(i,j) == 0.0)
		return 1;
		  _data[j][i] /= A._data[j][i];
		}

	  // return success
	  return 0;
	}

	//   C(is:ie,js:je) = A  (inserts values from A into a submatrix of C)
	//     is,ie,js,je negative  =>  offset from end of dimension (-1 == end)
	int Matrix::Insert(const Matrix& A, long int is, long int ie, long int js, long int je) {

	  // update is,ie,js,je if any are negative
	  is = (is < 0) ? is+_rows : is;
	  ie = (ie < 0) ? ie+_rows : ie;
	  js = (js < 0) ? js+_columnCount : js;
	  je = (je < 0) ? je+_columnCount : je;

	  // check that array sizes match
	  if (A._rows != (ie-is+1) || A._columnCount != (je-js+1)) {
		cerr << "Matrix::Insert error, matrix size mismatch\n";
		cerr << "  supplied Matrix is " << A._rows << " x " << A._columnCount 
		 << ", but requested submatrix is " << ie-is+1 << " x " 
		 << je-js+1 << endl;
		return 1;
	  }
	  // check for valid submatrix
	  if (is < 0 || is >= _rows) {
		cerr << "Matrix::Insert error, requested submatrix does not exist\n";
		cerr << "  illegal is = " << is << " (matrix has " << _rows << " rows)\n";
		return 1;
	  }
	  if (ie < 0 || ie >= _rows) {
		cerr << "Matrix::Insert error, requested submatrix does not exist\n";
		cerr << "  illegal ie = " << ie << " (matrix has " << _rows << " rows)\n";
		return 1;
	  }
	  if (js < 0 || js >= _columnCount) {
		cerr << "Matrix::Insert error, requested submatrix does not exist\n";
		cerr << "  illegal js = " << js << " (matrix has " << _columnCount << " columns)\n";
		return 1;
	  }
	  if (je < 0 || je >= _columnCount) {
		cerr << "Matrix::Insert error, requested submatrix does not exist\n";
		cerr << "  illegal je = " << je << " (matrix has " << _columnCount << " columns)\n";
		return 1;
	  }
	  if (ie < is || je < js) {
		cerr << "Matrix::Insert error, requested submatrix does not exist\n";
		cerr << "  upper index is below lower index: is = " << is << ", ie = " 
		 << ie << ", js = " << js << ", je = " << je << endl;
		return 1;
	  }
	  
	  // perform operation
	  for (Index j=0; j<A._columnCount; j++)  
		for (Index i=0; i<A._rows; i++)  
		  _data[j+js][i+is] = A._data[j][i];
	  
	  // return success
	  return 0;
	}

	//   C = a  (sets all entries of C to the scalar a)
	int Matrix::Constant(MathNumber a) {
	  for (Index j=0; j<_columnCount; j++)  
		for (Index i=0; i<_rows; i++)  
		  _data[j][i] = a;
	  return 0;
	}

	// C = C.^p (component-wise exponentiation of my data to the power p)
	int Matrix::Power(MathNumber p) {
	  for (Index j=0; j<_columnCount; j++)
		for (Index i=0; i<_rows; i++)
		  _data[j][i] = pow(_data[j][i], p);
	  return 0;
	}

	// Cij = |Cij| (component-wise absolute value of my data)
	int Matrix::Abs() {
	  for (Index j=0; j<_columnCount; j++)
		for (Index i=0; i<_rows; i++)
			_data[j][i] = std::abs(_data[j][i]);
	  return 0;
	}

	// Cij = Cji
	int Matrix::Trans() {

	  // perform operation in place if matrix is square
	  if (_rows == _columnCount) {
		for (Index i=0; i<_rows; i++)
		  for (Index j=0; j<i; j++)
		std::swap( _data[j][i], _data[i][j] );

	  // otherwise we need a new data array to ensure a clean transpose
	  } else {

		// create temporary data array
		vector< MathVector > newdata;
		newdata.resize(_rows);
		for (Index j=0; j<_rows; j++)
		  newdata[j].resize(_columnCount);

		// copy transposed data over 
		for (Index j=0; j<_columnCount; j++)
		  for (Index i=0; i<_rows; i++)
		newdata[i][j] = _data[j][i];

		// copy newdata values into existing data array
		_data = newdata;

		// swap matrix dimensions
		std::swap(_rows, _columnCount);
	  }
	  
	  // return success
	  return 0;
	}

	// computes the inverse of a nonsingular matrix 
	int Matrix::Inverse() {

	  // check that matrix sizes match
	  if (_rows != _columnCount) {
		cerr << "Inverse error, non-square matrix\n";
		cerr << "  Matrix is " << _rows << " x " << _columnCount << endl;
		return 1;
	  }

	  // create two temporary matrices for operation
	  Matrix B = Eye(_rows);
	  Matrix A = *this;

	  // call existing LinearSolve routine for computations
	  if (LinearSolve(A, *this, B) != 0)
		return 1;

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


	///// Scalar output operators on matrices /////

	// minimum entry in the matrix
	MathNumber Matrix::Min() const {
	  MathNumber mn=_data[0][0];
	  for (Index j=0; j<_columnCount; j++)
		for (Index i=0; i<_rows; i++)
		  mn = std::min( mn, _data[j][i] );
	  return mn;
	}

	// maximum entry in the matrix
	MathNumber Matrix::Max() const {
	  MathNumber mx=_data[0][0];
	  for (Index j=0; j<_columnCount; j++)
		for (Index i=0; i<_rows; i++)
		  mx = std::max( mx, _data[j][i] );
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
