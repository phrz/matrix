
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
		if (r*c != source.size()) {
			cerr << "Matrix constructor error: incompatible shape with vector length\n";
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
				cerr << "Matrix constructor error: rows in 2D vector must have the same length\n";
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


	///// Arithmetic operations defined on a given Mat

	// C = A*X
	int Matrix::Product(const Matrix& A, const Matrix& X) {

	  // check that array sizes match
	  if (A._rows != _rows || A._columnCount != X._rows || X._columnCount != _columnCount) {
		cerr << "Matrix::Product error, matrix size mismatch\n";
		cerr << "  Matrix 1 is " << A._rows << " x " << A._columnCount 
		 << "  Matrix 2 is " << X._rows << " x " << X._columnCount 
		 << ", output is " << _rows << " x " << _columnCount << endl;
		return 1;
	  }
	  
	  // perform operation
	  this->Constant(0.0);
	  for (Index k=0; k<_columnCount; k++) 
		for (Index j=0; j<A.columns(); j++)
		  for (Index i=0; i<A.rows(); i++) 
		(*this)(i,k) += A(i,j)*X(j,k);
		
	  // return success
	  return 0;

	}


	// C = a*A + b*B
	int Matrix::LinearSum(MathNumber a, const Matrix& A, MathNumber b, const Matrix& B) {

	  // check that array sizes match
	  if (A._rows != _rows || A._columnCount != _columnCount || 
		  B._rows != _rows || B._columnCount != _columnCount) {
		cerr << "Matrix::LinearSum error, matrix size mismatch\n";
		cerr << "  Matrix 1 is " << _rows << " x " << _columnCount 
		 << ",  Matrix 2 is " << A._rows << " x " << A._columnCount 
		 << ",  Matrix 3 is " << B._rows << " x " << B._columnCount << endl;
		return 1;
	  }
	  
	  // perform operation
	  for (Index j=0; j<B._columnCount; j++)
		for (Index i=0; i<B._rows; i++)
		  _data[j][i] = a*A._data[j][i] + b*B._data[j][i];
	  
	  // return success
	  return 0;
	}

	// C = C+a  (adds scalar a to my data)
	int Matrix::Add(MathNumber a) {
	  for (Index j=0; j<_columnCount; j++)
		for (Index i=0; i<_rows; i++)
		  _data[j][i] += a;
	  return 0;
	}

	// C = C.*A (component-wise multiply of my data by A)
	int Matrix::Multiply(const Matrix& A) {

	  // check that array sizes match
	  if (A._rows != _rows || A._columnCount != _columnCount) {
		cerr << "Matrix::Multiply error, matrix size mismatch\n";
		cerr << "  Matrix 1 is " << _rows << " x " << _columnCount 
		 << ",  Matrix 2 is " << A._rows << " x " << A._columnCount << endl;
		return 1;
	  }
	  
	  // perform operation
	  for (Index j=0; j<A._columnCount; j++)
		for (Index i=0; i<A._rows; i++)
		  _data[j][i] *= A._data[j][i];
	  
	  // return success
	  return 0;
	}

	// C = a*C  (scales my data by scalar a)
	int Matrix::Multiply(MathNumber a) {

	  // perform operation
	  for (Index j=0; j<_columnCount; j++)
		for (Index i=0; i<_rows; i++)
		  _data[j][i] *= a;
	  
	  // return success
	  return 0;
	}

	// C = C./A (component-wise division of my data by A)
	int Matrix::Divide(const Matrix& A) {

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
		  _data[j][i] = abs(_data[j][i]);
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

	// equivalence-checking operator
	bool Matrix::operator==(const Matrix& A) const {

	  // quick check for compatible sizes
	  if (A._rows != _rows || A._columnCount != _columnCount)
		return false;

	  // detailed check on values
	  bool equal = true;
	  for (Index j=0; j<_columnCount; j++)  
		for (Index i=0; i<_rows; i++)
		  equal &= (A._data[j][i] == _data[j][i]);
	  return equal;
	}


	//--- supplementary matrix arithmetic routines ---

	// column/row vector dot product, something else for matrices?
	MathNumber Dot(const Matrix& A, const Matrix& B) {
	  MathNumber sum=0.0;
	  if ((A.columns() != B.columns()) || (A.rows() != B.rows())) {
		cerr << "Dot error, matrix objects must be the same size\n";
	  } else {
		for (Index j=0; j<A.columns(); j++) 
		  for (Index i=0; i<A.rows(); i++)
		sum += A(i,j) * B(i,j);
	  }
	  return sum;
	}

	// matrix Frobenius norm (column/row vector 2-norm)
	MathNumber Norm(const Matrix& A) {
	  MathNumber sum=0.0;
	  for (Index j=0; j<A.columns(); j++) 
		for (Index i=0; i<A.rows(); i++)
		  sum += A(i,j)*A(i,j);
	  return sqrt(sum);
	}

	// matrix infinity norm (column vector infinity norm, row vector one norm)
	MathNumber InfNorm(const Matrix& A) {
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
	MathNumber OneNorm(const Matrix& A) {
	  MathNumber mx=0.0;
	  for (Index j=0; j<A.columns(); j++) {
		MathNumber sum=0.0;
		for (Index i=0; i<A.rows(); i++)
		  sum += std::abs(A(i,j));
		mx = std::max(mx,sum);
	  }
	  return mx;
	}


	//--- new matrix creation routines (C is the output, A and B are the operands) ---

	// C = A+B
	Matrix operator+(const Matrix& A, const Matrix& B) {

	  // check that array sizes match
	  if (B.rows() != A.rows() || B.columns() != A.columns()) {
		cerr << "Matrix operator+ error, matrix size mismatch\n";
		cerr << "  Matrix 1 is " << A.rows() << " x " << A.columns()
		 << ",  Matrix 2 is " << B.rows() << " x " << B.columns() << endl;
		return Matrix(0,0);
	  }

	  // create new Mat for output, and do operation
	  Matrix C(A.rows(),A.columns());
	  for (Index j=0; j<A.columns(); j++)
		for (Index i=0; i<A.rows(); i++)
		  C._data[j][i] = A._data[j][i] + B._data[j][i];

	  // return result
	  return C;
	}

	// C = A-B
	Matrix operator-(const Matrix& A, const Matrix& B) {

	  // check that array sizes match
	  if (B.rows() != A.rows() || B.columns() != A.columns()) {
		cerr << "Matrix operator- error, matrix size mismatch\n";
		cerr << "  Matrix 1 is " << A.rows() << " x " << A.columns()
		 << ",  Matrix 2 is " << B.rows() << " x " << B.columns() << endl;
		return Matrix(0,0);
	  }

	  // create new Mat for output, and do operation
	  Matrix C(A.rows(),A.columns());
	  for (Index j=0; j<A.columns(); j++)
		for (Index i=0; i<A.rows(); i++)
		  C._data[j][i] = A._data[j][i] - B._data[j][i];

	  // return result
	  return C;
	}

	// C = A*B
	Matrix operator*(const Matrix& A, const Matrix& B) {

	  // determine if either matrix is a scalar
	  bool A_scalar = ((A.rows()==1) && (A.columns()==1));
	  bool B_scalar = ((B.rows()==1) && (B.columns()==1));

	  // scalar-times-matrix
	  if (A_scalar) {

		// create new Matrix for output, and do operation
		Matrix C(B.rows(),B.columns());
		for (Index j=0; j<B.columns(); j++)
		  for (Index i=0; i<B.rows(); i++)
		C._data[j][i] = A._data[0][0] * B._data[j][i];
		return C;
		
	  // scalar-times-matrix
	  } else if (B_scalar) {
		
		// create new Mat for output, and do operation
		Matrix C(A.rows(),B.columns());
		for (Index j=0; j<A.columns(); j++)
		  for (Index i=0; i<A.rows(); i++)
		C._data[j][i] = A._data[j][i] * B._data[0][0];
		return C;
		
	  // normal matrix product
	  } else {

		// check that array sizes are acceptable
		if (B.rows() != A.columns()) {
		  cerr << "Matrix operator* error, inner dimension mismatch\n";
		  cerr << "  Matrix 1 is " << A.rows() << " x " << A.columns()
		   << ",  Matrix 2 is " << B.rows() << " x " << B.columns() << endl;
		  return Matrix(0,0);
		} else {
		
		  // create new Mat for output, and do operation
		  Matrix C(A.rows(),B.columns());
		  for (Index i=0; i<A.rows(); i++) 
		for (Index j=0; j<B.columns(); j++)
		  for (Index k=0; k<A.columns(); k++)
			C._data[j][i] += A._data[k][i] * B._data[j][k];
		  return C;
		}
	  }
	}


	// C = A*b
	Matrix operator*(const Matrix& A, const MathNumber b) {
	  Matrix C(A.rows(),A.columns());
	  for (Index j=0; j<A.columns(); j++)
		for (Index i=0; i<A.rows(); i++)
		  C._data[j][i] = A._data[j][i] * b;
	  return C;
	}

	// C = a*B
	Matrix operator*(const MathNumber a, const Matrix& B) {
	  Matrix C(B.rows(),B.columns());
	  for (Index j=0; j<B.columns(); j++)
		for (Index i=0; i<B.rows(); i++)
		  C._data[j][i] = a * B._data[j][i];
	  return C;
	}

	// create a new matrix of linearly spaced data
	Matrix Linspace(MathNumber a, MathNumber b, Index m, Index n) {
	  Matrix C(m,n);
	  MathNumber h = (b-a)/(m*n-1);
	  Index idx=0;
	  for (Index j=0; j<n; j++)
		for (Index i=0; i<m; i++)
		  C._data[j][i] = a + (idx++)*h;
	  return C;
	}

	// create a new column-vector matrix of logarithmically spaced data
	Matrix Logspace(MathNumber a, MathNumber b, Index m, Index n) {
	  Matrix C(m,n);
	  MathNumber h = (b-a)/(m*n-1);
	  Index idx=0;
	  for (Index j=0; j<n; j++)
		for (Index i=0; i<m; i++)
		  C._data[j][i] = pow(10.0, a + (idx++)*h);
	  return C;
	}

	// create a matrix with uniformly-distributed random numbers in [0,1]
	Matrix Random(Index m, Index n) {
	  Matrix C(m,n);
	  for (Index j=0; j<n; j++)
		for (Index i=0; i<m; i++)
		  C._data[j][i] = random() / (pow(2.0,31.0) - 1.0);
	  return C;
	}

	// create a new n by n identity matrix
	Matrix Eye(Index n) {
	  Matrix I(n,n);
	  for (Index i=0; i<n; i++)  I._data[i][i] = 1.0;
	  return I;
	}

	// creates a matrix from a specified input file
	Matrix MatrixRead(const char *infile) {
	  
	  // determine matrix size
	  Index _rows=0, _columnCount=0;
	  ifstream ifs;
	  string line;
	  ifs.open(infile);
	  while (getline(ifs, line)) {
		istringstream iss(line);   // convert line to stringstream
		float value;               // determine the number of columns on this row
		Index n=0;
		while (iss >> value)  n++;
		if ((n > 0) && (_rows == 0))  // first row, set _columnCount
		  _columnCount = n;
		if ((n > 0) && (n != _columnCount)) {  // later row, with bad number of columns
		  cerr << "MatrixRead() error, not all rows in file " << infile 
		   << " have the same number of cols, "
		   << n << " != " << _columnCount << endl;
		  return Matrix(0,0);
		}
		if (n > 0) _rows++;          // legal row, increment counter
	  }
	  ifs.close();

	  // create matrix of desired size
	  Matrix A(_rows,_columnCount);

	  // load matrix based on data from file
	  ifs.open(infile);   // reopen input file
	  for (Index i=0; i<_rows; i++) {
		getline( ifs, line ); 
		istringstream iss(line);   // convert line to stringstream
		for (Index j=0; j<_columnCount; j++) 
		  iss >> A._data[j][i];
	  }
	  ifs.close();

	  // return result
	  return A;
	}

	// computes the inverse of a nonsingular matrix A
	Matrix inverse(const Matrix& A) {
	  
	  // copy A into a new output matrix
	  Matrix X(A);

	  // call existing Inverse routine for computations
	  if (X.Inverse() != 0)
		cerr << "Inverse: error inverting matrix\n";

	  // return result
	  return X;
	}

	//--- linear algebra routines ---

	// backward substitution on the linear system U*X = B, filling in an existing Matrix X
	std::optional<Matrix> backSubstitution(const Matrix& U, const Matrix& B) {

		auto X = Matrix(U.rows(), B.columns());
		std::optional<Matrix> emptyResult;

		// check that matrix sizes match
		if (U.rows() != B.rows() || not U.isSquare()) {
			cerr << "BackSubstitution error, incompatible matrix/vector dimensions\n";
			cerr << "  Matrix is " << U.rows() << " x " << U.columns() 
			<< ",  rhs is " << B.rows() << " x " << B.columns()
			<< ",  solution is " << X.rows() << " x " << X.columns() << endl;
			return emptyResult;
		}

		// copy B into X
		X = B;

		// analyze matrix for magnitude
		MathNumber Umax = InfNorm(U);

		// perform column-oriented Backward Substitution algorithm
		for (Index j=U.rows()-1; j>=0; j--) {
			// check for nonzero matrix diagonal
			if (abs(U._data[j][j]) < STOL*Umax) {
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
	std::optional<Vector> backSubstitution(const Matrix& U, const MathVector& b) {
		
		auto x = Vector(U.rows());
		std::optional<Vector> emptyResult;
		
		// check that matrix sizes match
		if (U.rows() != b.size() || U.rows() != U.columns()) {
			cerr << "BackSubstitution error, incompatible matrix/vector dimensions\n";
			cerr << "  Matrix is " << U.rows() << " x " << U.columns() 
			<< ",  rhs is " << b.size() << " x 1"
			<< ",  solution is " << x.size() << " x 1\n";
			return emptyResult;
		}

		// copy b into x
		x = b;

		// analyze matrix for magnitude
		MathNumber Umax = InfNorm(U);

		// perform column-oriented Backward Substitution algorithm
		for (long int j=U.rows()-1; j>=0; j--) {
			// check for nonzero matrix diagonal
			if (abs(U._data[j][j]) < STOL * Umax) {
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
	std::optional<Matrix> forwardSubstitution(const Matrix& L, const Matrix& B) {

		auto X = Matrix(L.rows(), B.columns());
		std::optional<Matrix> emptyResult;

		// check that matrix sizes match
		if (L.rows() != B.rows() or not L.isSquare()) {
			cerr << "ForwardSubstitution error, illegal matrix/vector dimensions\n";
			cerr << "  Matrix is " << L.rows() << " x " << L.columns() 
			<< ",  rhs is " << B.rows() << " x " << B.columns()
			<< ",  solution is " << X.rows() << " x " << X.columns() << endl;
			return emptyResult;
		}

		// copy B into X
		X = B;

		// analyze matrix magnitude
		MathNumber Lmax = InfNorm(L);

		// perform column-oriented Forwards Substitution algorithm
		for (long int j=0; j<L.rows(); j++) {
			// check for nonzero matrix diagonal
			if (fabs(L._data[j][j]) < STOL*Lmax) {
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
	std::optional<Vector> forwardSubstitution(const Matrix& L, const Vector& b) {
		
		auto x = Vector(L.rows());
		std::optional<Vector> emptyResult;
		
		// check that matrix sizes match
		if (L.rows() != b.size() || L.rows() != L.columns()) {
			cerr << "ForwardSubstitution error, illegal matrix/vector dimensions\n";
			cerr << "  Matrix is " << L.rows() << " x " << L.columns() 
			<< ",  rhs is " << b.size() << " x 1"
			<< ",  solution is " << x.size() << " x 1\n";
			return emptyResult;
		}
		
		x = b;

		// analyze matrix for magnitude
		MathNumber Lmax = infNorm(L);
		
		// perform column-oriented Forwards Substitution algorithm
		for (long int j=0; j<L.rows(); j++) {

			// check for nonzero matrix diagonal
			if (abs(L._data[j][j]) < STOL*Lmax) {
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
		return std::optional<Vector>(x);
	}

	// solves a linear system A*X = B, returning Matrix X
	std::optional<Matrix> linearSolve(Matrix& A, Matrix& B) {
		
		auto X = Matrix(A.rows(), B.columns());
		std::optional<Matrix> emptyResult;
		
		// check that matrix sizes match
		if (A.rows() != B.rows() or not A.isSquare()) {
			cerr << "linearSolve: illegal matrix/vector dimensions\n";
			cerr << "  Matrix is " << A.rows() << " x " << A.columns() 
			<< ",  rhs is " << B.rows() << " x " << B.columns()
			<< ",  solution is " << X.rows() << " x " << X.columns() << endl;
			return emptyResult;
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
				if (fabs(A._data[k][i]) > std::abs(A._data[k][p])) {
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
			if (abs(A._data[k][k]) < STOL*Amax) {
				cerr << "linearSolve: numerically singular matrix!\n";
				return emptyResult;
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
		if (abs(A._data[A.rows()-1][A.rows()-1]) < STOL*Amax) {
			cerr << "linearSolve error: numerically singular matrix!\n";
			return emptyResult;
		}

		// perform Backward Substitution on result
		if (backSubstitution(A, X, B) != 0) {
			cerr << "linearSolve: error in backSubstitution call\n";
			return emptyResult;
		}

		// return success
		return std::optional<Matrix>(X);
	}

	// solves a linear system A*x = b, returning Vector x
	//    A and b are modified in this operation; x holds the result
	std::optional<Vector> linearSolve(Matrix& A, Vector& b) {
		
		auto x = Vector(A.columns());
		std::optional<Vector> emptyResult;
		
		// check that matrix sizes match
		if (A.rows() != b.size() or not A.isSquare()) {
			cerr << "linearSolve: illegal matrix/vector dimensions\n";
			cerr << "  Matrix is " << A.rows() << " x " << A.columns()
			<< ",  rhs is " << b.size() << " x 1"
			<< ",  solution is " << x.size() << " x 1\n";
			return emptyResult;
		}
		
		// create temporary variables
		Index i, j, k, p;
		MathNumber Amax;
		
		// determine magnitude of entries in A (for singularity
		// check later)
		Amax = InfNorm(A);
		
		// perform Gaussian elimination to convert A,B to an
		// upper-triangular system
		for (k = 0; k < A.rows() - 1; k++) {	// loop over diagonals
			
			// find the pivot row p
			p = k;
			
			for (i = k; i < A.rows(); i++) {
				if (abs(A._data[k][i]) > abs(A._data[k][p])) {
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
			if (abs(A._data[k][k]) < STOL * Amax) {
				cerr << "LinearSolve error: numerically singular matrix!\n";
				return 1;
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
		if (abs(A._data[A.rows() - 1][A.rows() - 1]) < STOL * Amax) {
			cerr << "linearSolve: numerically singular matrix!\n";
			return emptyResult;
		}
		// perform Backward Substitution on result
		if (BackSubstitution(A, x, b) != 0) {
			cerr << "LinearSolve: error in BackSubstitution call\n";
			return emptyResult;
		}
		
		return std::optional<Vector>(x);
	}
	
} // namespace PH
