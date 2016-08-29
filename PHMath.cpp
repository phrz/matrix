//
//  PHMath.cpp
//  matrix
//
//  Created by Paul Herz on 8/28/16.
//  Copyright © 2016 Paul Herz. All rights reserved.
//

#include "PHMath.h"


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
Matrix linSpace(MathNumber a, MathNumber b, Index m, Index n) {
	Matrix C(m,n);
	MathNumber h = (b-a)/(m*n-1);
	Index idx=0;
	for (Index j=0; j<n; j++)
		for (Index i=0; i<m; i++)
			C._data[j][i] = a + (idx++)*h;
	return C;
}

// create a new column-vector matrix of logarithmically spaced data
Matrix logSpace(MathNumber a, MathNumber b, Index m, Index n) {
	Matrix C(m,n);
	MathNumber h = (b-a)/(m*n-1);
	Index idx=0;
	for (Index j=0; j<n; j++)
		for (Index i=0; i<m; i++)
			C._data[j][i] = pow(10.0, a + (idx++)*h);
	return C;
}

// create a matrix with uniformly-distributed random numbers in [0,1]
Matrix randomMatrixOfSize(Index m, Index n) {
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

// backward substitution on the linear system U*X = B, filling in a Matrix X
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
std::optional<Vector> backSubstitution(const Matrix& U, const Vector& b) {
	
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
	if (std::abs(A._data[A.rows()-1][A.rows()-1]) < STOL*Amax) {
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
	if (std::abs(A._data[A.rows() - 1][A.rows() - 1]) < STOL * Amax) {
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
