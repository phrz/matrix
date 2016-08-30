//
//  PHMath.cpp
//  matrix
//
//  Created by Paul Herz on 8/28/16.
//  Copyright © 2016 Paul Herz. All rights reserved.
//

//#include <random>
//
//namespace PH {
//
//	// create a matrix with uniformly-distributed random numbers in [0,1]
//	Matrix randomMatrixOfSize(Index m, Index n) {
//		Matrix C(m,n);
//		for (Index j=0; j<n; j++)
//			for (Index i=0; i<m; i++)
//				C._data[j][i] = random() / (pow(2.0,31.0) - 1.0);
//		return C;
//	}	
//	
//	// create a new vector with uniformly-distributed random numbers in [0,1]
//	Vector randomVectorOfSize(Index n) {
//		if (n<1) {
//			throw new std::invalid_argument("randomVectorOfSize expects a vector size argument > 0.");
//		}
//		
//		// Create random device for a uniform integer
//		// distribution.
//		
//		// This is much like MATLAB's rand().
//		std::random_device randomDevice;
//		std::mt19937 generator(randomDevice());
//		std::uniform_real_distribution<> dist(0, 1);
//		
//		Vector result(n);
//		result.mapElements([dist, generator](MathNumber& element, Index i){
////			element = dist(generator);
//			element = std::rand();
//		});
//		
//		return result;
//	}
//	
//	
//	// standard matrix-vector product
//	Vector matrixVectorProduct(const Matrix& A, const Vector& v) {
//		
//		Vector result(A.rows(), 0.0);
//		
//		if (A.columns() != v.size()) {
//			throw new std::invalid_argument("matrix-vector product: incompatible matrix/vector sizes in A*v");
//		}
//		
//		for(Index column = 0; column < A.columns(); ++column) {
//			for(Index row = 0; row < A.rows(); ++row) {
//				 result += A(row,column) * v[column];
//			}
//		}
//		
//		return result;
//	}
//	
//	
//	Vector operator*(const Matrix& A, const Vector& v) {
//		return matrixVectorProduct(A,v);
//	}
//
//
//} // namespace PH
