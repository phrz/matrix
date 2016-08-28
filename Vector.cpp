//
//  Vector.cpp
//  matrix
//
//  Created by Paul Herz on 8/28/16.
//  Copyright © 2016 Paul Herz. All rights reserved.
//

#include "Vector.h"
#include "NotImplementedException.h"

namespace PH {
	
	// Serializable interface
	void Vector::serialize(std::ostream& output) const {
		
		// [CITE] http://stackoverflow.com/a/6693088/3592716
		// quick way to concatenate STL container with commas
		std::ostringstream ss;
		std::copy(_data.begin(), _data.end() - 1, std::ostream_iterator<int>(ss, ", "));
	}
	
	
	static Vector& Vector::deserialize(std::istream& input) {
		throw new NotImplementedException();
	}
	
	
	// standard matrix-vector product
	MathVector matrixVectorProduct(const Matrix& A, const MathVector& v) {
		MathVector res(A.rows(), 0.0);
		if (A.columns() != v.size()) {
			cerr << "matrixVectorProduct: incompatible matrix/vector sizes in A*v\n";
		} else {
			for (Index i=0; i<A.rows(); i++)
		  for (Index j=0; j<A.columns(); j++)
			  res[i] += A(i,j)*v[j];
		}
		return res;
	}
	
	
	MathVector operator*(const Matrix& A, const MathVector& v) {
		return matrixVectorProduct(A,v);
	}
	
	
	// inner product between two vectors
	MathNumber dot(const MathVector& v1, const MathVector& v2) {
		if (v1.size() != v2.size()) {
			cerr << "Dot: incompatible vector sizes\n";
			return 0.0;
		}
		MathNumber res = 0.0;
		for (Index i=0; i<v2.size(); i++)  res += v1[i]*v2[i];
		return res;
	}
	
	
	// vector 2-norm
	MathNumber norm(const MathVector& v) {
		MathNumber sum=0.0;
		for (Index i=0; i<v.size(); i++)  sum += v[i]*v[i];
		return sqrt(sum);
	}
	
	
	// vector infinity norm
	MathNumber infNorm(const MathVector& v) {
		MathNumber mx=0.0;
		for (Index i=0; i<v.size(); i++)
			mx = std::max(mx,std::abs(v[i]));
		return mx;
	}
	
	
	// vector one-norm
	MathNumber oneNorm(const MathVector& v) {
		MathNumber sum=0.0;
		for (Index i=0; i<v.size(); i++)  sum += std::abs(v[i]);
		return sum;
	}
	
	
	// create a new vector of linearly spaced data
	MathVector linSpace(MathNumber a, MathNumber b, Index n) {
		if (n<2) cerr << "Linspace::length must be > 1\n";
		MathVector v(n);
		MathNumber h = (b-a)/(n-1);
		for (Index i=0; i<n; i++)
			v[i] = a + i*h;
		return v;
	}
	
	
	// create a new vector of logarithmically spaced data
	MathVector logSpace(MathNumber a, MathNumber b, Index n) {
		if (n<2) cerr << "Logspace::length must be > 1\n";
		MathVector v(n);
		MathNumber h = (b-a)/(n-1);
		for (Index i=0; i<n; i++)
			v[i] = pow(10.0, a + i*h);
		return v;
	}
	
	
	// create a new vector with uniformly-distributed random numbers in [0,1]
	MathVector randomVectorOfSize(Index n) {
		if (n<1) cerr << "Random::length must be > 0\n";
		MathVector v(n);
		for (Index i=0; i<n; i++)
			v[i] = random() / (pow(2.0,31.0) - 1.0);
		return v;
	}
	
	
	// streaming output routine
	ostream& operator<<(ostream& os, const MathVector& v) {
		for (Index i=0; i<v.size(); i++)
			os << "  " << v[i];
		os << "\n";
		return os;
	}
	
	
	// extract routine for portions of vectors, y = x(is:ie)
	MathVector VecExtract(MathVector& x,
						  long int is, long int ie) {
		
		// update is,ie,js,je if any are negative
		is = (is < 0) ? is+x.size() : is;
		ie = (ie < 0) ? ie+x.size() : ie;
		
		// check that requested subvector exists
		if (is < 0 || is >= x.size()) {
			cerr << "VecExtract error, requested submatrix does not exist\n";
			cerr << "  illegal is = " << is << " (vector has " << x.size() << " entries)\n";
		}
		if (ie < 0 || ie >= x.size()) {
			cerr << "VecExtract error, requested submatrix does not exist\n";
			cerr << "  illegal ie = " << ie << " (matrix has " << x.size() << " entries)\n";
		}
		if (ie < is) {
			cerr << "VecExtract error, requested submatrix does not exist\n";
			cerr << "  upper index is below lower index: is = " << is << ", ie = "
		 << ie << endl;
		}
		
		// create new vector of desired size
		MathVector y(ie-is+1);
		
		// copy requested data
		for (Index i=is; i<=ie; i++)
			y[i-is] = x[i];
		
		// return object
		return y;
	}
	
	
	// insert routine for portions of vectors, x(is:ie) = y
	int VecInsert(MathVector& x, long int is,
				  long int ie, MathVector& y) {
		
		// update is,ie if any are negative
		is = (is < 0) ? is+x.size() : is;
		ie = (ie < 0) ? ie+x.size() : ie;
		
		// check that array sizes match
		if (y.size() != (ie-is+1)) {
			cerr << "VecInsert error, size mismatch\n    supplied vector has " << y.size()
			<< " entries, but requested subvector has " << ie-is+1 << " entries\n";
			return 1;
		}
		// check for valid subvector
		if (is < 0 || is >= x.size()) {
			cerr << "VecInsert error, requested subvector does not exist\n";
			cerr << "  illegal is = " << is << " (vector has " << x.size() << " entries)\n";
			return 1;
		}
		if (ie < 0 || ie >= x.size()) {
			cerr << "VecInsert error, requested subvector does not exist\n";
			cerr << "  illegal ie = " << ie << " (vector has " << x.size() << " entries)\n";
			return 1;
		}
		if (ie < is) {
			cerr << "VecInsert error, requested submatrix does not exist\n";
			cerr << "  upper index is below lower index: is = " << is << ", ie = "
		 << ie << endl;
			return 1;
		}
		
		// perform operation
		for (Index i=0; i<y.size(); i++)
			x[i+is] = y[i];
		
		// return success
		return 0;
	}
	
	
	// arithmetic operators
	MathVector& operator+=(MathVector& v, const MathNumber c) {
		for (Index i=0; i<v.size(); i++)
			v[i] += c;
		return v;
	}
	
	
	MathVector& operator+=(MathVector& v, const MathVector& w) {
		if (v.size() != w.size())
			cerr << "MathVector += error: incompatible vector sizes!";
		else
			for (Index i=0; i<v.size(); i++)
		  v[i] += w[i];
		return v;
	}
	
	
	MathVector& operator-=(MathVector& v, const MathNumber c) {
		for (Index i=0; i<v.size(); i++)
			v[i] -= c;
		return v;
	}
	
	
	MathVector& operator-=(MathVector& v, const MathVector& w) {
		if (v.size() != w.size())
			cerr << "MathVector -= error: incompatible vector sizes!";
		else
			for (Index i=0; i<v.size(); i++)
		  v[i] -= w[i];
		return v;
	}
	
	
	MathVector& operator*=(MathVector& v, const MathNumber c) {
		for (Index i=0; i<v.size(); i++)
			v[i] *= c;
		return v;
	}
	
	
	MathVector& operator*=(MathVector& v, const MathVector& w) {
		if (v.size() != w.size())
			cerr << "MathVector *= error: incompatible vector sizes!";
		else
			for (Index i=0; i<v.size(); i++)
		  v[i] *= w[i];
		return v;
	}
	
	
	MathVector& operator/=(MathVector& v, const MathNumber c) {
		for (Index i=0; i<v.size(); i++)
			v[i] /= c;
		return v;
	}
	
	
	MathVector& operator/=(MathVector& v, const MathVector& w) {
		if (v.size() != w.size())
			cerr << "MathVector /= error: incompatible vector sizes!";
		else
			for (Index i=0; i<v.size(); i++)
		  v[i] /= w[i];
		return v;
	}
	
	
	MathVector& operator^=(MathVector& v, const MathNumber c) {
		for (Index i=0; i<v.size(); i++)
			v[i] = pow(v[i], c);
		return v;
	}
	
	
	MathVector& operator^=(MathVector& v, const MathVector& w) {
		if (v.size() != w.size())
			cerr << "MathVector /= error: incompatible vector sizes!";
		else
			for (Index i=0; i<v.size(); i++)
		  v[i] = pow(v[i], w[i]);
		return v;
	}
	
	
	MathVector operator+(const MathVector& v, const MathNumber c) {
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = v[i] + c;
		return x;
	}
	
	
	MathVector operator+(const MathNumber c, const MathVector& v) {
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = v[i] + c;
		return x;
	}
	
	
	MathVector operator+(const MathVector& v, const MathVector& w) {
		if (v.size() != w.size()) {
			cerr << "MathVector + error: incompatible vector sizes!";
			return MathVector(0);
		}
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = v[i] + w[i];
		return x;
	}
	
	
	MathVector operator-(const MathVector& v, const MathNumber c) {
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = v[i] - c;
		return x;
	}
	
	
	MathVector operator-(const MathNumber c, const MathVector& v) {
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = c - v[i];
		return x;
	}
	
	
	MathVector operator-(const MathVector& v, const MathVector& w) {
		if (v.size() != w.size()) {
			cerr << "MathVector - error: incompatible vector sizes!";
			return MathVector(0);
		}
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = v[i] - w[i];
		return x;
	}
	
	
	MathVector operator*(const MathVector& v, const MathNumber c) {
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = v[i] * c;
		return x;
	}
	
	
	MathVector operator*(const MathNumber c, const MathVector& v) {
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = v[i] * c;
		return x;
	}
	
	
	MathVector operator*(const MathVector& v, const MathVector& w) {
		if (v.size() != w.size()) {
			cerr << "MathVector * error: incompatible vector sizes!";
			return MathVector(0);
		}
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = v[i] * w[i];
		return x;
	}
	
	
	MathVector operator/(const MathVector& v, const MathNumber c) {
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = v[i] / c;
		return x;
	}
	
	
	MathVector operator/(const MathNumber c, const MathVector& v) {
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = c / v[i];
		return x;
	}
	
	
	MathVector operator/(const MathVector& v, const MathVector& w) {
		if (v.size() != w.size()) {
			cerr << "MathVector / error: incompatible vector sizes!";
			return MathVector(0);
		}
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = v[i] / w[i];
		return x;
	}
	
	
	MathVector operator^(const MathVector& v, const MathNumber c) {
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = pow(v[i], c);
		return x;
	}
	
	
	MathVector operator^(const MathNumber c, const MathVector& v) {
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = pow(c, v[i]);
		return x;
	}
	
	
	MathVector operator^(const MathVector& v, const MathVector& w) {
		if (v.size() != w.size()) {
			cerr << "MathVector ^ error: incompatible vector sizes!";
			return MathVector(0);
		}
		MathVector x(v.size());
		for (Index i=0; i<v.size(); i++)
			x[i] = pow(v[i], w[i]);
		return x;
	}
	
} // namespace PH
