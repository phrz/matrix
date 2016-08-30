/* Daniel R. Reynolds
   SMU Mathematics
   19 June 2015 */

#include <stdlib.h>
#include <stdio.h>

#include "Matrix.h"
#include "Vector.h"

using namespace PH;


// Gram-Schmidt process for orthonormalizing a set of vectors
int GramSchmidt(Matrix& X) {

  // check that there is work to do
  if (X.columns() < 1)  return 0;

  // get entry magnitude (for linear dependence check)
	double Xmax = Matrix::infNorm(X);

  // normalize first column
	double colnorm = Vector::norm(Vector(X.row(0)));
  if (colnorm < 1.e-13*Xmax) {
    fprintf(stderr,"GramSchmidt error: vectors are linearly-dependent!\n");
    return 1;
  }
  X(0) *= 1.0/colnorm;

  // iterate over remaining vectors, performing Gram-Schmidt process
  for (int i=1; i<X.columns(); i++) {
    
    // subtract off portions in directions of existing basis vectors
    for (int j=0; j<i; j++)
		X(j) -= Matrix::dot(X(i), X(j)) * X(j);

    // normalize vector, checking for linear dependence
	  colnorm = Matrix::norm(X(i));
    if (colnorm < 1.e-13*Xmax) {
      fprintf(stderr,"GramSchmidt error: vectors are linearly-dependent!\n");
      return 1;
    }
    X(i) *= 1.0/colnorm;
  }

  // return success
  return 0;
}
