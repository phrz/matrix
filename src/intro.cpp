/* Daniel R. Reynolds
   SMU Mathematics
   2 July 2015 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Matrix.h"
using namespace std;


// Example routine to introduce our "Matrix" class
int main(int argc, char* argv[]) {

  cout << endl;

  cout << "\ncreating scalar variables (default to double in this class):\n";
  double a = 5.0;
  cout << "  double a = 5.0;\n";
  cout << "\n a: " << a << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\ncreating row vectors, editing entries, and printing to screen\n";
  Matrix y(1,3);
  y(0,0) = 1.0;
  y(1) = 4.0;
  y(2) = 6.0;
  cout << "  Matrix y(1,3);     // zero-valued matrix with 1 row, 3 columns\n";
  cout << "  y(0,0) = 1.0;      // access entries starting with 0, using () notation\n";
  cout << "  y(1) = 4.0;        // like Matlab, vectors also allow a single index accessor\n";
  cout << "  y(2) = 6.0;\n";
  cout << "  y.Write();         // the 'Write' routine will print to the screen\n";
  cout << "\n y:\n";
  y.Write();
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\ncreating column vectors and printing to file:\n";
  Matrix x(3);
  x(0,0) = 2.0;
  x(1) = 3.0;
  x.Write("x.txt");
  cout << "  Matrix x(3);       // single integer constructor, equivalent to 'Matrix x(3,1)'\n";
  cout << "  x(0,0) = 2.0;      // again, can access (row,column) entry, starting with 0\n";
  cout << "  x(1) = 3.0;        // again, vectors can also access with a single index\n";
  cout << "                     // NOTE: unset entries are zero-valued\n";
  cout << "  x.Write(\"x.txt\");  // if you supply a string to 'Write' it will create/write to that file\n";
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\ncreate Matrix from text file (Matlab/Numpy format):\n";
  Matrix x2 = MatrixRead("x.txt");
  cout << "  Matrix x2 = Read(\"x.txt\");\n";
  cout << "\n x2:\n" << x2;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nother ways to create matrices:\n";
  double Adata[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  Matrix A(2,3,Adata);
  cout << "  double Adata[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};   // standard C/C++ data array\n";
  cout << "  Matrix A(2,3,Adata);  // new 3x2 matrix, with values from a given array \n";
  cout << "  cout << A;            // streaming output for matrices/vectors also works\n";
  cout << "                        // NOTE: column-major ordering of matrix data\n";
  cout << "\n A:\n" << A;
  cout << " Press [enter] to continue\n";  cin.get();

  vector<double> Cve = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  Matrix C(3,2,Cve);
  cout << "  vector<double> Cve = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};   // C++ vector object\n";
  cout << "  Matrix C(3,2,Cve);  // new 3x2 matrix, with values from a given vector \n";
  cout << "\n C:\n" << C;
  cout << " Press [enter] to continue\n";  cin.get();
  
  cout << "\nchanging Matrix entries:\n";
  C(2,1) = 12.0;
  cout << "  C(2,1) = 12.0;      // can change existing matrix entries using ()\n";
  cout << "\n C:\n" << C;
  cout << " Press [enter] to continue\n";  cin.get();
  Cve[5] = -7.0;
  cout << "  Cve[5] = -7.0;      // changing original vector leaves Matrix unaffected\n";
  cout << "\n C:\n" << C;
  cout << " Press [enter] to continue\n";  cin.get();
  C[1][2] = -2.0;
  cout << "  C[1][2] = -2.0;     // can alternately change entries using [] (C/C++ style, column/row order)\n";
  cout << "\n C:\n" << C;
  cout << " Press [enter] to continue\n";  cin.get();
  C(5) = 6.0;
  cout << "  C(5) = 6.0;         // or we can use single-index with () to specify the overall entry\n";
  cout << "\n C:\n" << C;
  cout << " Press [enter] to continue\n";  cin.get();
  
  cout << "\ncreating matrices from strings:\n";
  Matrix D("1, 2, 3; 4, 5, 6");
  cout << "  Matrix D(\"1, 2, 3; 4, 5, 6\");  // new 2x3 matrix, specified by a Matlab-like string\n";
  cout << "\n D:\n" << D;
  cout << " Press [enter] to continue\n";  cin.get();
  
  cout << "\nmatrix arithmetic operations:\n";
  Matrix tmp0 = 10*D;
  cout << "  Matrix tmp0 = 10*D;      // new matrix created from scalar multiplication\n";
  cout << "\n tmp0:\n" << tmp0;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nvector arithmetic operations, first let's make some vectors:\n";
  Matrix u(1,3);
  u(0) = 1.0;  u(1) = 2.0;  u(2) = 3.0;
  Matrix v("4, 5, 6");
  Matrix w("7.0; 8.0; 9.0");
  cout << "  Matrix u(1,3);\n";
  cout << "  u(0) = 1.0;  u(1) = 2.0;  u(2) = 3.0;\n";
  cout << "  Matrix v(\"4, 5, 6\");\n";
  cout << "  Matrix w(\"7.0; 8.0; 9.0\");\n";
  cout << "\n u:\n" << u;
  cout << "\n v:\n" << v;
  cout << "\n w:\n" << w;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nnew matrix created by addition/subtraction:\n";
  cout << "  u-v:\n" << u-v;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nillegal matrix addition:\n";
  cout << "  u+w:\n";
  cout << u+w;
  cout << " Press [enter] to continue\n";  cin.get();
 
  cout << "\n Recall A:\n" << A << endl;
  cout << " Recall C:\n" << C << endl;
  cout << " Recall w:\n" << w << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nnew matrices created by matrix-vector and/or matrix/matrix multiplication:\n";
  cout << "  A*w:\n" << A*w;
  cout << " Press [enter] to continue\n";  cin.get();
  cout << "  A*C:\n" << A*C;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nadditional matrix creation operations:\n";
  Matrix tmp1 = A.T();
  cout << "  Matrix tmp1 = A.T();                   // matrix transpose\n";
  cout << "\n tmp1:\n" << tmp1;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix tmp2 = Linspace(0.0, 1.0, 5, 1);
  cout << "  Matrix tmp2 = Linspace(0.0, 1.0, 5, 1);    // linear span constructor (5x1 result)\n";
  cout << "\n tmp2:\n" << tmp2;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix tmp3 = Logspace(-1.0, 1.0, 1, 5);
  cout << "  Matrix tmp3 = Logspace(-1.0, 1.0, 1, 5);   // log10 span constructor (1x5 result)\n";
  cout << "\n tmp3:\n" << tmp3;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix tmp4 = Random(3,2);
  cout << "  Matrix tmp4 = Random(3,2);              // random matrix constructor\n";
  cout << "\n tmp4:\n" << tmp4;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix tmp5 = Eye(4);
  cout << "  Matrix tmp5 = Eye(4);                   // identity matrix constructor\n";
  cout << "\n tmp5:\n" << tmp5;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nwe can even work with submatrices:\n";
  cout << "\n Recall A:\n" << A << endl;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix tmp6 = A.Extract(0,-1,0,1);
  cout << "Let's create a new matrix out of a submatrix of A:\n\n";
  cout << "  Matrix tmp6 = A.Extract(0,-1,0,1);      // like Matlab's A(0:end,0:1), where\n";
  cout << "                                          // negative indices wrap around from end\n";
  cout << "\n tmp6:\n" << tmp6;
  cout << " Press [enter] to continue\n";  cin.get();
  cout << "Let's now modify tmp6:\n\n";
  tmp6 *= -5.0;
  cout << "  tmp6 *= -5.0;\n";
  cout << "\n tmp6:\n" << tmp6;
  cout << " Press [enter] to continue\n";  cin.get();
  A.Insert(tmp6,0,-1,-2,-1);
  cout << "Let's now insert tmp6 back into a different part of A:\n\n";
  cout << "  A.Insert(tmp6,0,1,-2,-1);               // negative values may be used here too\n";
  cout << "\n A:\n" << A;
  cout << " Press [enter] to continue\n";  cin.get();
  
  cout << "\nmatrix transformation operations (all matrices are assumed to be compatible sizes):\n";
  cout << "\n Recall u:\n" << u;
  cout << "\n Recall v:\n" << v;
  cout << "\n Recall w:\n" << w << endl;
  cout << " Press [enter] to continue\n";  cin.get();
  w.Trans();
  cout << "  w.Trans();                    // transpose-in-place (changes w),\n";
  cout << "                                // unlike T() that created a new Matrix\n";
  cout << "\n w:\n" << w;
  cout << " Press [enter] to continue\n";  cin.get();
  y.LinearSum(1.0, w, -1.0, v);
  cout << "  y.LinearSum(1.0, w, -1.0, v); // in-place linear combination  [y = 1w - 1v]\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  u.LinearSum(2.0, u, 1.0, v);
  cout << "  u.LinearSum(2.0, u, 1.0, v);  // we can also put the output as an input  [u = 2u + 1v]\n";
  cout << "\n u:\n" << u;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Add(u);
  cout << "  y.Add(u);           // adds u to y  [y = y + u]\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Subtract(u);
  cout << "  y.Subtract(u);      // subtracts u from y  [y = y - u]\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Add(5.0);
  cout << "  y.Add(5.0);         // adds 5 to all entries of y  [y = y + 5]\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Multiply(u);
  cout << "  y.Multiply(u);      // in-place product of y and u  [y = y .* u]\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Multiply(2.0);
  cout << "  y.Multiply(2.0);    // in-place scaling of y by 2  [y = 2 * y]\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Divide(u);
  cout << "  y.Divide(u);        // in-place quotient of y and u  [y = y ./ u]\n";
  cout << "                      // (no scalar equivalent, since Multiply(1.0/2.0) works just as well)\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Power(0.5);
  cout << "  y.Power(0.5);       // in-place exponentiation of y  [y = y .^ 0.5]\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Copy(u);
  cout << "  y.Copy(u);          // copies u into y  [y = u]\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Constant(-2.0);
  cout << "  y.Constant(-2.0);   // sets all entries to a constant  [y = -2*ones(size(y))]\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y.Abs();
  cout << "  y.Abs();            // takes the absolute value of y  [y = abs(y)]\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y += u;
  cout << "  y += u;             // shortcut for Add()\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y -= 2.0;
  cout << "  y -= 2.0;           // shortcut for Subtract()\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y *= v;
  cout << "  y *= v;             // shortcut for Multiply()\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y /= v;
  cout << "  y /= v;             // shortcut for Divide()\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y /= 3.0;
  cout << "  y /= 3.0;           // scalar /= operator\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y ^= 2.0;
  cout << "  y ^= 2.0;           // shortcut for Power()\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();
  y = 3.0;
  cout << "  y = 3.0;            // shortcut for Constant()\n";
  cout << "                      // (NOTE: there is no  y = u  'copy' operator!  While y = u\n";
  cout << "                      //  compiles/runs, it destroys/allocates memory (very slow))\n";
  cout << "\n y:\n" << y;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nMatrix/vector inquiry functions:\n";
  cout << "\nRecall A:\n" << A;
  cout << "\n and u:\n" << u;
  cout << "\n and v:\n" << v;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nmatrix dimensions:\n";
  cout << "  A.Rows():    " << A.Rows() << endl;
  cout << "  A.Columns(): " << A.Columns() << endl;
  cout << "  A.Size():    " << A.Size() << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nvector 2-norm, matrix Frobenius norm:\n";
  cout << "  Norm(v): " << Norm(v) << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nmatrix one-norm, column vector one-norm, row vector inf-norm:\n";
  cout << "  OneNorm(v): " << OneNorm(v) << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nmatrix inf-norm, column vector inf-norm, row vector one-norm:\n";
  cout << "  InfNorm(v): " << InfNorm(v) << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nequality comparison operator:\n";
  cout << "  u == v: " << (u == v) << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nsmallest entry of v:\n";
  cout << "  v.Min(): " << v.Min() << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nlargest entry of v:\n";
  cout << "  v.Max(): " << v.Max() << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\ndot-product of two vectors/matrices (must be the same shape):\n";
  cout << "  Dot(u,v): " << Dot(u,v) << endl;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nWe even have direct linear solvers (FwdSub, BackSub, GEPP):\n";
  Matrix E = Random(4,4);
  Matrix E2 = E;
  cout << "  Matrix E = Random(4,4);\n";
  cout << "  Matrix E2 = E;\n";
  cout << "\n E:\n" << E;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix b = Random(4,1);
  Matrix b2 = b;
  cout << "  Matrix b = Random(4,1);\n";
  cout << "  Matrix b2 = b;\n";
  cout << "\n b:\n" << b;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix s = LinearSolve(E,b);
  cout << "  Matrix s = LinearSolve(E,b);          // GEPP solver that allocates solution memory\n";
  cout << "                                        // (modifies E and b during solution process)\n";
  cout << "\n s:\n" << s;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix s2(4,1);
  LinearSolve(E2,s2,b2);
  cout << "  Matrix s2(4,1);\n";
  cout << "  LinearSolve(E2,s2,b2);                // GEPP solver that uses existing solution memory\n";
  cout << "                                        // (modifies E2 and b2 during solution process)\n";
  cout << "\n s2:\n" << s2 << endl;
  cout << " Press [enter] to continue\n";  cin.get();
  cout << "\nBackSubstitution() and ForwardSubstitution() work similarly\n";

  cout << "\nMatrices are stored as  vector< vector<double> >  objects, i.e. each column is a vector<double>,\n";
  cout << "so we can access these columns directly:\n";
  cout << "\nRecall A:\n" << A;
  cout << " Press [enter] to continue\n";  cin.get();
  cout << "cout << A[2];                            // column accessor (the reason A[2][1] works)\n";
  cout << "                                         // (NOTE: vectors print horizontally)\n";
  cout << A[2];
  cout << " Press [enter] to continue\n";  cin.get();
  vector<double> tmp7 = A.Column(1);
  cout << "  vector<double> tmp7 = A.Column(1);     // extract/copy column \n";
  cout << "\n tmp7:\n" << tmp7;
  cout << " Press [enter] to continue\n";  cin.get();
  vector<double> tmp8 = A.Row(0);
  cout << "  vector<double> tmp8 = A.Row(0);        // extract/copy row (copies entries from column vectors)\n";
  cout << "\n tmp8:\n" << tmp8;
  cout << " Press [enter] to continue\n";  cin.get();

  cout << "\nWe also provide a subset of this functionality for vector<double> objects:\n";
  vector<double> p = Linspace(0.0, 2.0, 5);
  vector<double> q = Logspace(-2.0, 2.0, 5);
  vector<double> r = Random(5);
  cout << "  vector<double> p = Linspace(0.0, 2.0, 5)\n";
  cout << "  vector<double> q = Logspace(-2.0, 2.0, 5)\n";
  cout << "  vector<double> r = Random(5)\n";
  cout << "\n p:\n" << p;
  cout << "\n q:\n" << q;
  cout << "\n r:\n" << r;
  cout << " Press [enter] to continue\n";  cin.get();
  cout << "  Norm(p): " << Norm(p) << endl;
  cout << "  OneNorm(q): " << OneNorm(q) << endl;
  cout << "  InfNorm(r): " << InfNorm(r) << endl;
  cout << "  Dot(p,q): " << Dot(p,q) << endl;
  cout << " Press [enter] to continue\n";  cin.get();
  vector<double> psub = VecExtract(p, 2, -1);
  cout << "  vector<double> psub = VecExtract(p, 2, -1)\n";
  cout << "\n psub:\n" << psub;
  cout << " Press [enter] to continue\n";  cin.get();
  VecInsert(q, 1, -2, psub);
  cout << "  VecInsert(q, 1, -2, psub);\n";
  cout << "\n q:\n" << q;
  cout << " Press [enter] to continue\n";  cin.get();
  Matrix F = Random(4,4);
  Matrix F2 = F;
  cout << "  Matrix F = Random(4,4);\n";
  cout << "  Matrix F2 = F;\n";
  cout << "\n F:\n" << F;
  cout << " Press [enter] to continue\n";  cin.get();
  vector<double> rhs = Random(4);
  vector<double> rhs2 = rhs;
  cout << "  vector<double> rhs = Random(4);\n";
  cout << "  vector<double> rhs2 = rhs;\n";
  cout << "\n rhs:\n" << rhs;
  cout << " Press [enter] to continue\n";  cin.get();
  vector<double> sol = LinearSolve(F,rhs);
  cout << "  vector<double> sol = LinearSolve(F,rhs);  // (modifies F and rhs during solution process)\n";
  cout << "\n sol:\n" << sol;
  cout << " Press [enter] to continue\n";  cin.get();
  vector<double> sol2(4);
  LinearSolve(F2,sol2,rhs2);
  cout << "  vector<double> sol2(4);\n";
  cout << "  LinearSolve(F2,sol2,rhs2);  // (modifies F2 and rhs2 during solution process)\n";
  cout << "\n sol2:\n" << sol2 << endl;
  
  return 0;
} // end main

