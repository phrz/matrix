# Matrix
by Paul Herz

I developed this library off of Dr. Reynolds' implementation. It has been approved as a drop-in replacement for that library in the course of the High Performance Scientific Computing class.

## Setup
To use the included `Matrix` and `Vector` classes in your project, simply drop the contents of the `src` folder here into your project. You may want to separate this code into a `lib` folder to denote that it is not your project code if using it in HPSC.

### Including in a Make project
1. Copy the contents of `src` in this repo, including all headers, into a location in your project.
2. Add `Matrix.cpp` and `Vector.cpp` to your build list (e.g. `g++ yourfile.cpp otherfile.cpp Matrix.cpp Vector.cpp`, including the relative path if they are in a different folder than your code.
3. If you've placed this code in a different folder than where you are running the compiler, add the include header folder flag (`-I /lib` if the files are in `lib`).
4. You will need to use at least the C++11 standard, specifically the GNU extension, using the flag `-std=gnu++11`. Your project may not work if it is using an older standard or a non-GNU extension standard. If you have trouble with a *newer* standard, file an issue and I will respond promptly.

### Including in an Xcode project
This example starts from absolute scratch in Xcode 9.0. This version changed the process of including files to be a little more complicated and restrictive. As such, it is difficult to actually store this code in a separate folder from your own. However, we can still use Xcode's Groups (pseudo-folders in the Xcode sidebar) to *visually* group them apart from your own code.

1. Open Xcode to the welcome screen, **Create a New Xcode Project**. If you don't get a welcome screen, go to `File > New > Project...`, or press `Shift + Cmd + N`.
2. Select **macOS**, **Command Line Tool** under the **Application** category. Click **Next**.
3. Enter a product name. Select **C++** as your language. Click **Next**.
4. Pick a location to save your project to. Be wary of the **Create Git repository on my Mac**. You may or may not want this checked. Xcode will create a folder named after your project here.
5. Open up a Finder window next to Xcode. In that window, browse to where you've downloaded this code, specifically to the `src` folder. Select all of the contents of this folder (not the folder itself) and drag into the Xcode file sidebar, dropping them directly above or below the `main.cpp` file. 
6. In the popup, select **Copy items if needed** and **Create groups**. Click **Finish**.
7. These files may clutter your future or existing code. In the Xcode file toolbar (*Not in Finder*), select all the files belonging to the Matrix library by clicking the first file, holding `Shift`, and clicking the last. Right click, and select **New Group from Selection**. Give it a name like `matrix` or `lib`, then press `Return`. If you don't want this group *inside* of your `ProjectName` group, you can carefully drag it above the `ProjectName` group on the sidebar in Xcode, making it a sibling group rather than a child group.
8. If you look in Finder and there is one copy of all the Matrix files inside /ProjectName/ProjectName and another in /ProjectName/ProjectName/GroupName, you will need to manually delete the old copies (in /ProjectName/ProjectName). This is a behavioral issue in Xcode 9.0.
9. Select your project from the file sidebar in Xcode. Select the **Build Settings** tab. Search `C++ Language Dialect` in the search bar below the tab bar.
10. Make sure the selection in the dropdown is *at least* **GNU++11 [-std=gnu++11]**. This library depends on at least C++11, and has been tested with only the C++11 standard with GNU extensions. Further testing is to be done.
11. Select the **Build Phases** tab now under the same project, and open the **Compile Sources** dropdown by clicking on the arrow next to it. Make sure `Matrix.cpp` and `Vector.cpp` are added. If they are not, click the plus and add them.
12. You should be able to add the following to your `main.cpp` file or any other file in the project:
```
#include "Matrix.h"
#include "Vector.h"
```
Include one, the other, or both as needed.

---

## Use
This project includes two classes, `PH::Matrix` and `PH::Vector`. The Matrix class represents a two-dimensional, `m,n` matrix of `doubles`. The Vector class is only for one-dimensional vectors of `doubles`. The reason this library has a separate Vector class while the original doesn't is to avoid accidentally providing a 2D matrix in places where only a vector suffices.

The `PH` namespace avoids colliding with anyone else's matrix code (including Reynolds'). You can reference these classes one of two ways:
```cpp
// the long way
#include "Matrix.h"
int main() {
	auto myMatrix = PH::Matrix();
}

// the short way
#include "Matrix.h"
using namespace PH;
int main() {
	auto myMatrix = Matrix();
}
```

The rest of this document will assume you're doing it the "quick" way, with the `using` statement.

### Type aliases
You will probably not see `doubles` mentioned directly in the code for values, or `size_t` values mentioned directly for dimensions/indexing. These types are aliased to `PH::MathNumber` and `PH::Index`, respectively. This is so they can be tweaked in one place rather than going through to edit them in 200 places.

### Creating basic matrices
You can create a Matrix one of many ways:

- A (0,0) matrix (presumably to be resized later) with `Matrix();`.
- A zero-filled (r,c) matrix: `Matrix(Index r, Index c);` (see Type aliases, above).
- With a C++ initializer list: `Matrix({{1,2,3},{4,5,6},{7,8,9}});` (for you MATLAB folks. Fails if misshapen).
- Again with an initializer list: `Matrix myMatrix = {{1,2,3},{4,5,6},{7,8,9}};`.

And in some more advanced ways, all explained in depth below:
- From a 1D `std::vector<double>` (`Matrix::fromArray<r,c>(Raw1DArray)`)
- From a 2D `std::vector< std::vector<double> >` (`Matrix(Raw2DArray)`)
- From another matrix (copy constructor, copy assignment constructor, move constructor, move assignment constructor)
- From a matrix previously saved to a file (`.loadFrom(filename)`) or as a string (`.serialize()`)

### Copying and moving matrices
the Matrix class has copy and move constructors, using the following syntaxes:
```cpp
auto myMatrix = otherMatrix;
auto myMatrix = Matrix(otherMatrix);
```
This will move the matrix if possible, otherwise it will perform a deep copy.

#### Copying rows and columns
Use `.copyRow(i)` and `.copyColumn(i)` to copy the i-th row and i-th column respectively, returning a `PH::Raw1DArray` (alias for `std::vector<double>`).

### Accessing elements in a Matrix
Unlike the library on which this one is based, we do not use the subscript (`[]`) operator. In the original library, the subscript operator exposed the underlying data structure, which was column-ordered. Whereas users would expect to access the element at row `r` and column `c` with `matrix[r][c]`, it would be the reverse. As such, the subscript operator has not merely been **changed**, as this would cause migrated code to behave unexpectedly — it has been removed and will return a `PH::NotImplementedException` when used.

This library provides exactly two accessors, not including automatic immutable versions for `const Matrix`:
```cpp
// get the element at row r and column c
myMatrix(r,c);

// get the element at linear index i (reading order)
myMatrix(i);
```

#### Taking and inserting submatrices
Use `.range(beginRow,beginColumn,endRow,endColumn)` to retrieve a rectangular "slice" of the matrix, i.e. a submatrix. Will fail if the dimensions are outside the range of the matrix, if the range goes backward, or if the resulting submatrix is `(0,0)`. **The submatrix is returned as a wholly new matrix rather than a reference.**

Use `.insert(source,beginRow,beginColumn,endRow,endColumn)` to replace the given rectangular submatrix with the given "source" matrix. This will fail for the same reasons as `range`, but additionally will fail if the source matrix does not exactly match the given submatrix area. You can first crop a matrix with `range` to the appropriate size and then `insert` it if you want to "clone" a submatrix from it to another matrix. **This operation overwrites the submatrix in the destination matrix, i.e. it is in-place/mutating.**

### Matrix dimensions
- Use `.isSquare()` to retrieve a boolean value representing if the matrix is square (rows = columns).
- `.size()` returns the number of cells (rows * columns).
- `.rows()` and `.columns()` returns the number of rows and columns, respectively.
- `.dimensions()` returns a `PH::Dimensions` struct, useful for comparing for equality with the same struct from another matrix. This feature is used internally.

#### Resizing arrays
Use `.resize(r,c)` to resize your array to size `(r,c)`. This will destroy elements if shrinking and create new cells if expanding. Value initialization depends on the behavior of `std::vector::resize` for the underlying data structures of the Matrix and may or may not be zero in newly created cells.

### String representation
To return a string representation of the Matrix, use `.str(int precision)`. If no argument is provided, default (full) precision is displayed. This returns an `std::string`, and **not** a `char*` ("C String").

To print a string representation directly, just use `std::cout << myMatrix`. This method uses the above `.str` method with precision `PH::Matrix.displayPrecision`.

### Saving and Loading
This is extremely easy in this library. To save to a format compatible with Numpy etc. (space-delimited columns, newline-delimited rows), use the `saveTo` method:
```cpp
// relative path
myMatrix.saveTo("../data/myMatrix.txt");

// absolute path
myMatrix.saveTo("/Users/paul/Desktop/myMatrix.txt");
```
This will throw a runtime error saying `Could not open file (/bad/path/to/file) to save to.` if you are providing a folder that does not exist or a path to which the program does not have access.

Loading is just as easy, assuming the same format.
```cpp
Matrix myMatrix;
myMatrix.loadFrom("../data/myMatrix.txt");
```

### Convenience initializers
The **identity matrix** is a square (n,n), zeroed-out array with ones in the diagonals. It has been called `eye` here to accomodate the refugees from MATLAB.
```cpp
auto myMatrix = Matrix::eye(4);
std::cout << myMatrix;
```
```
    1    0    0    0
    0    1    0    0
    0    0    1    0
    0    0    0    1
```

This library can generate a **random matrix** with values following a uniform distribution from 0.0...1.0. This is much more powerful than the random generator in MATLAB, which is pseudorandom, but not guaranteed to follow a uniform distribution.
```cpp
auto myMatrix = Matrix::random(3,5);
std::cout << myMatrix;
```
```
    0.0102902    0.261425    0.508335    0.262773    0.53815
    0.494359    0.270736    0.591562    0.222625    0.680818
    0.624379    0.198428    0.5823    0.0862867    0.297198
```

A **linear span** from values `a` to `b` can be generated. It will first generate a linear span across `r` times `c` values (number of cells in the matrix) and will proceed to fill the matrix in reading order, left-to-right, top-to-bottom. The syntax for a **logarithmic span** is identical, just with the `logSpace` function.
```cpp
auto myMatrix = Matrix::linSpace(0,0.8,3,3);
std::cout << myMatrix;
```
```
    0    0.1    0.2
    0.3    0.4    0.5
    0.6    0.7    0.8
```

### Converting STL vectors into a Matrix
#### A one dimensional STL vector
If you have a one dimensional `std::vector`, you can make it into a Matrix with the following code. You can use the `Raw1DArray` type alias to guarantee compatibility with this library even if the numeral type changes from `double`.
```cpp
std::vector<double> myArray = {1,2,3,4,5,6,7,8,9};
auto myMatrix = Matrix::fromArray<3,3>(myArray);
```
This will attempt to create a 3 row, 3 column Matrix and fill in the values from the array in reading order (left-to-right, top-to-bottom). It will throw an exception if the size of the array does not match the numbers in the `fromArray<r,c>` template. Because of the ambiguity as to the size of the desired matrix, there is not a direct-assignment alias for this method (you can't just say `Matrix myMatrix = myArray` when `myArray` is an `std::vector<double>`).

#### A two dimensional STL vector
If you have a 2D `std::vector`, it is even easier because the dimensional data is already provided. I recommend using the `Raw2DArray` type alias I provide to avoid typing too much.
```cpp
Raw2DArray my2DArray = {{1,2,3},{4,5,6},{7,8,9}};

// Method 1
Matrix myMatrix = my2DArray;
// Method 2
auto myMatrix = Matrix(my2DArray);
```
Notice that even though the `std::vector` is constructed with an initializer list, the `Matrix(Raw2DArray)` constructor is not the same as the above initializer list constructor, which lets you provide literal values directly. These methods will fail if the raw array is not consistenly shaped.

### Iterating matrices
Instead of writing cumbersome `for` loops and (even worse) nested `for` loops, the Matrix class provides several useful iterators.
```cpp
Matrix myMatrix = {{1,2,3},
                   {4,5,6},
                   {7,8,9}}
```
#### Iterating columns
There are two ways to iterate columns: one gives you access to the column through a `Raw1DArray` reference (`std::vector<double>`), and the other is immutable (`const`). Both run synchronously, and provide the column number as an parameter in the callback.
```cpp
// Iterate columns while editing them.
// This example sets every diagonal to 1.
myMatrix.mapColumns([&](Raw1DArray& column, Index index) {
	column[index] = 1;
});
``
``
#include <functional>
// Iterate columns without editing them.
// This example sums each column.
auto sums = Raw1DArray(3);
myMatrix.forEachColumn([&](const Raw1DArray& column, Index index) {
	// if you put 0 instead of 0.0, the answer will be an int.
	sums[index] = std::accumulate(column.begin(), column.end(), 0.0);
});
```
You will get a strange error if you get the types wrong in the above lambdas: note the addition of `const` in `forEachColumn`. 

**NOTE:** It's a pain to provide the types for the lambdas above. You can use `auto` instead, as long as you understand what types they are supposed to be:
```cpp
myMatrix.forEachColumn([&](auto column, auto index) {
	// column is a constant reference to an array
	// index is an unsigned integer
	sums[index] = std::accumulate(column.begin(), column.end(), 0.0);
});
```

#### Iterating elements
Iterating each element of the matrix is also easy, and comes in a mutable and immutable version.
```cpp
myMatrix.mapElements([&](auto element, Index row, Index column) {
	// set the cell to the sum of the row number and column number.
	// (why would you do this?)
	element = row + column;
});

// We call it "value" here because it's not a mutable reference,
// just an immutable reference we can retrieve and use.

// Here, we take an elementwise sum of the whole matrix.
// This is not the most creative or impressive use of this function.
auto sum = 0.0;
myMatrix.forEach([&](auto value, Index row, Index column) {
	sum += value;
});
```

#### Matrix arithmetic
Most arithmetic operations are divided into two categories: mutating (in-place) and non-mutating. Mutating methods will act upon the first matrix in the operation, whereas non-mutating methods will return a new matrix as a result of the operation. Where the difference is not made clear with just non-Latin operators (`+`,`-`,`/`,`*`, etc.), an explicit English name is used to avoid ambiguity.

Below, type names will be used where you would expect variable names for concision and clarity. All elementwise operations involving two matrices require them to be the *same size*. All mutating operations return a reference to the matrix being mutated (allowing for operation chaining), and all non-mutating operations return a new matrix (or a Vector or scalar, depending on the operation).

| Operation | Mutating | Syntax | Description |
|--------------------------------|----------|----------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `x*A+y*B` | Yes | `Matrix.linearSumInPlace (MathNumber,MathNumber,Matrix)` | First, multiplies each element of this matrix by the first number. Then adds each element of the second matrix to the corresponding element in this matrix, scaled by the second number. |
| **Elementwise Addition** | Yes | `Matrix += MathNumber` | Adds the number to each element of this matrix. |
| **Elementwise Sum** | Yes | `Matrix += Matrix` | Adds each element of the right matrix to each corresponding element of the left matrix. Only the left matrix is mutated. |
| **Elementwise Subtraction** | Yes | `Matrix -= MathNumber` | Subtracts the number from each element of this matrix. |
| | No | `Matrix - MathNumber` or `MathNumber - Matrix` | Returns a new matrix where each element has had the number subtracted, or where each element has been subtracted from the number (order-dependent). |
| **Elementwise Difference** | Yes | `Matrix -= Matrix` | Subtracts each element of the right matrix from each corresponding element of the left matrix. Only the left matrix is mutated. |
| | No | `Matrix - Matrix` | Returns a new matrix where each element is the difference between the corresponding elements in the two given matrices. |
| **Elementwise Product** | Yes | `Matrix.elementwiseMultiply` `(Matrix)` | *Not to be confused with the dot product or matrix multiplication*. Multiplies each element of this matrix by the corresponding element of the given matrix. |
| **Elementwise Multiplication** | Yes | `Matrix *= MathNumber` | Scales each element of this matrix by the number. |
| **Elementwise Division** | Yes | `Matrix.elementwiseDivide` `(Matrix)` | *Not to be confused with right or left division*. Divides each element of this matrix by the corresponding element of the given matrix. |
| **Elementwise Quotient** | Yes | `Matrix /= MathNumber` | Divides each element of this matrix by the number. |
| **Elementwise Power** | Yes | `Matrix.elementwisePower` `(MathNumber)` | *Not to be confused with actual matrix exponentials, which require square matrices.* Takes the power of each element of this matrix to the given number. |
|  | Yes | `Matrix.elementwisePower` `(Matrix)` | Takes the power of each element of this matrix to the corresponding element of the given matrix. |
| **Value Fill** | Yes | `Matrix = MathNumber` | Fills every element of the matrix with the number.|
| **Absolute Value** | Yes | `Matrix.absInPlace()` | Rewrites the matrix with absolute values.|
|  | No | `Matrix.abs()` | Returns a new matrix with the absolute values of each element of this matrix.|
| **Matrix Transpose** | Yes | `Matrix.transposeInPlace()` | *Only for square matrices.* Transposes the matrix in place.|
| | No | `Matrix.transpose()` | *Non-square matrices must use this.* Returns a new matrix that is the transposed `(c,r)` version of this `(r,c)` matrix.|
| **Matrix Inverse** | No | `Matrix.inverse()` | *Fails if the matrix is square or singular.* Returns a new matrix that is the inverse of this one.|
| **Backward Substitution** | No | `Matrix::backSubstitution` `(Matrix,Matrix)` | (static, literally called `Matrix::backSubstitution`). Given two matrices U and B, performs backward substitution on the linear system `U*X=B`, returning X as a new matrix, without mutating U or B. |
| | No | `Matrix::backSubstitution` `(Matrix,Vector)` | Given a matrix U and vector b, performs backward substitution on the linear system `U*x=b`, returning x as a new vector, without mutating U or b. |
| **Forward Substitution** | No | `Matrix::forwardSubstitution` `(Matrix,Matrix)` | (static, literally called `Matrix::forwardSubstitution`). Given two matrices L and B, performs forward substitution on the linear system `L*X=B`, returning X as a new matrix, without mutating L or B. |
| | No | `Matrix::forwardSubstitution` `(Matrix,Vector)` | Given a matrix L and vector b, performs backward substitution on the linear system `L*x=b`, returning x as a new vector, without mutating L or b. |
| **Linear Solution** | No | `Matrix::linearSolve` `(Matrix,Matrix)` | Given two matrices A and B, solves the linear system `A*X=B`, returning the matrix X holding the result. **A and B are modified in the course of this operation.**|
| | No | `Matrix::linearSolve` `(Matrix,Vector)` | Given a matrix A and vector b, solves the linear system `A*x=b`, returning the vector x holding the result. **A and b are modified in the course of this operation.**|
| **Dot Product** | No | `Matrix::dot(Matrix,Matrix)` | (note that this function is static, and is actually called verbatim as `Matrix::dot`) Takes the dot product of the two matrices. Returns a scalar. |
| **Matrix Multiplication** | No | `Matrix * Matrix` | Returns a new matrix as the result of matrix multiplication. *This is not the dot product or elementwise multiplication.*|
| **Matrix-Vector Product** | No | `Matrix * Vector` | Returns a new vector as the result of matrix-vector multiplication.|
| **Matrix Comparison** | No | `Matrix == Matrix` | Not as simple as it seems. Floating point numbers cannot be exactly compared, so the matrices are said to be equal if they are (1) the same size, (2) each pair of values have the same exponent, and (3) each pair's significands are within a tolerance of 1e-6.|

#### Scalar methods
| Name | Syntax | Description |
|-----------|--------|-------------|
|**Minimum**|`Matrix.min()`|Returns the minimum value in the matrix.|
|**Maximum**|`Matrix.max()`|Returns the maximum value in the matrix.|
|**Frobenius Norm**|`Matrix::norm(Matrix)`|Returns the scalar result of the matrix Frobenius norm. Acts as a 2-norm on column or row vectors represented as matrices.|
|**Infinity Norm**|`Matrix::infNorm(Matrix)`|Returns the scalar result of the matrix infinity norm. Acts as an infinity norm on column vectors, and a one-norm on row vectors represented as matrices.|
|**One Norm**|`Matrix::oneNorm(Matrix)`|Returns the scalar result of the matrix one norm. Acts as a one norm on column vectors, and an infinity norm on row vectors represented as matrices.|
