# Matrix
by Paul Herz

I developed this library off of Dr. Reynolds' implementation. It has been approved as a drop-in replacement for that library in the course of the High Performance Scientific Computing class.

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

### Copying and moving matrices:
the Matrix class has copy and move constructors, using the following syntaxes:
```cpp
auto myMatrix = otherMatrix;
auto myMatrix = Matrix(otherMatrix);
```
This will move the matrix if possible, otherwise it will perform a deep copy.

### Converting `std::vector` into a Matrix:
#### A one dimensional `std::vector<double>`
If you have a one dimensional `std::vector`, you can make it into a Matrix with the following code. You can use the `Raw1DArray` type alias to guarantee compatibility with this library even if the numeral type changes from `double`.
```cpp
std::vector<double> myArray = {1,2,3,4,5,6,7,8,9};
auto myMatrix = Matrix::fromArray<3,3>(myArray);
```
This will attempt to create a 3 row, 3 column Matrix and fill in the values from the array in reading order (left-to-right, top-to-bottom). It will throw an exception if the size of the array does not match the numbers in the `fromArray<r,c>` template. Because of the ambiguity as to the size of the desired matrix, there is not a direct-assignment alias for this method (you can't just say `Matrix myMatrix = myArray` when `myArray` is an `std::vector<double>`).

#### A two dimensional `std::vector< std::vector<double> >`
If you have a 2D `std::vector`, it is even easier because the dimensional data is already provided. I recommend using the `Raw2DArray` type alias I provide to avoid typing too much.
```cpp
Raw2DArray my2DArray = {{1,2,3},{4,5,6},{7,8,9}};

// Method 1
Matrix myMatrix = my2DArray;
// Method 2
auto myMatrix = Matrix(my2DArray);
```
Notice that even though the `std::vector` is constructed with an initializer list, the `Matrix(Raw2DArray)` constructor is not the same as the above initializer list constructor, which lets you provide literal values directly. These methods will fail if the raw array is not consistenly shaped.

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

### Serializing and deserializing (saving to and loading from a file)
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

## Adding to your project
To use the included `Matrix` and `Vector` classes in your project, simply drop the contents of the `src` folder here into your project. You may want to separate this code into a `lib` folder to denote that it is not your project code if using it in HPSC.

### Including into a traditional project (GCC, Makefiles)
1. Copy the contents of `src` in this repo, including all headers, into a location in your project.
2. Add `Matrix.cpp` and `Vector.cpp` to your build list (e.g. `g++ yourfile.cpp otherfile.cpp Matrix.cpp Vector.cpp`, including the relative path if they are in a different folder than your code.
3. If you've placed this code in a different folder than where you are running the compiler, add the include header folder flag (`-I /lib` if the files are in `lib`).
4. You will need to use at least the C++11 standard, specifically the GNU extension, using the flag `-std=gnu++11`. Your project may not work if it is using an older standard or a non-GNU extension standard. If you have trouble with a *newer* standard, file an issue and I will respond promptly.

### Including into an Xcode Project (Xcode 9.0)
This example starts from absolute scratch in Xcode 9.0. This version changed the process of including files to be a little more complicated and restrictive. As such, it is difficult to actually store this code in a separate folder from your own. However, we can still use Xcode's Groups (pseudo-folders in the Xcode sidebar) to *visually* group them apart from your own code.

1. Open Xcode to the welcome screen, and click **Create a New Xcode Project**. If you don't get a welcome screen when you open Xcode, go to `File > New > Project...` in the toolbar, or press `Shift + Cmd + N`.
2. Select **macOS** from the tabs at the top of the New Project wizard, then select **Command Line Tool** under the **Application** category (at the top of the list). Click **Next**.
3. Enter a product name, preferably without spaces. Select **C++** as your language. Click **Next**.
4. Pick a location to save your project to. Be wary of the **Create Git repository on my Mac**. You may or may not want this checked. Xcode will create a folder named after your project here, so there's no need to do so manually.
5. Open up a Finder window next to Xcode. In that window, browse to where you've downloaded/cloned this code, specifically to the `src` folder. Select all of the *contents* of this folder (not the folder itself) and *drag* the files into Xcode, dropping them directly above or below the `main.cpp` file. 
6. In the popup, select **Copy items if needed** and **Create groups**. Click **Finish**.
7. These files may clutter your future or existing code. You can split them off into a group on the sidebar, but it will not change their actual location in the folder. Select all the files belonging to the Matrix library by clicking the first file, holding `Shift`, and clicking the last. Right click, and select **New Group from Selection**. Name it something like `lib` or `matrix` to identify it as separate. Now all of the Matrix library files will be bundled into a folder and separated. If you don't want this group *inside* of your `ProjectName` group, you can carefully drag it above the `ProjectName` group on the sidebar in Xcode and the group as well as the folder will be moved to a *sibling* of the `ProjectName` group.
8. If you look in Finder and it appears that Xcode has *copied*, rather than *moved* the Matrix files in your project (there is one copy of all the files inside /ProjectName/ProjectName and another in /ProjectName/ProjectName/GroupName), you will need to manually delete the old copies (in /ProjectName/ProjectName). This is a behavioral issue in Xcode 9.0.
9. Select your project from the file sidebar in Xcode (it's the blueprint icon, at the top usually). Select the **Build Settings** tab. Search `C++ Language Dialect` in the search bar below the tab bar.
10. In the resulting setting that appears, make sure the selection is *at least* **GNU++11 [-std=gnu++11]**. This library depends on at least C++11, and has been tested with only the C++11 standard with GNU extensions. Further testing is to be done.
11. Select the **Build Phases** tab now under the same project, and open the **Compile Sources** dropdown by clicking on the arrow next to it. Make sure `Matrix.cpp` and `Vector.cpp` are added. If they are not, click the plus and add them.
12. You should be able to add the following to your `main.cpp` file or any other file in the project:
```
#include "Matrix.h"
```
Xcode should be able to find the header files.
