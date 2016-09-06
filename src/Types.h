//
//  Types.h
//  matrix
//
//  Created by Paul Herz on 8/30/16.
//  Copyright © 2016 Paul Herz. All rights reserved.
//

#ifndef Types_h
#define Types_h

// for accessing pi via M_PI
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <cmath>
#else
#include <cmath>
#endif

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>
#include <random>

// singularity tolerance
#define STOL 1.e-15

// the precision levels for serializing vs displaying
// numbers in Matrices and Vectors
static int fullPrecision = std::numeric_limits<double>::digits10 + 1;
static int displayPrecision = 6;

using  MathNumber = double;
using  Raw1DArray = std::vector<MathNumber>;
using  Raw2DArray = std::vector<Raw1DArray>;
using       Index = std::size_t;

#endif /* Types_h */
