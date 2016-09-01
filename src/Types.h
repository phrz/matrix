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

#include <vector>

using  MathNumber = double;
using  Raw1DArray = std::vector<MathNumber>;
using  Raw2DArray = std::vector<Raw1DArray>;
using       Index = std::size_t;

#endif /* Types_h */
