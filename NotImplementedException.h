//
//  NotImplementedException.cpp
//  matrix
//
//  Created by Paul Herz on 8/28/16.
//  Copyright © 2016 Paul Herz. All rights reserved.
//

#include <exception>

class NotImplementedException: std::logic_error {
public:
	NotImplementedException():
		std::logic_error("This function is not implemented.") {}
};
