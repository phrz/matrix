//
//  Serializable.h
//  matrix
//
//  Created by Paul Herz on 8/28/16.
//  Copyright © 2016 Paul Herz. All rights reserved.
//

#ifndef Serializable_h
#define Serializable_h

#include <iostream>

class Serializable {
public:
	Serializable() {}
	virtual Serializable() {}
	
	virtual void serialize(std::ostream& output) = 0;
	virtual void deserialize(std::istream& input) = 0;
};

#endif /* Serializable_h */
