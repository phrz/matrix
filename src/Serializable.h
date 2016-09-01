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
#include <iomanip>
#include <exception>
#include <fstream>

namespace PH {
	
	template<class SerialType>
	class Serializable {
	public:
		Serializable() {}
		virtual ~Serializable() {}
		
		virtual void serialize(std::ostream& output) const = 0;
		virtual SerialType deserialize(std::istream& input) = 0;
		
		virtual void saveTo(const std::string& path) const {
			std::ofstream file;
			file.open(path);
			
			if(not file.is_open()) {
				throw new std::runtime_error("Could not open file (" + path + ") to save to.");
			}
			
			serialize(file);
			// ofstream self-closes, guaranteed.
		}
		
		virtual SerialType loadFrom(const std::string& path) {
			std::ifstream file;
			file.open(path);
			
			if(not file.is_open()) {
				throw new std::runtime_error("Could not open file (" + path + ") to load from.");
			}
			
			return this->deserialize(file);
			// ifstream self-closes, guaranteed.
		}
	};
}

#endif /* Serializable_h */
