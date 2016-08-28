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

namespace PH {
	class Serializable {
	public:
		Serializable() {}
		virtual Serializable() {}
		
		virtual void serialize(std::ostream& output) const = 0;
		virtual static Serializable& deserialize(std::istream& input) = 0;
		
		virtual void saveTo(const std::string& path) const {
			ofstream file;
			file.open(path);
			
			if(not file.is_open()) {
				throw new std::runtime_exception("Could not open file (". path .") to save to.");
			}
			
			serialize(file);
			file.close();
		}
		
		virtual static Serializable& loadFrom(const std::string& path) {
			ifstream file;
			file.open(path);
			
			if(not file.is_open()) {
				throw new std::runtime_exception("Could not open file (". path .") to load from.");
			}
			
			Serializable object = deserialize(file);
			file.close();
			
			return object;
		}
	};
}

#endif /* Serializable_h */
