//	// create a new vector with uniformly-distributed random numbers in [0,1]
//	Vector randomVectorOfSize(Index n) {
//		if (n<1) {
//			throw new std::invalid_argument("randomVectorOfSize expects a vector size argument > 0.");
//		}
//		
//		// Create random device for a uniform integer
//		// distribution.
//		
//		// This is much like MATLAB's rand().
//		std::random_device randomDevice;
//		std::mt19937 generator(randomDevice());
//		std::uniform_real_distribution<> dist(0, 1);
//		
//		Vector result(n);
//		result.mapElements([&dist, &generator](MathNumber& element, Index i){
//			element = dist(generator);
//		});
//		
//		return result;
//	}
//
