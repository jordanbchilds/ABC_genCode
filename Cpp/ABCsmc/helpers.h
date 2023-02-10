#pragma once

auto getTime() { return std::chrono::high_resolution_clock::now(); };

std::default_random_engine generator;

std::uniform_real_distribution<double> unif01(0.0, 1.0);

template<typename T>
static inline double Lerp(T v0, T v1, T t){
	return (1 - t)*v0 + t*v1;
}

template<typename T>
static inline std::vector<T> Quantile(const std::vector<T>& inData, const std::vector<T>& probs, std::vector<T>& output){
	if ( inData.empty() ){
		return std::vector<T>();
	}
	if ( 1 == inData.size() ){
		return std::vector<T>(1, inData[0]);
	}

	std::vector<T> data = inData;
	std::sort(data.begin(), data.end());
	std::vector<T> quantiles;

	for (size_t i = 0; i < probs.size(); ++i){
		T poi = Lerp<T>(-0.5, data.size() - 0.5, probs[i]);

		size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
		size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));

		T datLeft = data.at(left);
		T datRight = data.at(right);

		T quantile = Lerp<T>(datLeft, datRight, poi - left);

		output[i] = quantile;
	}

	return output;
}

template<typename T>
static inline std::vector<std::vector<T>> colQuantile(const std::vector<std::vector<T>>& inData, const std::vector<T>& probs, std::vector<std::vector<T>>& output){
	if ( inData.empty() ){
		throw std::invalid_argument("colQuantiles requires non-empty input data of type std::vector<std::vector<T>>");
	}
	if ( 1 == inData.size() ){
		throw std::invalid_argument("colQuantile requires multi-dimensional datam Quantile can used for one-dimension data.");
	}

	for(std::size_t t=0; t<inData[0].size(); ++t){
		std::vector<T> colData(inData.size());
		for(std::size_t i=0; i<inData.size(); ++i){
			colData[i] = inData[i][t];
		}
		
		std::sort(colData.begin(), colData.end());
		std::vector<T> colQuantiles;
		
		for (size_t i = 0; i < probs.size(); ++i){
			T poi = Lerp<T>(-0.5, colData.size() - 0.5, probs[i]);

			size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
			size_t right = std::min(int64_t(std::ceil(poi)), int64_t(colData.size() - 1));

			T datLeft = colData.at(left);
			T datRight = colData.at(right);

			T quantile = Lerp<T>(datLeft, datRight, poi - left);
			output[t][i] = quantile;
		}
	}
	return output;
}

template<typename T>
float euclidean_dist(std::vector<T>& A, std::vector<T>& B){
	if( A.empty() | B.empty() ){
		throw std::invalid_argument("The matrices must be non-empty and of type std::vector<std::vector<T>>");
	}
	if( A.size() != B.size() ){
		throw std::invalid_argument("The matrices must be the same dimension");
	}
	
	float dist = 0.0;
	
	for(std::size_t i=0; i<A.size(); ++i){
		dist += sqrt(pow(A[i] - B[i], 2));
	}
	return dist;
}

template<typename T>
float euclidean_dist(std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& B){
	if( A.empty() | B.empty() ){
		throw std::invalid_argument("The matrices must be non-empty and of type std::vector<std::vector<T>>");
	}
	if( A.size() != B.size() ){
		throw std::invalid_argument("The matrices must be the same dimension");
	}
	for(std::size_t i=0; i<A.size(); ++i){
		if( A[i].size() != B[i].size() ){
			throw std::invalid_argument("The matrices must be the same dimension");
		}
	}
	
	float dist = 0.0;
	
	for(std::size_t i=0; i<A.size(); ++i){
		float sqSum = 0.0;
		for(std::size_t j=0; j<(A[0]).size(); ++j)
			sqSum += pow(A[i][j] - B[i][j], 2);
		
		dist += sqrt(sqSum);
	}
	return dist;
}

void colMeans(std::vector<double>& outMean, std::vector<std::vector<double>>& matrix){
	
	for(std::size_t i=0; i<outMean.size(); ++i){
		double sum = 0.0;
		for(std::size_t j=0; j<matrix.size(); ++j)
			sum += matrix[j][i];
		
		outMean[i] = sum / double(matrix.size());
	}
}
