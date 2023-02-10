#pragma once

void propose(std::vector<float>& thetaStar) {
	/*
	 Define some proposal distributions to draw from and store the
	 draws in the pre-defined thetaStar (std::vector)
	 */
}

double myDist( std::vector<double>& vecA, std::vector<double>& vecB ){
	if( vecA.size() != vecB.size() )
		throw std::invalid_argument("The vectors must be the same length");
	
	double sum = 0.0;
	for(std::size_t i=0; i<nTimes; ++i){
		double xx = 0.0;
		for(std::size_t j=0; j<nSpecies; ++j)
			xx += pow(vecA[i*nSpecies + j] - vecB[i*nSpecies + j], 2);
		
		sum += sqrt(xx);
	}
	return sum ;
}

