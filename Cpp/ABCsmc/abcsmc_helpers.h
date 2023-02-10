#pragma once

void propose(std::vector<float>& thetaStar) {
	/*
	 Define some proposal distributions to draw from and store the
	 draws in the pre-defined thetaStar (std::vector)
	 */
}

void perturb(std::vector<double>& theta){
	theta[0] += unif01(generator)*0.01 - 0.005 ;
	theta[1] += unif01(generator)*0.001 - 0.0005;
	theta[2] += unif01(generator)*0.01 - 0.005 ;
}

const double priorNorm = 2.0*2.0*0.1 ;
double priorDens(std::vector<double>& theta){
	int inRange = theta[0]>0.0 & theta[0]<2.0 & theta[1]>0.0 & theta[1]<0.1 & theta[2]>0.0 & theta[2]<2.0;
	return inRange / priorNorm ;
}
 
const double transNorm =  0.1*0.1*0.005 ;
double transDens(std::vector<double>& theta_old, std::vector<double> theta_new){
	int inRange = theta_new[0]>(theta_old[0]-0.005)  & theta_new[0]<(theta_old[0]+0.005) &            		  theta_new[1]>(theta_old[1]-0.0005) & theta_new[1]<(theta_old[1]+0.0005) & 		  		  theta_new[2]>(theta_old[2]-0.005)  & theta_new[2]<(theta_old[2]+0.005);
	return inRange / transNorm ;
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

double weightCalc(std::vector<double>& theta_new, std::vector<std::vector<double>>& theta_old, std::vector<double>& weights_old){
	
	double btm = 0.0;
	for(std::size_t i=0; i<theta_old.size(); ++i)
		btm += weights_old[i]*transDens(theta_old[i], theta_new);
	
	return priorDens(theta_new) / btm ;
}

int weightedSample(std::vector<double> weights){
	double norm = 0.0 ;
	for(std::size_t i=0; i<weights.size(); ++i)
		norm += weights[i];
	
	double rnd = unif01(generator)*norm;
	
	for(unsigned int i=0; i<weights.size(); ++i){
		rnd -= weights[i];
		if( rnd < 0 )
			return i ;
	}
	return -1;
}
