#pragma once


unsigned choose(unsigned n, unsigned k){
	if (k>n) return 0;
	if (k*2>n) k = n-k;
	if (k==0) return 1;
	if (k==1) return n;
	int result = n;
	for(int i=2; i<=(n-k); ++i){
		result *= (n-i+1);
		result /= i;
	}
	return result;
}

float rand_unif(float lower=0.0, float upper=1.0){
	float unif_01 = (float) rand() / (float) RAND_MAX ;
	return unif_01*(upper-lower) + lower;
}

float runif_01(){
	return (float) rand() / (float) RAND_MAX ;
}

double rand_exp(float lambda){
	return -1.0*log(1.0-runif_01())/lambda;
}

int randReact(float* weights, unsigned int n){
	double norm = 0.0 ;
	for(std::size_t i=0; i<n; ++i)
		norm += *(weights+i) ;
	
	double rnd = unif01(generator)*norm;
	
	for(unsigned int i=0; i<n; ++i){
		rnd -= *(weights + i) ;
		if( rnd < 0 )
			return i ;
	}
	return -1;
}

void gillespied(std::vector<double>::iterator begin, std::vector<double>::iterator end, std::vector<double>& rates){

	float Tmax = *(outTimes_ptr + nTimes - 1);
	
	std::vector<double>::iterator itr = begin;

	double xx[nSpecies];
	xx[0] = xInit[0]; xx[1] = xInit[1];
	
	int count = 0;
	float tt = 0.0;

	while( tt <= Tmax ){
		float hazards[nReactions];
		float haz_total = 0.0;
		for(int i=0; i<nReactions; ++i){
			float h_i = rates[i];
			for(int j=0; j<nSpecies; ++j)
				h_i *= choose(xx[j], PRE[i][j]);
			hazards[i] = h_i;
			haz_total += h_i;
		}
		
		tt += rand_exp(haz_total);
		if( tt >= *(outTimes_ptr + count) && itr!=end){
			for(std::size_t j=0; j<nSpecies; ++j){
				*itr = xx[j];
				itr++;
			}
			count += 1;
		}
		
		int r = rand_react(hazards, nReactions);
		for(std::size_t j=0; j<nSpecies; ++j)
			xx[j]  += STOI[r][j];
		
		if( haz_total <= 1e-10 ){
			while( itr != end){
				for(std::size_t j=0; j<nSpecies; ++j){
					*itr = xx[j];
					itr++;
				}
			}
			tt = 1e99;
		}
	}
}
