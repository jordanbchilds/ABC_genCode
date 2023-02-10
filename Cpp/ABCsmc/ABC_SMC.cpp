#include<cmath>
#include<iostream>
#include<fstream>
#include<random>
#include<algorithm>
#include<thread>
#include "helpers.h"

using namespace std;

const unsigned int nCores = 5;

const unsigned int nSpecies = 2;
const unsigned int nReactions = 3;
const double PRE[nReactions][nSpecies] = { {1,0}, {1,1}, {0,1} };
const double POST[nReactions][nSpecies] = { {2,0}, {0,2}, {0,0} };
const double STOI[nReactions][nSpecies] = { {1,0}, {-1,1}, {0,-1} };

const unsigned int nTimes = 20;
double outTimes[nTimes];
const double* outTimes_ptr = &outTimes[0];

const double xInit[nSpecies] = {50,100};
const unsigned int nParam = 3;

#include "gillespie.h"
#include "abc_helpers.h"

const vector<double> epsilon = {5000.0, 3000.0, 1500.0, 1000.0, 800.0, 700.0, 650.0};

int main () {
	
	int timeCount = 0;
	std::fstream timeFile("synTime.txt", std::ios_base::in);
	double t;
	while (timeFile >> t){
		outTimes[timeCount] = t;
		++timeCount;
	}
	timeFile.close();
	
	vector<double> dataSumm(nTimes*nSpecies) ;
	int dataCount = 0;
	std::fstream dataFile("synData.txt", std::ios_base::in);
	double x;
	while (dataFile >> x){
		dataSumm[dataCount] = x;
		++dataCount;
	}
	timeFile.close();
	
	const unsigned int nOut = 1e2;
	const unsigned int nSim = 10;
	
	vector<vector<double>> theta(nOut);
	vector<double> weights(nOut);
	for( size_t i=0; i<nOut; ++i ){
		theta[i] = propose();
		weights[i] = 1.0;
	}
	
	for(double ep : epsilon ){
		vector<vector<double>> theta_new(nOut);
		vector<double> weights_new(nOut);
		unsigned int nAccept = 0;
		
		while( nAccept < nOut ){
			vector<double> theta_star = theta[weightedSample(weights)];
			perturb(theta_star);

			if( priorDens(theta_star)>0.0 ){
				vector<vector<double>> rawData_star(nSim);
				for(size_t i=0; i<rawData_star.size(); ++i)
					rawData_star[i] = vector<double>(nTimes*nSpecies);
				
				vector<double> dataSumm_star(nTimes*nSpecies);
				
				for(size_t i=0; i<nSim; ++i)
					gillespied(rawData_star[i].begin(), rawData_star[i].end(), theta_star);
				
				colMeans(dataSumm_star, rawData_star);
				double dist = myDist(dataSumm, dataSumm_star);
				if( dist < ep ){
					cout << "ACCEPTED: " << nAccept+1 <<endl;
					theta_new[nAccept] = theta_star;
					weights_new[nAccept] = weightCalc(theta_star, theta, weights);
					nAccept++;
					
					cout << ep << ": ";
					for(size_t i=0; i<nParam; ++i)
						 cout << theta_star[i] << " ";
					cout << endl;
				}
			}
		}
		for(size_t i=0; i<nOut; ++i){
			weights[i] = weights_new[i];
			theta[i] = theta_new[i];
		}
		
		weights_new.clear();
		theta_new.clear();
	}
	
	std::ofstream outFile;
	outFile.open ("ABCSMC_output.txt");
	for(size_t i=0; i<nOut; ++i){
		for(size_t j=0; j<nParam; ++j)
			outFile << theta[i][j] << " ";
		outFile << "\n";
	}
	outFile.close();
	 
	return 0;
}
