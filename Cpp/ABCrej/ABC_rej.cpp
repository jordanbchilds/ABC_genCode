#include<cmath>
#include<iostream>
#include<fstream>
#include<random>
#include<algorithm>
#include "abc_helpers.h"
#include "helpers.h"
#include "simulator.h"

const float epsilon ;
// The threshold of acceptance
const unsigned int nOut ;
// nOut: the number of approximate posterior samples we wish to have

int main(){
	// read the data in as a single vector.
	unsigned int count = 0;
	std::vector<float> rawData;
	std::fstream dataFile ("data.txt", std::ios_base::in);
	float a;
	while (myfile >> a){
		rawData.push_back(a);
		++count;
	}
	dataFile.close();
	
	float myData[nTimes][nObs] ;
	// Do some manipulation to put the data into a two dimension array (or form you want)
	/*
		do stuff to myData
	*/
	
	std::vector<float> dataSumm = mySummary(myData);
	// dataSumm: summary of data used for ABC algorithm
	// Does not have to be std::vector
	
	/*
		define dataSumm using a function in abc_helpers e.g. quantiles
	*/
	
	myData.clear();
	// After defining dataSumm the raw data is not needed.

	const unsigned int nParam ;
	// nParam: the number of parameters we're trying to learn
	const unsigned int nSim ;
	// nSim: the numbe of simulations of the system done to compare to the data
	const int nTime ;
	// nTime: number of data points
	// assume each time point has the same number of observations

	const float* simOutput = new float[nSim][nTime];
	// simOutput: an array to store the output of the simulation.
	// This does not need to be kept from one iteration of the algorithm to the next.
	// Storing on heap because sometimes memory issues occured if this was large


	vector<float> thetaStar(nParam);
	// thetaStar: proposed parameter values
	
	vector<vector<float>> theta(nOut);
	// theta: array to store accepted parameter values
	const float* theta_ptr;

	unsigned int count = 0;
	// Number of accepted approximate posterior draws

	while( count < nOut ){
		propose_theta(thetaStar);
		// Propose parameter values.
		// proposal is defined in abc_helpers.h
		//
		for(int i=0; i<nSim; ++i){
			// simulate the system using function written in simulator.h
			// save the output to our simOutput array on the heap
			/*
			 simulate(simOutput[i], thetaSter);
			 */
			// Define the simulation function in simulator.h
		}
		
		vector<float> dataStar_summ = mySummary(simOutput);
		// calculate the summary statistics of the simulated data
		
		if( myDist(dataSumm, dataStar_summ) < epsilon ){
			theta[count] = thetaStar ;
			count++;
		}
		
	}
	
	delete [] simOutput;
	// delete the simulated data array, to prevent any weird leaks n stuff
	
	// Save output, accepted parameter values in a text file using a tab seperator.
	std::fstream outFile ("posterior.txt");
		for(int i=0; i<nOut; i++){
			for(int j=0; j<nParam-1; j++){
				outFile << theta[i][j] << "\t";
			}
			outFile << theta[i][nParam-1] << "\n";
		}
	dataFile.close();

}
