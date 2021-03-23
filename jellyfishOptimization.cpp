#include <iostream>
#include <vector>
#include <fstream>
#include <tuple>
#include <chrono>
#include <functional>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/lognormal_distribution.hpp>

#include "JellyfishClass.h"
#include "pVavInteractions.h"

using namespace std;

double cControl(double curTime, int maxIterations, double RNG){
	return fabs((1-curTime/maxIterations));
}
vector<tuple<double,double> > pVavSetBounds(vector<double>& inSpeciesCounts,vector<double>& inTimes);
int loadPvavInputs(vector<double>& speciesVector, vector<double>& initParameters, vector<double>&stoppingTimes, vector<tuple<double,double> >& bounds, string inFile);
void loadCovariance(vector<double>& outMeans, vector<vector<double> >& inMatrix, string dataPath);



int main(){
	//Initialize RNG generators
	//Main uniform generator for runtime
	boost::mt19937 generator;
	generator.seed(time(NULL));
	boost::uniform_real<> standardUniform(0,1);
	boost::variate_generator<boost::mt19937, boost::uniform_real<double> > mainUniform(generator,standardUniform);

	//Uniform generator for jellyfish motion
	boost::mt19937 jellyGenerator;
	jellyGenerator.seed(time(NULL));
	boost::variate_generator<boost::mt19937, boost::uniform_real<double> > jellyfishRNG(jellyGenerator,standardUniform);



	const int numParameters(6);
	const int numJellyfish(25);
	const int maxIterations(100);
	const int numSpecies(6);
	const double deltaT(0.2);
	const int numSamples(200);

	const int numRuns(100);


	vector<double> speciesVector(numSpecies,0);
	vector<double> resetVector=speciesVector;

	vector<double> initParameters={0.008,0.1,1.0,0.013904,0.05,0.07};
	vector<double> stoppingTimes={0,100};

	vector<tuple<double,double> > parameterBounds=pVavSetBounds(speciesVector,stoppingTimes);


	vector<double (*)(Jellyfish&,vector<double>&)> interactionPointers={dSyk,dVav,dSV,dpVav,dSHP1,dSHP1Vav};

	loadPvavInputs(speciesVector,initParameters,stoppingTimes,parameterBounds,"inputOne.txt");
	resetVector=speciesVector;

	vector<boost::mt19937> exNoiseEngines(3);
	vector<boost::variate_generator<boost::mt19937, boost::lognormal_distribution<double> > > extrinsicNoiseGenerators;
	vector<double> means;
	vector<vector<double> > inValues;
	loadCovariance(means, inValues,"outCov");
	for(int i=0;i<(int)exNoiseEngines.size();i++){
		exNoiseEngines[i].seed(time(NULL));
		double scaledMean(means[i]);
		double scaledSD(inValues[i][i]);
		boost::lognormal_distribution<> currentTest(scaledMean,scaledSD);
		boost::variate_generator<boost::mt19937,boost::lognormal_distribution<double> > createdEngine(exNoiseEngines[i],currentTest);
		extrinsicNoiseGenerators.push_back(createdEngine);
	}

	const int numTimepoints(stoppingTimes.size()-1);

	Jellyfish trueJellyfish=Jellyfish(initParameters,parameterBounds);
	trueJellyfish.interactionFunctions=interactionPointers;

	

	vector<vector<double> > trueData(numTimepoints,vector<double> (numSpecies,0));
	for(int stopTime=0;stopTime<numTimepoints;stopTime++){
		pVav_RungeKutta(trueJellyfish,speciesVector,stoppingTimes[stopTime],stoppingTimes[stopTime+1],deltaT);
		trueData[stopTime]=speciesVector;
	}
	bool extrinsicNoise(true);
	string outFileName;
	if(extrinsicNoise){
		outFileName="extrinsicNoiseResults.txt";
	}
	else{
		outFileName="deterministicResults.txt";
	}
	ofstream outFile(outFileName);
	outFile<<0<<" ";
	for(int i=0;i<(int)trueJellyfish.position.size();i++){
		outFile<<trueJellyfish.position[i]<<" ";
	}
	outFile<<endl;

	for(int curRun=0;curRun<numRuns;curRun++){

		vector<Jellyfish> jellyfishVector(numJellyfish);
		for(int jellyfish=0;jellyfish<numJellyfish;jellyfish++){
			vector<double> position(numParameters);
			for(int i=0;i<(int)position.size();i++){
				auto[lower,upper]=parameterBounds[i];
				position[i]=lower+mainUniform()*(upper-lower);
			}
			jellyfishVector[jellyfish]=Jellyfish(position,parameterBounds);
			jellyfishVector[jellyfish].interactionFunctions=interactionPointers;
		}


		int minIndex(0);
		for(int curTime=0;curTime<maxIterations;curTime++){

			
			double bestFitnessValue(1e13);
			vector<double> meanPosition(numParameters,0);
			vector<vector<double> > testData(numTimepoints,vector<double> (numSpecies,0));
			for(int i=0;i<(int)jellyfishVector.size();i++){
				if(extrinsicNoise){
					for(int sample=0;sample<numSamples;sample++){
						speciesVector=resetVector;
						speciesVector[0]=extrinsicNoiseGenerators[0]();
						speciesVector[1]=extrinsicNoiseGenerators[1]();
						speciesVector[4]=extrinsicNoiseGenerators[2]();
						for(int stopTime=0;stopTime<numTimepoints;stopTime++){
							pVav_RungeKutta(jellyfishVector[i],speciesVector,stoppingTimes[stopTime],stoppingTimes[stopTime+1],deltaT);
							for(int element=0;element<(int)speciesVector.size();element++){
								testData[stopTime][element]+=speciesVector[element]/numSamples;
							}
						}
					}
					jellyfishVector[i].currentFitnessValue=meanNoNoiseFitness(trueData,testData);
					for(int j=0;j<(int)jellyfishVector[i].position.size();j++){
						meanPosition[j]+=jellyfishVector[i].position[j]/jellyfishVector.size();
					}
				}
				else{
					speciesVector=resetVector;
					for(int stopTime=0;stopTime<numTimepoints;stopTime++){
						pVav_RungeKutta(jellyfishVector[i],speciesVector,stoppingTimes[stopTime],stoppingTimes[stopTime+1],deltaT);
						testData[stopTime]=speciesVector;
					}
					jellyfishVector[i].currentFitnessValue=meanNoNoiseFitness(trueData,testData);
					for(int j=0;j<(int)jellyfishVector[i].position.size();j++){
						meanPosition[j]+=jellyfishVector[i].position[j]/jellyfishVector.size();
					}
				}
			}

			for(int i=0;i<(int)jellyfishVector.size();i++){
				if(jellyfishVector[i].currentFitnessValue<bestFitnessValue){
					minIndex=i;
					bestFitnessValue=jellyfishVector[i].currentFitnessValue;
				}
			}

			for(int i=0;i<(int)jellyfishVector.size();i++){
				int alternateIndex(i);
				do{
					alternateIndex=generator()%jellyfishVector.size();
				}while(alternateIndex==i);
				jellyfishVector[i].motionControl=cControl(curTime,maxIterations,mainUniform());
				if(i!=minIndex){
					jellyfishVector[i].moveJellyfish(jellyfishVector[minIndex],meanPosition,jellyfishVector[alternateIndex],jellyfishRNG);
				}
			}
			
		}

		outFile<<jellyfishVector[minIndex].currentFitnessValue<<" ";
		for(int i=0;i<(int)jellyfishVector[minIndex].position.size();i++){
			outFile<<jellyfishVector[minIndex].position[i]<<" ";
		}
		outFile<<endl;





	}

	outFile.close();



	return 0;
}

int loadPvavInputs(vector<double>& speciesVector, vector<double>& initParameters, vector<double>&stoppingTimes, vector<tuple<double,double> >& bounds, string inFile){
	int errorOut(0);
	ifstream inData(inFile);
	int indexLoop(0);
	double doubleHold(0);
	double doubleHold2(0);
	int intHold(0);
	inData>>indexLoop;
	if(indexLoop!=speciesVector.size()){
		errorOut=1;
		return errorOut;
	}
	for(int i=0;i<indexLoop;i++){
		inData>>intHold;
		speciesVector[i]=intHold;
	}
	inData>>indexLoop;
	if(initParameters.size()!=indexLoop){
		errorOut=1;
		return errorOut;
	}
	for(int i=0;i<indexLoop;i++){
		inData>>doubleHold;
		initParameters[i]=doubleHold;
	}
	inData>>indexLoop;
	stoppingTimes.resize(indexLoop);
	for(int i=0;i<indexLoop;i++){
		inData>>doubleHold;
		stoppingTimes[i]=doubleHold;
	}
	bounds.resize(initParameters.size());
	for(int i=0;i<(int)bounds.size();i++){
		bounds[i]=make_tuple(initParameters[i]/10.,initParameters[i]*10.);
	}
	return errorOut;

}


vector<tuple<double,double> > pVavSetBounds(vector<double>& inSpeciesCounts,vector<double>& inTimes){
	double overallRateLimit(600); //10 per second
	vector<tuple<double,double> > outBounds(6);
	double maxTime(30);
	double k0Min(0), k0Max(0);
	k0Max=overallRateLimit/(inSpeciesCounts[0]*inSpeciesCounts[1]);
	k0Min=(1./maxTime)/(inSpeciesCounts[0]*inSpeciesCounts[1]);
	outBounds[0]=make_tuple(k0Min,k0Max);
	
	double k1Min(0), k1Max(0);
	double limitSykVav((min(inSpeciesCounts[0],inSpeciesCounts[1])));
	k1Max=overallRateLimit/limitSykVav;
	k1Min=(1./maxTime)/(limitSykVav);
	outBounds[1]=make_tuple(k1Min,k1Max);

	outBounds[2]=make_tuple(k1Min,k1Max);

	double k3Min(0), k3Max(0);
	k3Max=overallRateLimit/(inSpeciesCounts[1]*inSpeciesCounts[4]);
	k3Min=(1./maxTime)/(inSpeciesCounts[1]*inSpeciesCounts[4]);
	outBounds[3]=make_tuple(k3Min,k3Max);

	double k4Min(0), k4Max(0);
	double limitShpVav(min(inSpeciesCounts[1],inSpeciesCounts[4]));
	k4Max=overallRateLimit/limitShpVav;
	k4Min=(1./maxTime)/limitShpVav;
	outBounds[4]=make_tuple(k4Min,k4Max);

	outBounds[5]=make_tuple(k4Min,k4Max);

	return outBounds;
}

void loadCovariance(vector<double>& outMeans, vector<vector<double> >& inMatrix, string dataPath){
    int inSize(0);
    ifstream inData(dataPath+".txt");
	if(!inData.good()){
		cout<<"failed to load data for covariance matrix"<<endl;
	}
    inData>>inSize;
    outMeans.resize(inSize);
    for(int i=0;i<inSize;i++){
        inData>>outMeans[i];
    }
    vector<vector<double> > interMatrix(inSize,vector<double> (inSize,0));
    for(int i=0;i<(int)interMatrix.size();i++){
        for(int j=0;j<(int)interMatrix[i].size();j++){
            inData>>interMatrix[i][j];
        }
    }
    inData.close();
    inMatrix=interMatrix;
}
