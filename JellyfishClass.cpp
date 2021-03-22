#include <vector>
#include <tuple>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "JellyfishClass.h"


Jellyfish::Jellyfish(vector<double> inPosition, vector<tuple<double,double> > inBounds){
	this->position=inPosition;
	this->bounds=inBounds;
	this->motionControl=0;
	this->motionType=0;
	this->beta=3;
	this->gamma=0.1;
	this->currentFitnessValue=0;
};

Jellyfish::Jellyfish(){
	this->position={4,2};
	this->bounds={make_tuple(2,3)};
	this->motionControl=0;
	this->motionType=0;
	this->beta=3;
	this->gamma=0.1;
	this->currentFitnessValue=0;
}

Jellyfish::~Jellyfish(){
};

void Jellyfish::moveJellyfish(Jellyfish& bestJellyfish, vector<double>& meanPosition, Jellyfish& alternateJellyfish,boost::variate_generator<boost::mt19937, boost::uniform_real<double> >& RNG){
	if(motionControl<0.5){
		if(RNG()>(1-motionControl)){
			this->typeAMotion(RNG);
		}
		else{
			this->typeBMotion(alternateJellyfish,RNG);
		}
	}
	else{
		this->currentMotion(bestJellyfish, meanPosition,RNG);
	}
	this->checkBoundaries();
};

void Jellyfish::typeAMotion(boost::variate_generator<boost::mt19937, boost::uniform_real<double> >& RNG){
	double randValue(RNG());
	for(int i=0;i<(int)position.size();i++){
		auto[lower,upper]=this->bounds[i];
		position[i]+=(this->gamma)*(upper-lower)*randValue;
	}
};

void Jellyfish::typeBMotion(Jellyfish& alternateJellyfish,boost::variate_generator<boost::mt19937, boost::uniform_real<double> >& RNG){
	double randValue(RNG());
	for(int i=0;i<(int)position.size();i++){
		if(alternateJellyfish.currentFitnessValue<this->currentFitnessValue){
			position[i]+=randValue*(alternateJellyfish.position[i]-this->position[i]);
		}
		else{
			position[i]+=randValue*(this->position[i]-alternateJellyfish.position[i]);
		}
	}
};

void Jellyfish::currentMotion(Jellyfish& bestJellyfish, vector<double>& meanPosition,boost::variate_generator<boost::mt19937, boost::uniform_real<double> >& RNG){
	vector<double> Trend(meanPosition.size(),0);
	double randValue(RNG());
	for(int i=0;i<(int)Trend.size();i++){
		Trend[i]=this->beta*randValue*meanPosition[i];
	}
	randValue=RNG();
	for(int i=0;i<(int)position.size();i++){
		position[i]+=randValue*(bestJellyfish.position[i]-Trend[i]);
	}
};

void Jellyfish::checkBoundaries(){
	for(int i=0;i<(int)position.size();i++){
		auto[lower,upper]=this->bounds[i];
		if(position[i]>upper){
			position[i]=position[i]-upper+lower;
		}
		if(position[i]<lower){
			position[i]=position[i]-lower+upper;
		}
	}
};


