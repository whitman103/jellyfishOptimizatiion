#ifndef JELLYFISH_H
#define JELLYFISH_H

#include <vector>
#include <tuple>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>


using namespace std;

class Jellyfish{
	public: 
	Jellyfish(vector<double> inPosition, vector<tuple<double,double> > inBounds);
	Jellyfish();
	~Jellyfish();

	vector<double> position;
	vector<tuple<double,double> > bounds;
	int motionType;
	double motionControl;
	double currentFitnessValue;
	double gamma, beta;

	vector<double (*)(Jellyfish&,vector<double>&)> interactionFunctions;

	void moveJellyfish(Jellyfish& bestJellyfish, vector<double>& meanPosition, Jellyfish& alternateJellyfish, boost::variate_generator<boost::mt19937, boost::uniform_real<double> >& RNG);

	void typeAMotion(boost::variate_generator<boost::mt19937, boost::uniform_real<double> >& RNG);
	void typeBMotion(Jellyfish& alternativeJellyfish,boost::variate_generator<boost::mt19937, boost::uniform_real<double> >& RNG);
	void currentMotion(Jellyfish& bestJellyfish, vector<double>& meanPosition,boost::variate_generator<boost::mt19937, boost::uniform_real<double> >& RNG);
	void checkBoundaries();

};




#endif 