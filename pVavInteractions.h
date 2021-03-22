#ifndef PVAVINT 
#define PVAVINT

#include <vector>
#include "JellyfishClass.h"
using namespace std;

void pVav_RungeKutta(Jellyfish& currentParticle, vector<double>& speciesVec, double currentTime, double stoppingTime, double deltaT);

double dSyk(Jellyfish& inParticle, vector<double> &currentSpecies);

double dVav(Jellyfish& inParticle, vector<double> &currentSpecies);

double dSV(Jellyfish& inParticle, vector<double> &currentSpecies);

double dpVav(Jellyfish& inParticle, vector<double> &currentSpecies);

double dSHP1(Jellyfish& inParticle, vector<double> &currentSpecies);

double dSHP1Vav(Jellyfish& inParticle, vector<double> &currentSpecies);

double meanNoNoiseFitness(vector<vector<double> >& trueData, vector<vector<double> >& testData);


#endif