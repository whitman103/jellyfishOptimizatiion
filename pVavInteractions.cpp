#include <vector>

#include "pVavInteractions.h"
#include "JellyfishClass.h"

using namespace std;

void pVav_RungeKutta(Jellyfish& currentParticle, vector<double>& speciesVec, double currentTime, double stoppingTime, double deltaT){
	int numSpecies=speciesVec.size();
    int n=(int)(((stoppingTime-currentTime))/(deltaT));

    vector<double (*)(Jellyfish&,vector<double>&)> interactionPointer=currentParticle.interactionFunctions;

    for(int t=0;t<n;t++){
        vector<double> interSpecies(numSpecies,0);
        vector<double> k1(numSpecies,0);
        for(int i=0;i<(int)k1.size();i++){
            k1[i]=deltaT*interactionPointer[i](currentParticle,speciesVec);
			interSpecies[i]=speciesVec[i]+k1[i]/2.;
        }
        vector<double> k2(numSpecies,0);
        for(int i=0;i<(int)k2.size();i++){
            k2[i]=deltaT*interactionPointer[i](currentParticle,interSpecies);
			interSpecies[i]=speciesVec[i]+k2[i]/2.;
        }
        vector<double> k3(numSpecies,0);
        for(int i=0;i<(int)k3.size();i++){
            k3[i]=deltaT*interactionPointer[i](currentParticle,interSpecies);
			interSpecies[i]=speciesVec[i]+k3[i];
        }
        vector<double> k4(numSpecies,0);
        for(int i=0;i<(int)k4.size();i++){
            k4[i]=deltaT*interactionPointer[i](currentParticle,interSpecies);
        }
        for(int i=0;i<numSpecies;i++){
            speciesVec[i]+=(k1[i]/6.+k2[i]/3.+k3[i]/3.+k4[i]/6.);
        }
    }
};

double dSyk(Jellyfish& inParticle, vector<double>& currentSpecies){
	double complexForm(-1.*inParticle.position[0]*currentSpecies[0]*currentSpecies[1]);
	double complexBreak(inParticle.position[1]*currentSpecies[2]);
	double phosphor(inParticle.position[2]*currentSpecies[2]);
	return (complexForm+complexBreak+phosphor);
}

double dVav(Jellyfish& inParticle, vector<double>& currentSpecies){
	double complexForm(-1.*inParticle.position[0]*currentSpecies[0]*currentSpecies[1]);
	double complexBreak(inParticle.position[1]*currentSpecies[2]);
	double ShpBreak(inParticle.position[5]*currentSpecies[5]);
	return (complexForm+complexBreak+ShpBreak);
}
	
double dSV(Jellyfish& inParticle, vector<double>& currentSpecies){
	double complexForm(inParticle.position[0]*currentSpecies[0]*currentSpecies[1]);
	double complexBreak(-1.*inParticle.position[1]*currentSpecies[2]);
	double phosphor(-1.*inParticle.position[2]*currentSpecies[2]);
	return (complexForm+complexBreak+phosphor);
}
	
double dpVav(Jellyfish& inParticle, vector<double>& currentSpecies){
	double phosphor(inParticle.position[2]*currentSpecies[2]);
	double complexForm(-1.*inParticle.position[3]*currentSpecies[3]*currentSpecies[4]);
	double complexBreak(inParticle.position[4]*currentSpecies[5]);
	return (phosphor+complexForm+complexBreak);
}
	
double dSHP1(Jellyfish& inParticle, vector<double>& currentSpecies){
	double complexForm(-1.*inParticle.position[3]*currentSpecies[4]*currentSpecies[3]);
	double complexBreak(inParticle.position[4]*currentSpecies[5]);
	double phosphor(inParticle.position[5]*currentSpecies[5]);
	return (phosphor+complexBreak+complexForm);
}
	
double dSHP1Vav(Jellyfish& inParticle, vector<double>& currentSpecies){
	double complexForm(inParticle.position[3]*currentSpecies[3]*currentSpecies[4]);
	double complexBreak(-1.*inParticle.position[4]*currentSpecies[5]);
	double phosphor(-1.*inParticle.position[5]*currentSpecies[5]);
	return (complexForm+complexBreak+phosphor);
}

double meanNoNoiseFitness(vector<vector<double> >& trueData, vector<vector<double> >& testData){
	double outValue(0);
	for(int i=0;i<(int)trueData.size();i++){
		for(int j=0;j<(int)trueData[i].size();j++){
			outValue+=pow(trueData[i][j]-testData[i][j],2);
		}
	}
	return outValue;
}

