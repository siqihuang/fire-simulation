// Written by Peter Kutz.

#include "constants.h"

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {4, 4, 1};
#else
//const int theDim[3] = {12, 12, 1};
const int theDim[3] = {30, 30, 3};
//const int theDim[3] = {20, 20, 10};
//const int theDim[3] = {50, 50, 5};
#endif

const double theCellSize = 0.5;

const double theAirDensity = 1.0;

//const double theBuoyancyAlpha = 0.08; // Gravity's effect on the smoke particles.
const double theBuoyancyAlpha = 1.98; // Gravity's effect on the smoke particles.
const double theBuoyancyBeta = 0.05; // Buoyancy's effect due to temperature difference.	
//const double theBuoyancyAmbientTemperature = 0.0; // Ambient temperature.
const double theBuoyancyAmbientTemperature = 3.0; // Ambient temperature.
const double theVorticityEpsilon = 0.10;

//fire simulation

const double radius=5*theCellSize+0.1;
const glm::dvec3 fireCenter=glm::dvec3(theDim[0]/2, 2, 1); 
const double S=0.2;
const double fuelDensity=1;
const double gasDensity=0.2;


const double Tmax = 35.0;   
const double Tignition = 20.0; 
const double Trise = Tmax - Tignition; 
const double CoolT = 2.0;  
const double Tair= 0.0;  

const double AdvectYk = 1.0; //constant in equation 16 
const double Y_threshold = 0.9; 

const double a=3*theCellSize,b=6*theCellSize,c=3*theCellSize;
