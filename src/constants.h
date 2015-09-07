// Written by Peter Kutz.

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <glm/glm.hpp>

#ifndef __MINMAX_DEFINED
#  define max(a,b)    (((a) > (b)) ? (a) : (b))
#  define min(a,b)    (((a) < (b)) ? (a) : (b))
#endif

#define LERP(a,b,t) (1-t)*a + t*b

// Don't modify the values of these here.
// Modify the values of these in Constants.cpp instead.
extern const int theMillisecondsPerFrame;
extern const int theDim[3];
extern const double theCellSize;
extern const double theAirDensity;
extern const double theBuoyancyAlpha;
extern const double theBuoyancyBeta;	
extern const double theBuoyancyAmbientTemperature;
extern const double theVorticityEpsilon;

//fire simulation
extern const double radius;
extern const glm::dvec3 fireCenter;
extern const double S;
extern const double fuelDensity,gasDensity;

extern const double Tmax; 
extern const double Tignition; 
extern const double Trise; 
extern const double CoolT; 
extern const double Tair; 
extern const double AdvectYk; //constant in equation 16 
extern const double Y_threshold; //when Y reach to Y_threshold, T---->Tmax

extern const double a,b,c;

#endif