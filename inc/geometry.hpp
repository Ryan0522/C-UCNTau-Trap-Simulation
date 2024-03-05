#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

/**
 * @see geometry.cpp
 * Calculates the cross product of two functions.
*/
std::vector<double> cross(std::vector<double> a, std::vector<double> b);

/**
 * @see geometry.cpp
 * Normalizes a vector.
*/
void normalize(std::vector<double> &a);

/**
 * @see geometry.cpp
 * Working on it.
*/
double zOffDipCalc(double t);

/**
 * @see geometry.cpp
 * Simulates the reflection of the neutron off the detector.
*/
void reflect(std::vector<double> &state, std::vector<double> norm, std::vector<double> tang);

/**
 * @see geometry.cpp
 * Simulates the reflection of the neutron off the detector.
*/
bool checkDagHit(double x, double y, double z, double zOff);

/**
 * @see geometry.cpp
 * Checks if the point 3D is hit by the detector.
*/
int checkClean(std::vector<double> state, std::vector<double> prevState, double cleanHeight);

/**
 * @see geometry.cpp
 * Calculates the distance between the bottom edge of the detector and the Halbach array.
*/
double calcDagZeta(double x, double y, double z, double zOff);

/**
 * @see geometry.cpp
 * Checks if a point in 3D space hit a specific part of the house.
*/
bool checkHouseHitLow(double x, double y, double z, double zOff);

/**
 * @see geometry.cpp
 * Checks if a point in 3D space hit a specific part of the house.
*/
bool checkHouseHitHigh(double x, double y, double z, double zOff);

/**
 * @see geometry.cpp
 * Initializes a state close to ref for further analysis in perturbation of initial condition.
 * The specific adjustments made to prepare the pair vector is through the momentum components. 
*/
std::vector<double> initializeLyapState(std::vector<double> ref);

/**
 * @see geometry.cpp
 * Adjusts the elements of pair vector so it gets closer to ref. 
*/
void resetStates(std::vector<double> ref, std::vector<double> &pair);

/**
 * @see geometry.cpp
 * Calculates the distance between two vectors.
*/
double distance(std::vector<double> ref, std::vector<double> pair);


#endif /* GEOMETRY_H */