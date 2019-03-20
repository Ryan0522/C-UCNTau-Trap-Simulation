#ifndef TRACKUCN_H
#define TRACKUCN_H

#include <vector>

extern "C" {
    #include "../inc/fields_nate.h"
}

#define NRECORDS 50

typedef struct noabsResult {
    double energy;
    double theta;
    float times[NRECORDS];
    float ePerps[NRECORDS];
    float zetas[NRECORDS];
} noabsResult;

typedef struct fixedResult {
    double energy;
    double theta;
    double t;
    double settlingT;
    double ePerp;
    double x;
    double y;
    double z;
    double zOff;
    int nHit;
    int nHitHouseLow;
    int nHitHouseHigh;
    double eStart;
    double deathTime;
} fixedResult;

typedef struct cleanResult {
    float energy;
    float theta;
    float t;
    float x;
    float y;
    float z;
    int code;
} cleanResult;

fixedResult fixedEffDaggerHitTime(std::vector<double> state, double dt, trace tr);
fixedResult fixedEffDaggerHitTime_reinsert(std::vector<double> state, double dt, trace tr);
fixedResult fixedEffDaggerHitTime_PSE(std::vector<double> state, double dt, trace tr);
cleanResult cleanTime(std::vector<double> state, double dt, trace tr);
noabsResult noabsCleanTime(std::vector<double> state, double dt, trace tr);
noabsResult daggerHitTimes(std::vector<double> state, double dt, trace tr);

void trackAndPrint(std::vector<double> state, double dt, trace tr);

#endif /* TRACKUCN_H */