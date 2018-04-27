#include "../inc/trackUCN.hpp"
#include "../inc/constants.h"
#include "../inc/symplectic.hpp"
#include "../inc/geometry.hpp"
#include "../inc/quant_refl.hpp"
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <assert.h>

#include "../setup.h"

extern "C" {
    #include "../inc/fields_nate.h"
    #include "../inc/xorshift.h"
}

noabsResult daggerHitTimes(std::vector<double> state, double dt, trace tr) {
    std::vector<double> tang = {0.0, 0.0, 1.0};
    std::vector<double> normPlus = {0.0, 1.0, 0.0};
    std::vector<double> normMinus = {0.0, -1.0, 0.0};
    noabsResult res;
    res.theta = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]));

    double t = 0;
    potential(&state[0], &state[1], &state[2], &(res.energy), &t, &tr);
    res.energy = res.energy - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2*MASS_N);
    
    double settlingTime;
    do {
        settlingTime = -70*log(nextU01());
    } while(settlingTime >= 150);
    
    settlingTime = settlingTime + CLEANINGTIME;

    std::vector<double> prevState(6);
    int numSteps = settlingTime/dt;
    double energy;
    for(int i = 0; i < numSteps; i++) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        if((prevState[2] < -1.5+0.38 && state[2] > -1.5+0.38 && state[1] > 0) || (prevState[2] > -1.5+0.38 && state[2] < -1.5+0.38 && state[1] > 0)) { //cleaned
            res.energy = -1.0;
            return res;
        }
        t = t + dt;
    }
    
    int nHit = 0;
    int nHitHouseLow = 0;
    int nHitHouseHigh = 0;
    while(true) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        t = t + dt;
        if(t >= 3000.0) {
            for(int i = nHit; i < NRECORDS; i++) {
                res.times[i] = t - settlingTime;
                res.ePerps[i] = 0.0;
                res.zetas[i] = 0.0;
            }
            return res;
        }
        
        if((prevState[2] < -1.5+0.38 && state[2] > -1.5+0.38 && state[1] > 0) || (prevState[2] > -1.5+0.38 && state[2] < -1.5+0.38 && state[1] > 0)) { //cleaned
            for(int i = nHit; i < NRECORDS; i++) {
                res.times[i] = t - settlingTime;
                res.ePerps[i] = 0.0;
                res.zetas[i] = 0.0;
            }
            return res;
        }

        if((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) {
            double fracTravel = fabs(prevState[1])/(fabs(state[1]) + fabs(prevState[1]));
            double predX = prevState[0] + fracTravel * (state[0] - prevState[0]);
            double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);
            
            double zOff = zOffDipCalc(t - settlingTime);
            
            if(checkDagHit(predX, 0.0, predZ, zOff)) {
                res.times[nHit] = t - settlingTime;
                res.ePerps[nHit] = state[4]*state[4]/(2*MASS_N);
                res.zetas[nHit] = calcDagZeta(predX, 0.0, predZ, zOff);
                nHit += 1;
                if(nHit >= NRECORDS) {
                    break;
                }
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitLow(predX, 0.0, predZ, zOff)) {
                nHitHouseLow += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitHigh(predX, 0.0, predZ, zOff)) {
                nHitHouseHigh += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
        }
    }
    return res;
}

fixedResult fixedEffDaggerHitTime(std::vector<double> state, double dt, trace tr) {
    std::vector<double> tang = {0.0, 0.0, 1.0};
    std::vector<double> normPlus = {0.0, 1.0, 0.0};
    std::vector<double> normMinus = {0.0, -1.0, 0.0};
    fixedResult res;
    res.theta = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]));
    double t = 0;
    potential(&state[0], &state[1], &state[2], &(res.eStart), &t, &tr);
    res.eStart = res.eStart - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2*MASS_N);
    
    double deathTime = -877.7*log(nextU01());
    
    double settlingTime;
    do {
        settlingTime = -70*log(nextU01());
    } while(settlingTime >= 150);
    
    settlingTime = settlingTime + CLEANINGTIME;
    
    res.settlingT = settlingTime;
    
//    if(deathTime < FIRSTDIPTIME) {
//        res.energy = res.eStart;
//        res.t = t;
//        res.ePerp = state[4]*state[4]/(2*MASS_N);
//        res.x = state[0];
//        res.y = state[1];
//        res.z = state[2];
//        res.zOff = -1;
//        res.nHit = 0;
//        res.nHitHouseLow = 0;
//        res.nHitHouseHigh = 0;
//        res.deathTime = deathTime;
//        return res;
//    }

    std::vector<double> prevState(6);
    int numSteps = settlingTime/dt;
    double energy;
    for(int i = 0; i < numSteps; i++) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        if((prevState[2] < -1.5+0.38 && state[2] > -1.5+0.38 && state[1] > 0) || (prevState[2] > -1.5+0.38 && state[2] < -1.5+0.38 && state[1] > 0)) { //cleaned
            res.energy = energy;
            res.t = t - settlingTime;
            res.ePerp = state[5]*state[5]/(2*MASS_N);
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -2;
            res.nHit = 0;
            res.nHitHouseLow = 0;
            res.nHitHouseHigh = 0;
            res.deathTime = deathTime;
            return res;
        }
        t = t + dt;
    }
    
    int nHit = 0;
    int nHitHouseLow = 0;
    int nHitHouseHigh = 0;
    while(true) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        t = t + dt;
        if(t - settlingTime > deathTime) {
            res.energy = energy;
            res.t = t - settlingTime;
            res.ePerp = state[4]*state[4]/(2*MASS_N);
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -1;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
        if((prevState[2] < -1.5+0.38+0.05 && state[2] > -1.5+0.38+0.05 && state[1] > 0) || (prevState[2] > -1.5+0.38+0.05 && state[2] < -1.5+0.38+0.05 && state[1] > 0)) { //cleaned
            res.energy = energy;
            res.t = t - settlingTime;
            res.ePerp = state[5]*state[5]/(2*MASS_N);
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -3;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
        if(isnan(energy)) {
            res.energy = energy;
            res.t = t - settlingTime;
            res.ePerp = 0.0;
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -4;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
        if((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) {
            double fracTravel = fabs(prevState[1])/(fabs(state[1]) + fabs(prevState[1]));
            double predX = prevState[0] + fracTravel * (state[0] - prevState[0]);
            double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);
            
            double zOff = zOffDipCalc(t - settlingTime);
            
            if(checkDagHit(predX, 0.0, predZ, zOff)) {
                nHit += 1;
                if(absorbMultilayer(state[4]*state[4]/(2*MASS_N), BTHICK, predX, 0.0, predZ, zOff)) {
                    res.energy = energy;
                    res.t = t - settlingTime;
                    res.ePerp = state[4]*state[4]/(2*MASS_N);
                    res.x = predX;
                    res.y = 0.0;
                    res.z = predZ;
                    res.zOff = zOff;
                    res.nHit = nHit;
                    res.nHitHouseLow = nHitHouseLow;
                    res.nHitHouseHigh = nHitHouseHigh;
                    res.deathTime = deathTime;
                    return res;
                }
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitLow(predX, 0.0, predZ, zOff)) {
                nHitHouseLow += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitHigh(predX, 0.0, predZ, zOff)) {
                nHitHouseHigh += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
        }
    }
}

cleanResult cleanTime(std::vector<double> state, double dt, trace tr){
    std::vector<double> tang = {0.0, 0.0, 1.0};
    std::vector<double> normPlus = {0.0, 1.0, 0.0};
    std::vector<double> normMinus = {0.0, -1.0, 0.0};
    cleanResult res;
    res.theta = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]));
    double t = 0;
    
//    double deathTime = -877.7*log(nextU01());
    
    double settlingTime;
    do {
        settlingTime = -70*log(nextU01());
    } while(settlingTime >= 150);
    
    settlingTime = settlingTime + CLEANINGTIME;

    std::vector<double> prevState(6);
    int numSteps = settlingTime/dt;
    double energy;
    for(int i = 0; i < numSteps; i++) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        if((prevState[2] < -1.5+CLEANINGHEIGHT && state[2] > -1.5+CLEANINGHEIGHT && state[1] > 0) || (prevState[2] > -1.5+CLEANINGHEIGHT && state[2] < -1.5+CLEANINGHEIGHT && state[1] > 0)) { //cleaned
            res.energy = energy;
            res.t = t;
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.code = -2;
            return res;
        }
        t = t + dt;
    }
    
    while(true) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        t = t + dt;
//        if(t - settlingTime > deathTime) {
//            res.energy = energy;
//            res.t = t-settlingTime;
//            res.x = state[0];
//            res.y = state[1];
//            res.z = state[2];
//            res.code = -1;
//            return res;
//        }
        if((prevState[2] < -1.5+RAISEDCLEANINGHEIGHT && state[2] > -1.5+RAISEDCLEANINGHEIGHT && state[1] > 0) || (prevState[2] > -1.5+RAISEDCLEANINGHEIGHT && state[2] < -1.5+RAISEDCLEANINGHEIGHT && state[1] > 0)) { //cleaned
            res.energy = energy;
            res.t = t-settlingTime;
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.code = -3;
            return res;
        }
        if(isnan(energy)) {
            res.energy = energy;
            res.t = t-settlingTime;
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.code = -4;
            return res;
        }
        if((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) {
            double fracTravel = fabs(prevState[1])/(fabs(state[1]) + fabs(prevState[1]));
            double predX = prevState[0] + fracTravel * (state[0] - prevState[0]);
            double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);
            
            double zOff = zOffDipCalc(t - settlingTime);
            
            if(checkDagHit(predX, 0.0, predZ, zOff)) {
                if(absorbMultilayer(state[4]*state[4]/(2*MASS_N), BTHICK, predX, 0.0, predZ, zOff)) {
                    res.energy = energy;
                    res.t = t-settlingTime;
                    res.x = state[0];
                    res.y = state[1];
                    res.z = state[2];
                    res.code = 0;
                    return res;
                }
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitLow(predX, 0.0, predZ, zOff)) {
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitHigh(predX, 0.0, predZ, zOff)) {
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
        }
    }
}
