#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include "../inc/constants.h"
#include "../inc/track_gen.hpp"
#include "../setup.h"

extern "C" {
    #include "../inc/xorshift.h"
    #include "../inc/fields_nate.h"
}

std::vector<double> randomPointTrapEdE(trace tr) {
    std::vector<double> state(6);
    double maxEnergy = GRAV*MASS_N*0.3444;
    double maxP = sqrt(2*MASS_N*maxEnergy);
    
    double t = 0.0;    
    
    double energy;
    while(true) {
        energy = maxEnergy * nextU01();
        if(energy < 0.5/JTONEV) {
            continue;
        }
        if(nextU01() < energy/maxEnergy) {
            break;
        }
    }
    
    double totalU;
    do {
        state[2] = -1.464413669130002;
        state[0] = nextU01()*0.15 - 0.075;
        state[1] = nextU01()*0.15 - 0.075;
        potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
        totalU = totalU - MINU;
    } while(totalU >= energy);
    
    double targetP = sqrt(2.0*MASS_N*(energy - totalU));
    
    double u1 = nextU01();
    double u2 = nextU01();
    double theta = asin(sqrt(u1));
    double phi = 2 * M_PI * u2;
    
    state[3] = sin(theta)*cos(phi);
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
    
    double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    
    state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
    
    return state;
}

std::vector<double> randomPointTrapOptimum(trace tr) {
    std::vector<double> state(6);
    double maxEnergy = GRAV*MASS_N*0.3444;
    double maxP = sqrt(2*MASS_N*maxEnergy);
    
    double t = 0.0;    
    
    double energy;
    while(true) {
        energy = maxEnergy * nextU01();
        if(energy < ECUT/JTONEV) {
            continue;
        }
        if(nextU01() < pow(energy/maxEnergy, EPOW)) {
            break;
        }
    }
    
    double totalU;
    do {
        state[2] = -1.464413669130002;
        state[0] = nextU01()*0.15 - 0.075;
        state[1] = nextU01()*0.15 - 0.075;
        potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
        totalU = totalU - MINU;
    } while(totalU >= energy);
    
    double targetP = sqrt(2.0*MASS_N*(energy - totalU));
    
    double theta;
    while(true) {
        double u1 = nextU01();
        theta = asin(sqrt(u1));
        if(nextU01() < pow(cos(theta), THETAPOW)) {
            break;
        }
    }
    
    double u2 = nextU01();
    double phi = 2 * M_PI * u2;
    
    state[3] = sin(theta)*cos(phi);
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
    
    double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    
    state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
    
    return state;
}

std::vector<double> randomPointTrapOptimumCleanable(trace tr) {
    std::vector<double> state(6);
    double maxEnergy = GRAV*MASS_N*0.45;
    double maxP = sqrt(2*MASS_N*maxEnergy);
    
    double t = 0.0;    
    
    double energy;
    while(true) {
        energy = maxEnergy * nextU01();
        if(energy < ECUT/JTONEV) {
            continue;
        }
        if(nextU01() < pow(energy/maxEnergy, EPOW)) {
            break;
        }
    }
    
    double totalU;
    do {
        state[2] = -1.464413669130002;
        state[0] = nextU01()*0.15 - 0.075;
        state[1] = nextU01()*0.15 - 0.075;
        potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
        totalU = totalU - MINU;
    } while(totalU >= energy);
    
    double targetP = sqrt(2.0*MASS_N*(energy - totalU));
    
    double theta;
    while(true) {
        double u1 = nextU01();
        theta = asin(sqrt(u1));
        if(nextU01() < pow(cos(theta), THETAPOW)) {
            break;
        }
    }
    
    double u2 = nextU01();
    double phi = 2 * M_PI * u2;
    
    state[3] = sin(theta)*cos(phi);
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
    
    double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    
    state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
    
    return state;
}

std::vector<double> randomPointTrapEdECleanable(trace tr) {
    std::vector<double> state(6);
    double maxEnergy = GRAV*MASS_N*0.45;
    double maxP = sqrt(2*MASS_N*maxEnergy);
    
    double t = 0.0;    
    
    double energy;
    while(true) {
        energy = maxEnergy * nextU01();
        if(energy < 0.5/JTONEV) {
            continue;
        }
        if(nextU01() < energy/maxEnergy) {
            break;
        }
    }
    
    double totalU;
    do {
        state[2] = -1.464413669130002;
        state[0] = nextU01()*0.15 - 0.075;
        state[1] = nextU01()*0.15 - 0.075;
        potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
        totalU = totalU - MINU;
    } while(totalU >= energy);
    
    double targetP = sqrt(2.0*MASS_N*(energy - totalU));
    
    double u1 = nextU01();
    double u2 = nextU01();
    double theta = asin(sqrt(u1));
    double phi = 2 * M_PI * u2;
    
    state[3] = sin(theta)*cos(phi);
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
    
    double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    
    state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
    
    return state;
}

std::vector<double> randomPointTrapOptimumOnlyCleanable(trace tr) {
    std::vector<double> state(6);
    double maxEnergy = GRAV*MASS_N*0.45;
    double maxP = sqrt(2*MASS_N*maxEnergy);
    
    double t = 0.0;    
    
    double energy;
    while(true) {
        energy = maxEnergy * nextU01();
        if(energy < ECLEAN) {
            continue;
        }
        if(nextU01() < pow(energy/maxEnergy, EPOW)) {
            break;
        }
    }
    
    double totalU;
    do {
        state[2] = -1.464413669130002;
        state[0] = nextU01()*0.15 - 0.075;
        state[1] = nextU01()*0.15 - 0.075;
        potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
        totalU = totalU - MINU;
    } while(totalU >= energy);
    
    double targetP = sqrt(2.0*MASS_N*(energy - totalU));
    
    double theta;
    while(true) {
        double u1 = nextU01();
        theta = asin(sqrt(u1));
        if(nextU01() < pow(cos(theta), THETAPOW)) {
            break;
        }
    }
    
    double u2 = nextU01();
    double phi = 2 * M_PI * u2;
    
    state[3] = sin(theta)*cos(phi);
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
    
    double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    
    state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
    
    return state;
}
