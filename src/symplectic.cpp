#include "../inc/symplectic.hpp"
#include "../inc/constants.h"
#include <vector>

extern "C" {
    #include "../inc/fields_nate.h"
}

void symplecticStep(std::vector<double> &state, double deltaT, double &energy) {
    double fx, fy, fz, totalU, t, freq;
    t = 0.0;
    freq = 0.0;
    int n;
    n = 0;
    
    const double a[4] = {
        .5153528374311229364e0,
        -.085782019412973646e0,
        .4415830236164665242e0,
        .1288461583653841854e0,
    };
    const double b[4] = {
        .1344961992774310892e0,
        -.2248198030794208058e0,
        .7563200005156682911e0,
        .3340036032863214255e0,
    };
    
    force(&state[0], &state[1], &state[2], &fx, &fy, &fz, &totalU, &t, &freq);
    energy = totalU - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2.0*MASS_N);

    state[3] = state[3] + b[n]*fx*deltaT;
    state[4] = state[4] + b[n]*fy*deltaT;
    state[5] = state[5] + b[n]*fz*deltaT;
    state[0] = state[0] + a[n]*state[3]*deltaT/MASS_N;
    state[1] = state[1] + a[n]*state[4]*deltaT/MASS_N;
    state[2] = state[2] + a[n]*state[5]*deltaT/MASS_N;
    
    for(n = 1; n < 4; n++) {
        force(&state[0], &state[1], &state[2], &fx, &fy, &fz, &totalU, &t, &freq);
        state[3] = state[3] + b[n]*fx*deltaT;
        state[4] = state[4] + b[n]*fy*deltaT;
        state[5] = state[5] + b[n]*fz*deltaT;
        state[0] = state[0] + a[n]*state[3]*deltaT/MASS_N;
        state[1] = state[1] + a[n]*state[4]*deltaT/MASS_N;
        state[2] = state[2] + a[n]*state[5]*deltaT/MASS_N;
    }
}