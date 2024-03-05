#ifndef LYAP_H
#define LYAP_H

#include <vector>

typedef struct lyapResult {
    float eStart;
    float eEnd;
    float theta;
    float lce;
} lyapResult;

extern "C" {
    #include "../inc/fields_nate.h"
}

/**
 * @see lyap.cpp
 * Calculates the Lyaponov Exponent of the neutron motion given parameters.
*/
lyapResult calcLyap(std::vector<double> ref, double dt, trace tr, double tStart = 0);

#endif /* LYAP_H */