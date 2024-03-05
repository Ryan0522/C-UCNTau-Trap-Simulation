#ifndef SYMPLECTIC_H
#define SYMPLECTIC_H

#include <vector>

extern "C" {
    #include "../inc/fields_nate.h"
}

/**
 * @see symplectic.cpp
 * Advances the system's state by one inetegration step using symplectic integration scheme.
*/
void symplecticStep(std::vector<double> &state, double deltaT, double &energy, double t, trace tr);

#endif /* SYMPLECTIC_H */