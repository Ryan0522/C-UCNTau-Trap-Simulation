#ifndef SYMPLECTIC_H
#define SYMPLECTIC_H

#include <vector>

extern "C" {
    #include "../inc/fields_nate.h"
}

void symplecticStep(std::vector<double> &state, double deltaT, double &energy, double t, trace tr);

#endif /* SYMPLECTIC_H */