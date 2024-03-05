#ifndef TRACK_GEN_H
#define TRACK_GEN_H

#include <vector>
extern "C" {
    #include "../inc/fields_nate.h"
}

std::vector<double> randomPointTrapOptimum(trace tr);
std::vector<double> randomPointTrapEdE(trace tr);
std::vector<double> randomPointTrapOptimumCleanable(trace tr);
std::vector<double> randomPointTrapEdECleanable(trace tr);
std::vector<double> randomPointTrapOptimumOnlyCleanable(trace tr);

#endif /* TRACK_GEN_H */