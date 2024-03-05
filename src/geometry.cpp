#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <assert.h>
#include "../setup.h"
#include "../inc/constants.h"
#include "../inc/geometry.hpp"

extern "C" {
    #include "../inc/xorshift.h"
}

/**
 * Calculates the cross product of two functions.
 * @param a First vector in the cross product.
 * @param b Second vector in the cross product.
 * @return A vector of doubles as the outcome of the cross product.
*/
std::vector<double> cross(std::vector<double> a, std::vector<double> b) {
    std::vector<double> res(3);
    res[0] = a[1]*b[2]-a[2]*b[1];
    res[1] = a[0]*b[2]-a[2]*b[0];
    res[2] = a[0]*b[1]-a[1]*b[0];
    return res;
}

/**
 * Normalizes a vector.
 * @param a Vector to noramlize
*/
void normalize(std::vector<double> &a) {
    assert(a.size()==3);
    double len = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    a[0] = a[0]/len;
    a[1] = a[1]/len;
    a[2] = a[2]/len;
}

/**
 * Working on definition.
 * @param t working...
 * @return working...
*/
double zOffDipCalc(double t) {
    double acc = 6.0;
    double vel = 1.6;
    double spr = 51200/2.0;
    
    double holdT = HOLDTIME;
    
    double dipHeights[NDIPS] = HEIGHTS; //3 dip
    double dipEnds[NDIPS] = ENDTIMES; //3 dip
    
    int i;
    
    if(t > dipEnds[NDIPS-1]) {
        return 0.01;
    }
    
    for(i = 0; i < NDIPS; i++) {
        if(dipEnds[i] > t) {
            break;
        }
    }
    
    if(i == 0) {
        return dipHeights[0];
    }
    
    double target = dipHeights[i];
    double start = dipHeights[i-1];
    
    double distance = fabs(target - start)*1000000;
    
    double time = (distance - vel*(vel/acc)*spr)/(vel*spr);
    
    double sign = dipHeights[i] - dipHeights[i-1] > 0 ? 1 : -1;
    sign = dipHeights[i] == dipHeights[i-1] ? 0 : sign;
    
    if(time > 0) {
//        printf("Flat\n");
        if((t - dipEnds[i-1]) < (vel/acc)) {
            return dipHeights[i-1]
                + sign*(
                    acc*spr*(t - dipEnds[i-1])*(t - dipEnds[i-1])/2.0
                )/1000000;
        }
        else if((t - dipEnds[i-1]) >= (vel/acc) && (t - dipEnds[i-1]) < ((vel/acc) + time)) {
            return dipHeights[i-1]
                + sign*(
                    (vel/acc)*vel*spr/2.0 + (t - dipEnds[i-1] - vel/acc)*vel*spr
                )/1000000;
        }
        else if((t - dipEnds[i-1]) >= ((vel/acc) + time) && (t - dipEnds[i-1]) < (2*(vel/acc) + time)){
            return dipHeights[i-1]
                + sign*(
                    (vel/acc)*vel*spr/2.0
                    + time*vel*spr
                    + vel*spr*(t - dipEnds[i-1] - time - vel/acc)
                    - acc*spr*(t - dipEnds[i-1] - time - vel/acc)*(t - dipEnds[i-1] - time - vel/acc)/2.0
                )/1000000;
        }
        else {
            return dipHeights[i];
        }
    }
    else {
        double halfTime = sqrt(distance/(acc*spr));
//        printf("UnFlat\n");
        if(t - dipEnds[i-1] < halfTime) {
            return dipHeights[i-1]
                + sign*(
                    acc*spr*(t - dipEnds[i-1])*(t - dipEnds[i-1])/2.0
                )/1000000;
        }
        else if(t - dipEnds[i-1] >= halfTime && t - dipEnds[i-1] < 2*halfTime){
            return dipHeights[i-1]
                + (dipHeights[i]-dipHeights[i-1])/2
                + sign*(
                    (acc*spr*halfTime*(t - dipEnds[i-1] - halfTime))
                    - acc*spr*(t - dipEnds[i-1] - halfTime)*(t - dipEnds[i-1] - halfTime)/2.0
                )/1000000;
        }
        else {
            return dipHeights[i];
        }
    }
}

//double zOffDipCalc(double t) {
//    double z;
//    int i;
//    
//    double holdT = HOLDTIME;
//    double speed;
//    
//    double dipHeights[NDIPS] = HEIGHTS; //3 dip
//    double dipEnds[NDIPS] = ENDTIMES; //3 dip
//    
//    if(t > dipEnds[NDIPS-1]) {
//        return 0.01;
//    }
//    
//    for(i = 0; i < NDIPS; i++) {
//        if(dipEnds[i] > t) {
//            break;
//        }
//    }
//    
//    if(i == 0) {
//        return dipHeights[0];
//    }
//    
//    if(dipHeights[i-1] > dipHeights[i]) {
//        speed = -0.49/13.0;
//        z = dipHeights[i-1] + speed*(t-dipEnds[i-1]);
//        return z > dipHeights[i] ? z : dipHeights[i];
//    }
//    else {
//        speed = 0.49/13.0;
//        z = dipHeights[i-1] + speed*(t-dipEnds[i-1]);
//        return z < dipHeights[i] ? z : dipHeights[i];
//    }
//}

/**
 * Simulates the reflection of the neutron off the detector.
 * @param state state of the neutron
 * @param norm normal vector to the surface of detector.
 * @param tang tangent vector to the surface of detector.
*/
void reflect(std::vector<double> &state, std::vector<double> norm, std::vector<double> tang) {
    double pTarget = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    std::vector<double> tangPrime = cross(norm, tang);
    
    double u1 = nextU01();
    double u2 = nextU01();
    
    double theta = asin(sqrt(u1));
    double phi = 2 * M_PI * u2;
    
    double pN = cos(theta);
    double pT = sin(theta)*cos(phi);
    double pTprime = sin(theta)*sin(phi);
    std::vector<double> newPdir(3);
    newPdir[0] = pN*norm[0] + pT*tang[0] + pTprime*tangPrime[0];
    newPdir[1] = pN*norm[1] + pT*tang[1] + pTprime*tangPrime[1];
    newPdir[2] = pN*norm[2] + pT*tang[2] + pTprime*tangPrime[2];
    
    double pLen = sqrt(newPdir[0]*newPdir[0] + newPdir[1]*newPdir[1] + newPdir[2]*newPdir[2]);
    
    state[3] = newPdir[0] * pTarget / pLen;
    state[4] = newPdir[1] * pTarget / pLen;
    state[5] = newPdir[2] * pTarget / pLen;
}

/**
 * Checks if the point 3D is hit by the detector. (?)
 * @param x x-coordinate.
 * @param y y-coordinate.
 * @param z z-coordinate.
 * @param zOff working...
 * @return bool indicating whether the point is hit.
*/
bool checkDagHit(double x, double y, double z, double zOff) {
    double zeta;
    if(x > 0) {
        zeta = 0.5 - sqrt(x*x + pow(fabs(z - zOff) - 1.0, 2));
    }
    else {
        zeta = 1.0 - sqrt(x*x + pow(fabs(z - zOff) - 0.5, 2));
    }
    if((x > -0.3524) && (x < 0.0476) && (zeta > 0.0) && (z < (-1.5 + zOff + 0.2))) {
        return true;
    }
    return false;
}

/**
 * Working on it
 * @param state working...
 * @param prevState working...
 * @param clearnHeight working...
 * @return working...
*/
int checkClean(std::vector<double> state, std::vector<double> prevState, double cleanHeight) {
    if((prevState[2] < -1.5 + cleanHeight && state[2] > -1.5 + cleanHeight) || (prevState[2] > -1.5 + cleanHeight && state[2] < -1.5 + cleanHeight)) {
        if(state[1] > 0) {
            return 1;
        }
        if(state[1] > -(0.218041/2 + 0.335121 + 0.3556) && state[1] < -(0.218041/2 + 0.335121) && state[0] > (0.115529/2 + 0.212841 - 0.6604) && state[0] < (0.115529/2 + 0.212841)) {
            return 2;
        }
    }
    return 0;
}

/**
 * Calculates the distance between the bottom edge of the detector and the Halbach array. (?)
 * @param x x-coordinate.  
 * @param y y-coordinate.
 * @param z z-coordinate.
 * @return zeta of the detector
*/
double calcDagZeta(double x, double y, double z, double zOff) {
    double zeta;
    if(x > 0) {
        zeta = 0.5 - sqrt(x*x + pow(fabs(z - zOff) - 1.0, 2));
    }
    else {
        zeta = 1.0 - sqrt(x*x + pow(fabs(z - zOff) - 0.5, 2));
    }
    return zeta;
}

/**
 * Checks if a point in 3D space hit a specific part of the house. (?)
 * @param x x-coordinate.
 * @param y y-coordinate.
 * @param z z-coordinate.
 * @param zOff working...
 * @return bool indicating whether the point hits.
*/
bool checkHouseHitLow(double x, double y, double z, double zOff) {
    if(z >= (-1.5 + zOff + 0.2) && z < (-1.5 + zOff + 0.2 + 0.14478) && fabs(x + 0.1524) < (0.40 + 2.0179*(z + 1.5 - zOff - 0.2))/2.0) {
        return true;
    }
    return false;
}

/**
 * Checks if a point in 3D space hit a specific part of the house. (?)
 * @param x x-coordinate.
 * @param y y-coordinate.
 * @param z z-coordinate.
 * @param zOff working...
 * @return bool indicating whether the point hits.
*/
bool checkHouseHitHigh(double x, double y, double z, double zOff) {
    if(z >= (-1.5 + zOff + 0.2 + 0.14478) && z < (-1.5 + zOff + 0.2 + 0.2667) && fabs(x + 0.1524) < 0.69215/2.0) {
        return true;
    }
    return false;
}

/**
 * Initializes a state close to ref for further analysis in perturbation of initial condition. 
 * The specific adjustments made to prepare the pair vector is through the momentum components.
 * @param ref reference vector.
 * @return the initialized state with perturbation.
*/
std::vector<double> initializeLyapState(std::vector<double> ref) {
    double lenP = sqrt(ref[3]*ref[3] + ref[4]*ref[4] + ref[5]*ref[5]);
    
    std::vector<double> pair(6);
    std::vector<double> randomP(3);
    std::vector<double> paraP(3);
    
    randomP[0] = nextU01();
    randomP[1] = nextU01();
    randomP[2] = nextU01();
    
    paraP[0] = ref[3];
    paraP[1] = ref[4];
    paraP[2] = ref[5];
    
    normalize(randomP);
    normalize(paraP);
    
    std::vector<double> perpP = cross(paraP, randomP);
    
    normalize(perpP);
        
    double alpha = EPSILON*EPSILON*PSCALE*PSCALE/(lenP*2);
    double beta = sqrt(EPSILON*EPSILON*PSCALE*PSCALE - (EPSILON*EPSILON*PSCALE*PSCALE/(2*lenP))*(EPSILON*EPSILON*PSCALE*PSCALE/(2*lenP)));
    
    pair[0] = ref[0];
    pair[1] = ref[1];
    pair[2] = ref[2];
    
    pair[3] = ref[3] - alpha*paraP[0] + beta*perpP[0];
    pair[4] = ref[4] - alpha*paraP[1] + beta*perpP[1];
    pair[5] = ref[5] - alpha*paraP[2] + beta*perpP[2];
    
    return pair;
}

/**
 * Adjusts the elements of pair vector so it gets closer to ref. (?)
 * @param ref reference vector.
 * @param pair second vector.
*/
void resetStates(std::vector<double> ref, std::vector<double> &pair) {
    std::vector<double> dist(6);
    dist[0] = pair[0]-ref[0];
    dist[1] = pair[1]-ref[1];
    dist[2] = pair[2]-ref[2];
    dist[3] = pair[3]-ref[3];
    dist[4] = pair[4]-ref[4];
    dist[5] = pair[5]-ref[5];
    
    double scaling = EPSILON/distance(ref, pair);
    
    pair[0] = ref[0] + scaling*dist[0];
    pair[1] = ref[1] + scaling*dist[1];
    pair[2] = ref[2] + scaling*dist[2];
    pair[3] = ref[3] + scaling*dist[3];
    pair[4] = ref[4] + scaling*dist[4];
    pair[5] = ref[5] + scaling*dist[5];
}

/**
 * Calculates the distance between two vectors.
 * 
 * @param ref reference vector.
 * @param pair second vector.
 * @return A double as the distance between ref and pair.
*/
double distance(std::vector<double> ref, std::vector<double> pair) {
    double dist = 0;
    dist += (pair[0]-ref[0])*(pair[0]-ref[0])/(XSCALE*XSCALE);
    dist += (pair[1]-ref[1])*(pair[1]-ref[1])/(YSCALE*YSCALE);
    dist += (pair[2]-ref[2])*(pair[2]-ref[2])/(ZSCALE*ZSCALE);
    dist += (pair[3]-ref[3])*(pair[3]-ref[3])/(PSCALE*PSCALE);
    dist += (pair[4]-ref[4])*(pair[4]-ref[4])/(PSCALE*PSCALE);
    dist += (pair[5]-ref[5])*(pair[5]-ref[5])/(PSCALE*PSCALE);
    
    dist = sqrt(dist);
    return dist;
}
