#ifndef UCNTFIELDS_F_H
#define UCNTFIELDS_F_H

#define SAMPDT 0.0004

typedef struct trace {
    double* x;
    double* y;
    double* z;
    int num;
} trace;

/**
 * @see fields_nate.c
 * Simulate movement and/or motion based on recorded traces.
*/
void shift(double *x, double *y, double *z, double t, trace* tr);

/**
 * @see fields_nate.c
 * Computes force experineced by a point with coordinates (x_in, y_in, z_in).
*/
void force(double *x_in, double *y_in, double *z_in, double *fx, double *fy, double *fz, double *totalU, double* t, trace* tr);

/**
 * @see fields_nate.c
 * Computes the total magnetic field strength at a given point (x_in, y_in, z_in).
*/
void fieldstrength(double *x_in, double *y_in, double *z_in, double *totalB, double* t, trace* tr);

/**
 * @see fields_nate.c
 * Computes the total potential energy at a given point (x_in, y_in, z_in).
*/
void potential(double *x_in, double *y_in, double *z_in, double *totalU, double* t, trace* tr);

#endif /* UCNTFIELDS_F_H */
