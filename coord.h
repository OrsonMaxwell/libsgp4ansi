/*
 * transform.h - Coordinate transformation routines for libsgp4ansi.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#ifndef COORD_H_
#define COORD_H_

// Return unity multiplier with the sign of the argument
int
signof(double);

// Transform position and velocity vectors from TEME to ECEF frame of reference
void
teme2ecef(vec3*, vec3*, double, vec3*, vec3*);

// Transform ECEF position to geodetic latitude, longitude, and altitude
void
ecef2latlonalt(vec3*, double, unsigned int, double, vec3*);

// Transform geodetic latitude, longitude, and altitude to ECEF position vector
void
latlonalt2ecef(vec3*, vec3*);

// Calculate range between observer and satellite from their ECEF vectors
double
ecef2range(vec3*, vec3*);

void rv_tradec
(
vec3* rijk, vec3* vijk, vec3* rsijk,
int direct,
double* rho, double* trtasc, double* tdecl,
double* drho, double* dtrtasc, double* dtdecl
);

#endif /* COORD_H_ */
