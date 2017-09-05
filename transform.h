/*
 * transform.h - Coordinate and time transformation routines for libsgp4ansi.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#ifndef TRANSFORM_H_
#define TRANSFORM_H_

// Return unity multiplier with the sign of the argument
int
signof(double);

// Find magnitude of a 3D vector
double
mag(vect*);

// Copy a 3D vector
void
copyvect(vect*, vect*);

// Add two 3D vectors with coefficients
void
addvect(double, vect*, double, vect*, vect*);

// Dot product of two 3D vectors
double
dot(vect*, vect*);

// Convert year and fractional day to unix time
time_t
fractday2unix(unsigned int, double);

// Convert unix time to Julian date
double
unix2jul(time_t*, unsigned int);

// Convert Julian date to Greenwich Sidereal Time
double
jul2gst(double);

// Transform position and velocity vectors from TEME to ECEF frame of reference
void
teme2ecef(vect*, vect*, double, vect*, vect*);

// Transform ECEF position to geodetic latitude, longitude, and altitude
void
ecef2latlonalt(vect*, double, unsigned int, double, vect*);

// Transform geodetic latitude, longitude, and altitude to ECEF position vector
void
latlonalt2ecef(vect*, vect*);

// Calculate range between observer and satellite from their ECEF vectors
double
ecef2range(vect*, vect*);

void rv_tradec
(
vect* rijk, vect* vijk, vect* rsijk,
int direct,
double* rho, double* trtasc, double* tdecl,
double* drho, double* dtrtasc, double* dtdecl
);

#endif /* TRANSFORM_H_ */
