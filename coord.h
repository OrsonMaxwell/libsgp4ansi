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

// Solve Kepler's equation for known true anomaly
double
kepler_newton
(
  const double ecc,
  const double nu
);

// Get classical orbital elements from TEME vectors
extern int
teme2coe
(
  vec3*, vec3*, double*, double*, double*, double*, double*, double*, double*,
  double*, double*, double*, double*
);

// Transform position and velocity vectors from TEME to ECEF frame of reference
void
teme2ecef
(
  const vec3*  posteme,
  const vec3*  velteme,
        double julian,
        vec3*  posecef,
        vec3*  velecef
);

// Transform ECEF position to geodetic latitude, longitude, and altitude
void
ecef2latlonalt
(
  const vec3* posecef,
        vec3* latlonalt
);
// Transform geodetic latitude, longitude, and altitude to ECEF position vector
void
latlonalt2ecef
(
  const vec3* latlonalt,
        vec3* posecef
);

// Calculate range between observer and satellite from their ECEF vectors
double
ecef2range
(
  const vec3* obsposecef,
  const vec3* satposecef
);

#endif /* COORD_H_ */
