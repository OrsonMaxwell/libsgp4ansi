/*
 * coord.h - Coordinate transformation routines for libsgp4ansi.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 * Fundamentals of Astrodynamics and Applications, D. Vallado, Second Edition
 * Astronomical Algorithms, Jean Meeus
 *
 * Copyright � 2017 Orson J. Maxwell. Please see LICENSE for details.
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
vec3
ecef2geo
(
  const vec3* posecef
);

// Transform geodetic latitude, longitude, and altitude to ECEF position vector
vec3
geo2ecef
(
  const vec3* latlonalt
);

// Find azimuth, elevation and range from ECEF vectors
vec3
ecef2azelrng
(
  const vec3* op,
  const vec3* dp
);

// Convert equatorial vector to az-el-rng vector from ground st-n at given time
vec3
eq2azelrng
(
  const vec3*  radecrv,
  const vec3*  obs_geo,
        time_t timestamp,
        float  time_ms
);

// Transform position vector from equatorial to TEME frame
vec3
eq2teme
(
  const vec3*  radecrv
);

#endif /* COORD_H_ */
