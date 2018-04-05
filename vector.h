/*
 * vector.h - Vector math routines for libsgp4ansi.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 * Fundamentals of Astrodynamics and Applications, D. Vallado, Second Edition
 * Astronomical Algorithms, Jean Meeus
 * 1980 IAU Theory of nutation
 *
 * Copyright (c) 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#ifndef VECTOR_H_
#define VECTOR_H_

// Find magnitude of a 3D vector
double
vec3_mag
(
  const vec3* arg
);

// Copy a 3D vector
vec3
vec3_copy
(
  const vec3* from
);

// Add two 3D vectors with coefficients
vec3
vec3_add
(
        double c1,
  const vec3*  vect1,
        double c2,
  const vec3*  vect2
);

// Dot product of two 3D vectors
double
vec3_dot
(
  const vec3* vect1,
  const vec3* vect2
);

// Crossing of two 3D vectors
vec3
vec3_cross
(
  const vec3* vect1,
  const vec3* vect2
);

// Planar angle between two 3D vectors
double
vec3_angle
(
  const vec3* vect1,
  const vec3* vect2
);

#endif /* VECTOR_H_ */
