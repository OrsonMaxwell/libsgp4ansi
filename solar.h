/*
 * solar.h - Solar system bodies position routines for libsgp4ansi.
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

#ifndef SOLAR_H_
#define SOLAR_H_

// Find julian century from a given time and date
double
unix2century
(
  time_t  time,
  float   time_ms
);

// Determine Earth nutation parametres on a given date
void
nutation
(
  double  T,
  double* D,
  double* M,
  double* Mdot,
  double* F,
  double* Omega,
  double* L,
  double* Ldot
);

// Find (coarse) position of the Sun at given Julian time in equatorial frame
vec3
solar_pos
(
  time_t  time,
  float   time_ms,
  double* lambda
);

// Find (coarse) position of the Moon at given Julian time in equatorial frame
vec3
lunar_pos
(
  time_t  time,
  float   time_ms,
  double* lambda
);

#endif /* SOLAR_H_ */
