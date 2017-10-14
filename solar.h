/*
 * solar.h - Solar system bodies position routines for libsgp4ansi.
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

#ifndef SOLAR_H_
#define SOLAR_H_

// Find (coarse) position of the Sun at given Julian time in equatorial frame
vec3
solar_pos
(
  time_t time,
  float  time_ms
);

// Find (coarse) position of the Moon at given Julian time in equatorial frame
vec3
lunar_pos
(
  time_t time,
  float  time_ms
);

#endif /* SOLAR_H_ */
