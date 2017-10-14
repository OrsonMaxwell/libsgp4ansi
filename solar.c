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

#include <math.h>

#include "libsgp4ansi.h"
#include "epoch.h"

/*
 * Find (coarse) position of the Sun at given Julian time in equatorial frame
 *
 * Inputs:  time    - Unix time
 *          time_ms - Millisecond portion of time
 * Returns: azelrng - Equatorial coordinates vector (rad, rad, au)
 */
vec3
solar_pos
(
  time_t time,
  float  time_ms
)
{
  vec3 azelrng;

  // Julian century
  double T  = (unix2jul(time, time_ms) - J2000) / 36525;
  double T2 = T * T;
  double T3 = T * T * T;

  // Geometric mean longitude
  double Lo = 280.46645 + 36000.76983 * T + 0.0003032 * T2;

  // Mean anomaly
  double M = 357.52910 + 35999.0503 * T - 0.0001559 * T2 - 0.00000048 * T3;

  // Earth orbit eccentricity
  double e = 0.016708617 - 0.000042037 * T - 0.0000001236 * T2;

  // Sun center
  double C = (1.9146 - 0.004817 * T - 0.000014 * T2) * sin(M * DEG2RAD)
           + (0.019993 - 0.000101 * T) * sin(2 * M * DEG2RAD)
           + 0.00029 * sin(3 * M * DEG2RAD);

  // True longitude
  double Theta = Lo + C;

  // True anomaly
  double v = M + C;

  // Nutation and abberation correction factor
  double Omega = 125.04 - 1934.136 * T;

  // Apparent longitude
  double lambda = Theta - 0.00569 - 0.00478 * sin(Omega * DEG2RAD);

  // Mean oliquity of the ecliptic
  double epsilono = 23.43929 - 0.01300417 * T - 1.638889e-7 * T2
                  + 5.036111e-7 * T3;

  // True obliquity of the ecliptic
  double epsilon = epsilono + 0.00256 * cos(Omega * DEG2RAD);

  // Apparent right ascension and declination
  double alpha = atan2(cos((epsilon) * DEG2RAD) * sin(lambda * DEG2RAD),
                       cos(lambda * DEG2RAD));
  double delta = asin(sin((epsilon) * DEG2RAD) * sin(lambda * DEG2RAD));

  // Radius vector from Earth to Sun, au
  double R = (1.000001018 * (1 - e * e)) / (1 + e * cos(v * DEG2RAD));

  azelrng.az  = (alpha < 0)?(alpha + TAU):(alpha);
  azelrng.el  = delta;
  azelrng.rng = R;

  return azelrng;
}

/*
 * Find (coarse) position of the Moon at given Julian time in equatorial frame
 *
 * Inputs:  time    - Unix time
 *          time_ms - Millisecond portion of time
 * Returns: azelrng - Equatorial coordinates vector (rad, rad, au)
 */
vec3
solar_pos
(
  time_t time,
  float  time_ms
)
{
  vec3 azelrng;

  // Julian century
  double T  = (unix2jul(time, time_ms) - J2000) / 36525;
  double T2 = T * T;
  double T3 = T * T * T;


  return azelrng;
}
