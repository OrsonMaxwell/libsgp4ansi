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
  double T2 = T  * T;
  double T3 = T2 * T;

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
lunar_pos
(
  time_t time,
  float  time_ms
)
{
  vec3 azelrng;

  // Julian century
  double T  = (unix2jul(time, time_ms) - J2000) / 36525;
  double T2 = T  * T;
  double T3 = T2 * T;
  double T4 = T3 * T;

  // Mean longitude of the Moon
  double Ldot = 218.3164591 + 481267.88134236 * T - 0.0013268 * T2
              + T3 / 538841 - T4 / 65194000;

  // Mean elongation of the Moon
  double D = 297.8502042 + 445267.1115168 * T - 0.00163 * T2 + T3 / 545868
           - T4 / 113065000;

  // Mean anomaly of the Sun
  double M = 357.5291092 + 35999.0502909 * T - 0.0001536 * T2 - T3 / 24490000;

  // Mean nomaly of the Moon
  double Mdot = 134.9634114 + 477198.8676313 * T + 0.008997 * T2 + T2 / 69699
              - T4 / 14712000;

  // Argument of latitude of the Moon
  double F = 93.2720993 + 483202.0175273 * T - 0.0034029 * T2 - T3 / 3526000
           + T4 / 863310000;

  // Aux arguments
  double A1 = 119.75 + 131.849    * T;
  double A2 = 53.09  + 479264.29  * T;
  double A3 = 313.45 + 481266.484 * T;

  // Eccentricity correction coefficient
  double E = 1 - 0.002516 * T - 0.0000074 * T2;

  // Table periodic terms
  double Sigmal = 0;
  double Sigmar = 0;
  double Sigmab = 0;

  Sigmal +=  6288774 * sin(Mdot);
  Sigmal +=  1274027 * sin(2 * D - Mdot);
  Sigmal +=  658314  * sin(2 * D);
  Sigmal +=  213618  * sin(2 * Mdot);
  Sigmal += -185116  * sin(M);
  Sigmal += -114332  * sin(F);
  Sigmal +=  58793   * sin(2 * D - 2 * Mdot);
  Sigmal +=  57066   * sin(2 * D - M - Mdot);
  Sigmal +=  53322   * sin(2 * D + Mdot);
  Sigmal +=  45758   * sin(2 * D - M);
  Sigmal += -40923   * sin(M - Mdot);
  Sigmal += -34720   * sin(D);
  Sigmal += -30383   * sin(M + Mdot);
  Sigmal +=  15327   * sin(2 * D - 2 * F);
  Sigmal += -12528   * sin(Mdot + 2 * F);
  Sigmal +=  10980   * sin(Mdot - 2 * F);
  Sigmal +=  10675   * sin(4 * D - Mdot);
  Sigmal +=  10034   * sin(3 * Mdot);
  Sigmal +=  8548    * sin(4 * D - 2 * Mdot);
  Sigmal += -7888    * sin(2 * D + M - Mdot);
  Sigmal += -6766    * sin(2 * D + M);
  Sigmal += -5163    * sin(D - Mdot);
  Sigmal +=  4987    * sin(D + M);
  Sigmal +=  4036    * sin(2 * D - M + Mdot);
  Sigmal +=  3994    * sin(2 * D + 2 * Mdot);
  Sigmal +=  3861    * sin(4 * D);
  Sigmal +=  3665    * sin(2 * D - 3 * Mdot);
  Sigmal += -2689    * sin(M + 2 * Mdot);
  Sigmal += -2602    * sin(2 * D - Mdot + 2 * F);
  Sigmal +=  2390    * sin(2 * D - M - 2 * Mdot);
  Sigmal += -2348    * sin(D + Mdot);
  Sigmal +=  2236    * sin(2 * D - 2 * M);
  Sigmal += -2120    * sin(M + 2 * Mdot);
  Sigmal += -2069    * sin(2 * M);
  Sigmal +=  2048    * sin(2 * D - 2 * M - Mdot);
  Sigmal += -1773    * sin(2 * D + Mdot - 2 * F);
  Sigmal += -1595    * sin(2 * D + 2 * F);
  Sigmal +=  1215    * sin(4 * D - M - Mdot);
  Sigmal += -1110    * sin(2 * Mdot + 2 * F);
  Sigmal += -892     * sin(3 * D - Mdot);
  Sigmal += -810     * sin(2 * D + M + Mdot);
  Sigmal +=  759     * sin(4 * D - M - 2 * Mdot);
  Sigmal += -713     * sin(2 * M - Mdot);
  Sigmal += -700     * sin(2 * D + 2 * M - Mdot);
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();
//  Sigmal +=  1     * sin();


  // Right ascension, declination and range
  double lambda = Ldot + Sigmal / 1000000;
  double beta   = Sigmab / 1000000;
  double Delta  = 385000.56 + Sigmar / 1000;

  // Mean equatorial parallax
  double pi = asin(RE / Delta);

  // Nutation in longitude
  double dpsi = 0.00461;

  // Apparent longitude
  double alambda = lambda + dpsi;

  // Nutation and abberation correction factor
  double Omega = 125.04 - 1934.136 * T;

  // Mean oliquity of the ecliptic
  double epsilono = 23.43929 - 0.01300417 * T - 1.638889e-7 * T2
                  + 5.036111e-7 * T3;

  // True obliquity of the ecliptic
  double epsilon = epsilono + 0.00256 * cos(Omega * DEG2RAD);

  // Apparent right ascension, declination and radius vector
  azelrng.ra  = atan2(sin(alambda * DEG2RAD) * cos(epsilon * DEG2RAD)
                    - tan(beta * DEG2RAD) * sin(epsilon * DEG2RAD),
                      cos(alambda * DEG2RAD));
  azelrng.dec = sin(beta * DEG2RAD) * cos(epsilon * DEG2RAD)
              + cos(beta * DEG2RAD) * sin(epsilon * DEG2RAD)
              * sin(alambda * DEG2RAD);
  azelrng.rv  = Delta / AU;

  return azelrng;
}
