/*
 * solar.c - Solar system bodies position routines for libsgp4ansi.
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

#include <math.h>
#include <stdio.h>

#include "libsgp4ansi.h"
#include "epoch.h"

/*
 * Find julian century from a given time and date
 *
 * Inputs:  time    - Unix time
 *          time_ms - Millisecond portion of time
 * Outputs: None
 * Returns: Julian century
 */
__attribute__((always_inline))
inline double
unix2century
(
  time_t  time,
  float   time_ms
)
{
  return (unix2jul(time, time_ms) - J2000) / 36525;
}

/*
 * Determine Earth nutation parametres on a given date
 *
 * Inputs:  T        - Julian century
 * Outputs: D        - Mean elongation of the Moon from the Sun, rad
 *          M        - Mean anomaly of the Sun (Earth), rad
 *          Mdot     - Mean anomaly of the Moon, ead
 *          F        - Moon's argument of latitude, rad
 *          Omega    - Longitude of asc.node of Moon's mean orbit, rad
 *          Ldot     - Mean longitude of the Moon, rad
 *          dpsi     - Nutation in longitude, rad
 *          depsilon - Nutation in obliquity, rad
 * Returns: None
 */
void
nutation
(
  double  T,
  double* D,
  double* M,
  double* Mdot,
  double* F,
  double* Omega,
  double* Ldot,
  double* dpsi,
  double* depsilon
)
{
  // For convenience
  double T2 = pow(T, 2);
  double T3 = pow(T, 3);
  double T4 = pow(T, 4);

  // Mean longitude of the Moon
  *Ldot     = 218.3164591 + 481267.88134236 * T - 0.0013268 * T2
            + T3 / 538841 - T4 / 65194000;
  *Ldot     = (*Ldot < 0) ? (fmod(*Ldot, 360) + 360) : (fmod(*Ldot, 360));
  printf("L':         %.6f\n", *Ldot);
  *Ldot    *= DEG2RAD;

  // Mean elongation of the Moon
  *D        = 297.8502042 + 445267.1115168 * T - 0.00163 * T2 + T3 / 545868
            - T4 / 113065000;
  *D        = (*D < 0) ? (fmod(*D, 360) + 360) : (fmod(*D, 360));
  printf("D:          %.6f\n", *D);
  *D       *= DEG2RAD;

  // Mean anomaly of the Sun
  *M        = 357.5291092 + 35999.0502909 * T - 0.0001536 * T2 - T3 / 24490000;
  *M        = (*M < 0) ? (fmod(*M, 360) + 360) : (fmod(*M, 360));
  printf("M:          %.6f\n", *M);
  *M       *= DEG2RAD;

  // Mean nomaly of the Moon
  *Mdot      = 134.9634114 + 477198.8676313 * T + 0.008997 * T2 + T2 / 69699
              - T4 / 14712000;
  *Mdot     = (*Mdot < 0) ? (fmod(*Mdot, 360) + 360) : (fmod(*Mdot, 360));
  printf("M':         %.6f\n", *Mdot);
  *Mdot     *= DEG2RAD;

  // Argument of latitude of the Moon
  *F         = 93.2720993 + 483202.0175273 * T - 0.0034029 * T2 - T3 / 3526000
             + T4 / 863310000;
  *F        = (*F < 0) ? (fmod(*F, 360) + 360) : (fmod(*F, 360));
  printf("F:          %.6f\n", *F);
  *F *= DEG2RAD;

  // Latitude of the ascending node of lunar mean orbit on the ecliptic
  // measured from mean equinox
  *Omega     = 125.04452 - 1934.136261 * T + 0.0020708 * T2 + T3 / 450000;
  *Omega    = (*Omega < 0) ? (fmod(*Omega, 360) + 360) : (fmod(*Omega, 360));
  printf("Om:         %.2f\n", *Omega);

  // Most significant terms
  *dpsi      = 0;
  *dpsi     += (-171996 - 174.2 * T) * sin(*Omega);
  *dpsi     += ( -13187 -   1.6 * T) * sin(-2 * *D + 2 * *F + 2 * *Omega);
  *dpsi     += (  -2274 -   0.2 * T) * sin(2 * *F + 2 * *Omega);
  *dpsi     += (   2062 +   0.2 * T) * sin(2 * *Omega);
  *dpsi     += (   1426 -   3.4 * T) * sin(*M);
  *dpsi     += (    712 +   0.1 * T) * sin(*Mdot);
  *dpsi     += (   -517 +   1.2 * T) * sin(-2 * *D + *M + 2 * *F + 2 * *Omega);
  *dpsi     += (   -386 -   0.4 * T) * sin(2 * *F + *Omega);
  *dpsi     += (   -301            ) * sin(*Mdot + 2 * *F + 2 * *Omega);
  *dpsi     += (   +217 -   0.5 * T) * sin(-2 * *D -*M + 2 * *F + 2 * *Omega);

  *depsilon = 0;
  *depsilon += (92025 + 8.9 * T) * cos(*Omega);
  *depsilon += ( 5736 - 3.1 * T) * cos(-2 * *D + 2 * *F + 2 * *Omega);
  *depsilon += (  977 - 0.5 * T) * cos(2 * *F + 2 * *Omega);
  *depsilon += ( -895 + 0.5 * T) * cos(2 * *Omega);
  *depsilon += (  +54 - 0.1 * T) * cos(*M);
  *depsilon += (   -7          ) * cos(*Mdot);
  *depsilon += (  224 - 0.6 * T) * cos(-2 * *D + *M + 2 * *F + 2 * *Omega);
  *depsilon += (  200          ) * cos(2 * *F + *Omega);
  *depsilon += (  129 - 0.1 * T) * cos(*Mdot + 2 * *F + 2 * *Omega);
  *depsilon += (  -95 + 0.3 * T) * cos(-2 * *D -*M + 2 * *F + 2 * *Omega);
  printf("dp:         %.6f\n", fmod(*dpsi / 36000000, 360));
  printf("de:         %.6f\n", fmod(*depsilon / 36000000, 360));
  // Convert from 0.0001" to radians
  *dpsi       = (*dpsi < 0) ? (*dpsi / 36000000 + 360) : (*dpsi / 36000000);
  *depsilon   = (*depsilon < 0) ? (*depsilon / 36000000 + 360)
                                : (*depsilon / 36000000);
  *dpsi      *= DEG2RAD;
  *depsilon  *= DEG2RAD;
}

/*
 * Find (coarse) position of the Sun at given Julian time in geocentric
 * equatorial frame
 *
 * Inputs:  time    - Unix time
 *          time_ms - Millisecond portion of time
 * Outputs: lambda  - Apparent geocentric longitude of the Sun, rad
 * Returns: redecrv - Equatorial coordinates vector (rad, rad, km)
 */
vec3
solar_pos
(
  time_t  time,
  float   time_ms,
  double* lambda
)
{
  vec3 radecrv;

  // Julian century
  double T  = unix2century(time, time_ms);
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
  double alambda  = Theta - 0.00569 - 0.00478 * sin(Omega * DEG2RAD);
  alambda *= DEG2RAD;

  if (lambda != NULL)
  {
    *lambda = alambda;
  }

  // Mean oliquity of the ecliptic
  double epsilono = 23.43929 - 0.01300417 * T - 1.638889e-7 * T2
                  + 5.036111e-7 * T3;

  // True obliquity of the ecliptic
  double epsilon = epsilono + 0.00256 * cos(Omega * DEG2RAD);

  // Apparent right ascension and declination
  double alpha = atan2(cos((epsilon) * DEG2RAD) * sin(alambda),
                       cos(alambda));
  double delta = asin(sin((epsilon) * DEG2RAD) * sin(alambda));

  // Radius vector from Earth to Sun, au
  double R = (1.000001018 * (1 - e * e)) / (1 + e * cos(v * DEG2RAD));

  radecrv.az  = (alpha < 0)?(alpha + TAU):(alpha);
  radecrv.el  = delta;
  radecrv.rng = R * AU;

  return radecrv;
}

/*
 * Find (coarse) position of the Moon at given Julian time in geocentric
 * equatorial frame
 *
 * Inputs:  time    - Unix time
 *          time_ms - Millisecond portion of time
 * Outputs: lambda  - Apparent geocentric longitude of the Moon, rad
 * Returns: radecrv - Equatorial coordinates vector (rad, rad, km)
 */
vec3
lunar_pos
(
  time_t  time,
  float   time_ms,
  double* lambda
)
{
  vec3 redecrv;

  printf("JD:        %lf\n", unix2jul(time, time_ms));
  // Julian century
  double T  = unix2century(time, time_ms);;
  double T2 = pow(T, 2);
  double T3 = pow(T, 3);
  double T4 = pow(T, 4);

  printf("T:          %.16f\n", T);

  double D, M, Mdot, F, Omega, Ldot, dpsi, depsilon;
  nutation(T, &D, &M, &Mdot, &F, &Omega, &Ldot, &dpsi, &depsilon);

  // Aux arguments
  double A1 = 119.75 + 131.849    * T;
  double A2 = 53.09  + 479264.29  * T;
  double A3 = 313.45 + 481266.484 * T;
  printf("A1:         %.2f\n", fmod(A1, 360));
  printf("A2:         %.2f\n", fmod(A2, 360) + 360);
  printf("A3:         %.2f\n", fmod(A3, 360) + 360);
  A1       *= DEG2RAD;
  A2       *= DEG2RAD;
  A3       *= DEG2RAD;

  // Eccentricity correction coefficient
  double E  = 1 - 0.002516 * T - 0.0000074 * T2;
  double E2 = pow(E, 2);
  printf("E:          %lf\n", E);

  // Table periodic terms
  double Sigmal = 0;
  double Sigmar = 0;
  double Sigmab = 0;

  Sigmal +=  6288774  * sin(Mdot);
  Sigmal +=  1274027  * sin(2 * D - Mdot);
  Sigmal +=  658314   * sin(2 * D);
  Sigmal +=  213618   * sin(2 * Mdot);
  Sigmal += -185116   * sin(M) * E;
  Sigmal += -114332   * sin(2 * F);
  Sigmal +=  58793    * sin(2 * D - 2 * Mdot);
  Sigmal +=  57066    * sin(2 * D - M - Mdot) * E;
  Sigmal +=  53322    * sin(2 * D + Mdot);
  Sigmal +=  45758    * sin(2 * D - M) * E;
  Sigmal += -40923    * sin(M - Mdot) * E;
  Sigmal += -34720    * sin(D);
  Sigmal += -30383    * sin(M + Mdot) * E;
  Sigmal +=  15327    * sin(2 * D - 2 * F);
  Sigmal += -12528    * sin(Mdot + 2 * F);
  Sigmal +=  10980    * sin(Mdot - 2 * F);
  Sigmal +=  10675    * sin(4 * D - Mdot);
  Sigmal +=  10034    * sin(3 * Mdot);
  Sigmal +=  8548     * sin(4 * D - 2 * Mdot);
  Sigmal += -7888     * sin(2 * D + M - Mdot) * E;
  Sigmal += -6766     * sin(2 * D + M) * E;
  Sigmal += -5163     * sin(D - Mdot);
  Sigmal +=  4987     * sin(D + M) * E;
  Sigmal +=  4036     * sin(2 * D - M + Mdot) * E;
  Sigmal +=  3994     * sin(2 * D + 2 * Mdot);
  Sigmal +=  3861     * sin(4 * D);
  Sigmal +=  3665     * sin(2 * D - 3 * Mdot);
  Sigmal += -2689     * sin(M - 2 * Mdot) * E;
  Sigmal += -2602     * sin(2 * D - Mdot + 2 * F);
  Sigmal +=  2390     * sin(2 * D - M - 2 * Mdot) * E;
  Sigmal += -2348     * sin(D + Mdot);
  Sigmal +=  2236     * sin(2 * D - 2 * M) * E2;
  Sigmal += -2120     * sin(M + 2 * Mdot) * E;
  Sigmal += -2069     * sin(2 * M) * E2;
  Sigmal +=  2048     * sin(2 * D - 2 * M - Mdot) * E2;
  Sigmal += -1773     * sin(2 * D + Mdot - 2 * F);
  Sigmal += -1595     * sin(2 * D + 2 * F);
  Sigmal +=  1215     * sin(4 * D - M - Mdot) * E;
  Sigmal += -1110     * sin(2 * Mdot + 2 * F);
  Sigmal += -892      * sin(3 * D - Mdot);
  Sigmal += -810      * sin(2 * D + M + Mdot) * E;
  Sigmal +=  759      * sin(4 * D - M - 2 * Mdot) * E;
  Sigmal += -713      * sin(2 * M - Mdot) * E2;
  Sigmal += -700      * sin(2 * D + 2 * M - Mdot) * E2;
  Sigmal +=  691      * sin(2 * D + M - 2 * Mdot) * E;
  Sigmal +=  596      * sin(2 * D - M - 2 * F) * E;
  Sigmal +=  549      * sin(4 * D + Mdot);
  Sigmal +=  537      * sin(4 * Mdot);
  Sigmal +=  520      * sin(4 * D - M) * E;
  Sigmal += -487      * sin(D - 2 * Mdot);
  Sigmal += -399      * sin(2 * D + M - 2 * F) * E;
  Sigmal += -381      * sin(2 * Mdot - 2 * F);
  Sigmal +=  351      * sin(D + M + Mdot) * E;
  Sigmal += -340      * sin(3 * D - 2 * Mdot);
  Sigmal +=  330      * sin(4 * D - 3 * Mdot);
  Sigmal +=  327      * sin(2 * D - M + 2 * Mdot * E);
  Sigmal += -323      * sin(2 * M + Mdot) * E2;
  Sigmal +=  299      * sin(D + M - Mdot) * E;
  Sigmal +=  294      * sin(2 * D + 3 * Mdot);

  Sigmal += 3958 * sin(A1);
  Sigmal += 1962 * sin(Ldot - F);
  Sigmal += 318  * sin(A2);

  Sigmar += -20905355 * cos(Mdot);
  Sigmar += -3699111  * cos(2 * D - Mdot);
  Sigmar += -2955968  * cos(2 * D);
  Sigmar += -569925   * cos(2 * Mdot);
  Sigmar +=  48888    * cos(M) * E;
  Sigmar += -3149     * cos(2 * F);
  Sigmar +=  246158   * cos(2 * D - 2 * Mdot);
  Sigmar += -152138   * cos(2 * D - M - Mdot) * E;
  Sigmar += -170733   * cos(2 * D + Mdot);
  Sigmar += -204586   * cos(2 * D - M) * E;
  Sigmar += -129620   * cos(M - Mdot) * E;
  Sigmar +=  108743   * cos(D);
  Sigmar +=  104755   * cos(M + Mdot) * E;
  Sigmar +=  10321    * cos(2 * D - 2 * F);
  Sigmar +=  79661    * cos(Mdot - 2 * F);
  Sigmar += -34782    * cos(4 * D - Mdot);
  Sigmar += -23210    * cos(3 * Mdot);
  Sigmar += -21636    * cos(4 * D - 2 * Mdot);
  Sigmar +=  24208    * cos(2 * D + M - Mdot) * E;
  Sigmar +=  30824    * cos(2 * D + M) * E;
  Sigmar += -8379     * cos(D - Mdot);
  Sigmar += -16675    * cos(D + M) * E;
  Sigmar += -12831    * cos(2 * D - M + Mdot) * E;
  Sigmar += -10445    * cos(2 * D + 2 * Mdot);
  Sigmar += -11650    * cos(4 * D);
  Sigmar +=  14403    * cos(2 * D - 3 * Mdot);
  Sigmar += -7003     * cos(M - 2 * Mdot) * E;
  Sigmar +=  10056    * cos(2 * D - M - 2 * Mdot) * E;
  Sigmar +=  6322     * cos(D + Mdot);
  Sigmar += -9884     * cos(2 * D - 2 * M) * E2;
  Sigmar +=  5751     * cos(M + 2 * Mdot) * E;
  Sigmar += -4950     * cos(2 * D - 2 * M - Mdot) * E2;
  Sigmar +=  4130     * cos(2 * D + Mdot - 2 * F);
  Sigmar += -3958     * cos(4 * D - M - Mdot) * E;
  Sigmar +=  3258     * cos(3 * D - Mdot);
  Sigmar +=  2616     * cos(2 * D + M + Mdot) * E;
  Sigmar += -1897     * cos(4 * D - M - 2 * Mdot) * E;
  Sigmar += -2117     * cos(2 * M - Mdot) * E2;
  Sigmar +=  2354     * cos(2 * D + 2 * M - Mdot) * E2;
  Sigmar += -1423     * cos(4 * D + Mdot);
  Sigmar += -1117     * cos(4 * Mdot);
  Sigmar += -1571     * cos(4 * D - M * E);
  Sigmar += -1739     * cos(D - 2 * Mdot);
  Sigmar += -4421     * cos(2 * Mdot - 2 * F);
  Sigmar +=  1165     * cos(2 * M + Mdot) * E2;
  Sigmar +=  8752     * cos(2 * D - Mdot - 2 * F);

  Sigmab +=  5128122  * sin(F);
  Sigmab +=  280602   * sin(Mdot + F);
  Sigmab +=  277693   * sin(Mdot - F);
  Sigmab +=  173237   * sin(2 * D - F);
  Sigmab +=  55413    * sin(2 * D - Mdot + F);
  Sigmab +=  46271    * sin(2 * D - Mdot - F);
  Sigmab +=  32573    * sin(2 * D + F);
  Sigmab +=  17198    * sin(2 * Mdot + F);
  Sigmab +=  9266     * sin(2 * D + Mdot - F);
  Sigmab +=  8822     * sin(2 * Mdot - F);
  Sigmab +=  8216     * sin(2 * D - M - F) * E;
  Sigmab +=  4326     * sin(2 * D - 2 * Mdot - F);
  Sigmab +=  4200     * sin(2 * D + Mdot +F);
  Sigmab += -3359     * sin(2 * D + M - F) * E;
  Sigmab +=  2463     * sin(2 * D - M - Mdot + F) * E;
  Sigmab +=  2211     * sin(2 * D - M + F) * E;
  Sigmab +=  2065     * sin(2 * D - M - Mdot - F) * E;
  Sigmab += -1870     * sin(M - Mdot - F) * E;
  Sigmab +=  1828     * sin(4 * D - Mdot - F);
  Sigmab += -1794     * sin(M + F) * E;
  Sigmab += -1749     * sin(3 * F);
  Sigmab += -1565     * sin(M - Mdot + F) * E;
  Sigmab += -1491     * sin(D + F);
  Sigmab += -1475     * sin(M + Mdot + F) * E;
  Sigmab += -1410     * sin(M + Mdot - F) * E;
  Sigmab += -1344     * sin(M - F) * E;
  Sigmab += -1335     * sin(D - F);
  Sigmab +=  1107     * sin(3 * Mdot + F);
  Sigmab +=  1021     * sin(4 * D - F);
  Sigmab +=  833      * sin(4 * D - Mdot + F);
  Sigmab +=  777      * sin(Mdot - 3 * F);
  Sigmab +=  671      * sin(4 * D - 2 * Mdot + F);
  Sigmab +=  607      * sin(2 * D - 3 * F);
  Sigmab +=  596      * sin(2 * D + 2 * Mdot - F);
  Sigmab +=  491      * sin(2 * D - M + Mdot - F) * E;
  Sigmab += -451      * sin(2 * D - 2 * Mdot + F);
  Sigmab +=  439      * sin(3 * Mdot - F);
  Sigmab +=  422      * sin(2 * D + 2 * Mdot + F);
  Sigmab +=  421      * sin(2 * D - 3 * Mdot - F);
  Sigmab += -366      * sin(2 * D + M - Mdot + F) * E;
  Sigmab += -351      * sin(2 * D + M + F) * E;
  Sigmab +=  331      * sin(4 * D + F);
  Sigmab +=  315      * sin(2 * D - M + Mdot + F) * E;
  Sigmab +=  302      * sin(2 * D - 2 * M - F) * E2;
  Sigmab += -283      * sin(Mdot + 3 * F);
  Sigmab += -229      * sin(2 * D + M + Mdot - F) * E;
  Sigmab +=  223      * sin(D + M - F);
  Sigmab +=  223      * sin(D + M + F);
  Sigmab += -220      * sin(M - 2 * Mdot - F) * E;
  Sigmab += -220      * sin(2 * D + M - Mdot - F) * E;
  Sigmab += -185      * sin(D + Mdot + F);
  Sigmab +=  181      * sin(2 * D - M - 2 * Mdot - F) * E;
  Sigmab += -177      * sin(M + 2 * Mdot + F) * E;
  Sigmab +=  176      * sin(4 * D - 2 * Mdot - F);
  Sigmab +=  166      * sin(4 * D - M - Mdot - F) * E;
  Sigmab += -164      * sin(D + Mdot - F);
  Sigmab +=  132      * sin(4 * D + Mdot - F);
  Sigmab += -119      * sin(D - Mdot - F);
  Sigmab +=  115      * sin(4 * D - M - F) * E;
  Sigmab +=  107      * sin(2 * D - 2 * M + F) * E2;

  Sigmab += -2235 * sin(Ldot);
  Sigmab +=  382  * sin(A3);
  Sigmab +=  175  * sin(A1 - F);
  Sigmab +=  175  * sin(A1 + F);
  Sigmab +=  127  * sin(Ldot - Mdot);
  Sigmab += -115  * sin(Ldot + Mdot);

  printf("El:         %lf\n", Sigmal); // -1127527
  printf("Eb:         %lf\n", Sigmab); // -3229127
  printf("Er:         %lf\n", Sigmar); // -16590675

  // Right ascension, declination and range
  double ra     = Ldot * RAD2DEG + Sigmal / 1000000;
  printf("ra:         %lf\n", fmod(ra, 360));
  double beta   = Sigmab / 1000000;
  printf("b:          %lf\n", fmod(beta, 360));
  double Delta  = 385000.56 + Sigmar / 1000;
  printf("d:          %lf\n", Delta);

  // Apparent longitude
  double alambda  = ra + dpsi * RAD2DEG;
  printf("alambda:    %lf\n", fmod(alambda, 360));
  alambda *= DEG2RAD;

  if (lambda != NULL)
  {
    *lambda = alambda;
  }

  // Mean oliquity of the ecliptic
  double epsilono = 23.43929 - 0.01300417 * T - 1.638889e-7 * T2
                  + 5.036111e-7 * T3;

  // True obliquity of the ecliptic
  double epsilon = epsilono + depsilon * RAD2DEG;
  printf("epsilon:    %lf\n", fmod(epsilon, 360));

  // Apparent right ascension, declination and radius vector
  redecrv.ra  = atan2(sin(alambda) * cos(epsilon * DEG2RAD)
                    - tan(beta * DEG2RAD) * sin(epsilon * DEG2RAD),
                      cos(alambda));
  printf("GC RA:      %lf\n", redecrv.ra * RAD2DEG);
  redecrv.dec = asin(sin(beta * DEG2RAD) * cos(epsilon * DEG2RAD)
              + cos(beta * DEG2RAD) * sin(epsilon * DEG2RAD)
              * sin(alambda));
  printf("GC Dec:     %lf\n", redecrv.dec * RAD2DEG);
  redecrv.rv  = Delta;

  return redecrv;
}
