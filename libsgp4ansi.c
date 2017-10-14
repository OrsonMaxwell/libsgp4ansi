/*
 * libsgp4ansi.c - an ANSI C-11 SGP4/SDP4 implementation library for sgp4ansid.
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

#include <ctype.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "libsgp4ansi.h"
#include "epoch.h"
#include "const.h"
#include "coord.h"
#include "vector.h"
#include "solar.h"

// ************************************************************************* //
//                                VERSION                                    //
// ************************************************************************* //

const char libsgp4ansi_version[] = VERSION;

// ************************************************************************* //
//                            PRIVATE PROTOTYPES                             //
// ************************************************************************* //

// Trim leading and trailing whitespaces from a string
static size_t
str_trim
(
  const char*  str,
        size_t len,
        char*  out
);

// Recursively zero in on an AOS or LOS event down to 1 sec resolution
time_t
find_aos_los
(
        sat*         s,
  const vec3*        obs_geo,
  const time_t       start_time,
        unsigned int delta_t,
        double       horizon,
        bool         is_aos
);

// ************************************************************************* //
//                             PRIVATE FUNCTIONS                             //
// ************************************************************************* //

/*
 * Trim leading and trailing whitespaces from a string
 *
 * Inputs:  len   - Length of the input in interest
 *          str   - String to trim
 * Outputs: out   - Trimmed and null terminated string, truncated to size
 * Returns: Any   - Length of resulting string
 *         -1     - Failure to read TLE lines
 */
size_t
str_trim
(
  const char*  str,
        size_t len,
        char*  out
)
{
  if(len == 0)
  {
    return 0;
  }

  const char *end;
  size_t out_size;

  // Trim leading space
  while(isspace((unsigned char)* str)) str++;

  if(*str == 0)  // Whitespace only
  {
    *out = 0;
    return 1;
  }

  // Trim trailing spaces
  end = str + strlen(str) - 1;
  while(end > str && isspace((unsigned char)* end)) end--;
  end++;

  // Set output size to minimum of trimmed string length and buffer size minus 1
  out_size = ((end - str) < len - 1) ? (end - str) : len - 1;

  // Copy trimmed string to out
  memcpy(out, str, out_size);
  out[out_size] = 0; // Terminate

  return out_size;
}

// ************************************************************************* //
//                                    API                                    //
// ************************************************************************* //

/*
 * Initialize SGP4/SDP4 orbit model from a raw NORAD TLE lines
 *
 * Inputs:  tlestr1 - TLE line 1
 *          tlestr2 - TLE line 2
 * Outputs: sat     - orbit struct pointer with full orbital data
 * Returns: 0       - Success
 *         -1       - Failure to read TLE lines
 *
 * Calls: orbit_init
 */
int
sat_load_tle
(
  const char* tlestr0,
  const char* tlestr1,
  const char* tlestr2,
  sat* s
)
{
  char str0[130];
  char str1[130];
  char str2[130];

  strcpy(str0, tlestr0);
  strcpy(str1, tlestr1);
  strcpy(str2, tlestr2);

  // Pre-formatting the raw TLE input
  for (char j = 10; j <= 15; j++)
    if (str1[j] == ' ')
      str1[j] = '_';

  if (str1[44] != ' ')
    str1[43] = str1[44];

  str1[44] = '.';
  if (str1[7] == ' ')
    str1[7] = 'U';

  if (str1[9] == ' ')
    str1[9] = '.';

  for (char j = 45; j <= 49; j++)
    if (str1[j] == ' ')
      str1[j] = '0';

  if (str1[51] == ' ')
    str1[51] = '0';

  if (str1[53] != ' ')
    str1[52] = str1[53];

  str1[53] = '.';
  str2[25] = '.';

  for (char j = 26; j <= 32; j++)
    if (str2[j] == ' ')
      str2[j] = '0';

  if (str1[62] == ' ')
    str1[62] = '0';

  if (str1[68] == ' ')
    str1[68] = '0';

  str_trim(str0, 24, s->name);

  struct tm epoch_tm, prop_tm;

  int cardnum, epochyr, nexp, Bexp, checksum, ephem_type, elset_number;
  double nddot, Bstar, epochdays;

  int retval = sscanf(str1,"%2d %5u %1c %8s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6d ",
         &cardnum, &s->norad_number, &s->sec_class, s->int_designator, &epochyr,
         &epochdays,&s->mean_motion_dt2, &nddot, &nexp, &Bstar,
         &Bexp, &ephem_type, &elset_number);

  if (retval != 13)
  {
    return -1;
  }

  if (str2[52] == ' ') // check for minus sign
    retval = sscanf(str2,"%2d %5u %9lf %9lf %8lf %9lf %9lf %10lf %6d \n",
           &cardnum,&s->norad_number, &s->inclination, &s->right_asc_node, &s->eccentricity, &s->argument_perigee,
           &s->mean_anomaly, &s->mean_motion, &s->orbit_number);
  else
    retval = sscanf(str2,"%2d %5u %9lf %9lf %8lf %9lf %9lf %11lf %6d \n",
           &cardnum,&s->norad_number, &s->inclination, &s->right_asc_node, &s->eccentricity, &s->argument_perigee,
           &s->mean_anomaly, &s->mean_motion, &s->orbit_number);

  if (retval != 9)
  {
    return -1;
  }

  if (fractday2unix(epochyr, epochdays, &s->epoch, &s->epoch_ms) != 0)
  {
    return -1;
  }

  s->julian_epoch = unix2jul(s->epoch, s->epoch_ms);

  s->mean_motion_ddt6 = nddot * pow(10, nexp);
  s->Bstar = Bstar * pow(10, Bexp);

  // Converting from TLE to SGP4 units (minutes, radians and kilometres)
  s->inclination      *= DEG2RAD;
  s->right_asc_node   *= DEG2RAD;
  s->argument_perigee *= DEG2RAD;
  s->mean_anomaly     *= DEG2RAD;
  s->mean_motion      /= RPD2RADPM;
  s->mean_motion_dt2  /= (RPD2RADPM * 1440);
  s->mean_motion_ddt6 /= (RPD2RADPM * 1440 * 1440);

  return sat_init(s);
}

/*
 * Initialize SGP4/SDP4 orbit model from a raw NORAD TLE parametres
 *
 * Inputs:  tlestr1 - TLE line 1
 *          tlestr2 - TLE line 2
 * Outputs: sat     - orbit struct pointer with full orbital data
 * Returns: 0       - Success
 *         -1       - Failure to read TLE lines
 *
 * Calls: orbit_init
 */
int
sat_load_params
(
  const char         name[25],
        char         sec_class,
  const char         int_designator[9],
        time_t       epoch,
        float        epoch_ms,
        double       mean_motion_dt2,
        double       mean_motion_ddt6,
        double       Bstar,
        double       inclination,
        double       right_asc_node,
        double       eccentricity,
        double       argument_perigee,
        double       mean_anomaly,
        double       mean_motion,
        unsigned int norad_number,
        unsigned int orbit_number,
        sat*         s
)
{
  strcpy(s->name, name);
  s->sec_class    = sec_class;
  strcpy(s->int_designator, int_designator);

  s->epoch        = epoch;
  s->epoch_ms     = epoch_ms;
  s->julian_epoch = unix2jul(s->epoch, s->epoch_ms);

  // Converting from TLE to SGP4 units (minutes, radians and kilometres)
  s->mean_motion_dt2  = mean_motion_dt2  / (RPD2RADPM * 1440);
  s->mean_motion_ddt6 = mean_motion_ddt6 / (RPD2RADPM * 1440 * 1440);
  s->Bstar            = Bstar;
  s->inclination      = inclination      * DEG2RAD;
  s->right_asc_node   = right_asc_node   * DEG2RAD;
  s->argument_perigee = argument_perigee * DEG2RAD;
  s->mean_anomaly     = mean_anomaly     * DEG2RAD;
  s->mean_motion      = mean_motion      / RPD2RADPM;
  s->norad_number     = norad_number;
  s->orbit_number     = orbit_number;

  return sat_init(s);
}

/*
 * Expand SGP4/SDP4 orbit elements from a sat record containing NORAD TLE data
 *
 * Inputs:  s - sat struct pointer to a record containing NORAD TLE data
 * Outputs: s - sat struct pointer with full orbital data
 * Returns: 0 - Success
 *         -1 - Eccentricity out of range [0, 1)
 *         -2 - Inclination out of range [0, pi]
 */
int
sat_init
(
  sat* s
)
{
#ifdef MATH_TRACE
  printf("---------------------------------------- i1\n");
  printf("mean_motion_dt2  %+.15e\n", s->mean_motion_dt2);
  printf("mean_motion_ddt6 %+.15e\n", s->mean_motion_ddt6);
  printf("Bstar            %+.15e\n", s->Bstar);
  printf("inclination      %+.15e\n", s->inclination);
  printf("right_asc_node   %+.15e\n", s->right_asc_node);
  printf("eccentricity     %+.15e\n", s->eccentricity);
  printf("argument_perigee %+.15e\n", s->argument_perigee);
  printf("mean_anomaly     %+.15e\n", s->mean_anomaly);
  printf("mean_motion      %+.15e\n", s->mean_motion);
#endif

  if ((s->eccentricity < 0) || (s->eccentricity > 0.999999))
  {
    return -1;
  }
  if ((s->inclination < 0) || (s->inclination > PI))
  {
    return -2;
  }

  // Recovering mean motion and semi-major axis
  double aa      = pow(XKE / s->mean_motion, TWOTHIRD);
  double cosio   = cos(s->inclination);
  double sinio   = sin(s->inclination);
  double eo2     = pow(s->eccentricity, 2);
  double theta2  = pow(cosio, 2);
  double x3th2m1 = 3 * theta2 - 1;
      s->x1mth2  = 1.0 - theta2;
      s->x1m5th2 = 1 - 5 * theta2;
      s->con41   = -s->x1m5th2 - 2 * theta2;
      s->x7thm1  = 7 * theta2 - 1;
  double betao2  = 1 - eo2;
  double betao   = sqrt(betao2);

  // Un-Kozai mean motion
  double delta1  = 0.75 * J2 * (3.0 * theta2 - 1.0) / (betao * betao2);
  double a0      = aa * (1.0 - pow(delta1 / pow(aa, 2), 2)
                 - delta1 / pow(aa, 2) * (1.0 / 3.0 + 134.0
                 * pow(delta1 / pow(aa, 2), 2) / 81.0));
  double delta0  = delta1 / pow(a0, 2);
  s->xnodp       = s->mean_motion / (1 + delta0);
  s->aodp        = pow(XKE / s->xnodp, TWOTHIRD);
  s->perigee_alt = (s->aodp * (1 - s->eccentricity)) * RE - RE;
  s->period      = TAU / s->xnodp;

#ifdef MATH_TRACE
  printf("---------------------------------------- i2\n");
  printf("aa      %+.15e\n", aa);
  printf("cosio   %+.15e\n", cosio);
  printf("sinio   %+.15e\n", sinio);
  printf("eo2     %+.15e\n", eo2);
  printf("theta2  %+.15e\n", theta2);
  printf("betao2  %+.15e\n", betao2);
  printf("betao   %+.15e\n", betao);
  printf("delta1  %+.15e\n", delta1);
  printf("a0      %+.15e\n", a0);
  printf("delta0  %+.15e\n", delta0);
  printf("xnodp   %+.15e\n", s->xnodp);
  printf("aodp    %+.15e\n", s->aodp);
  printf("x1m5th2 %+.15e\n", s->x1m5th2);
  printf("con41   %+.15e\n", s->con41);
#endif

  if (s->period >= 225)
  {
    s->is_deep_space      = true;
    s->use_simple_model   = true;
  }
  else
  {
    s->is_deep_space      = false;
    s->use_simple_model   = false;

    if (s->perigee_alt < 220)
    {
      s->use_simple_model = true;
    }
  }

  double qoms24 = pow(42 / RE, 4);
  double sfour  = S;

  // Adjust sfour and qoms4 for different perigee altitudes
  if (s->perigee_alt < 156)
  {
    sfour   = s->perigee_alt - 78;
    if (s->perigee_alt < 98)
    {
      sfour = 20;
    }

    qoms24  = pow((120 - sfour) / RE, 4);
    sfour   = sfour / RE + 1; // Get back to AE units
  }

  // Compute common constants
  double pinv2   = 1 / (pow(s->aodp, 2) * pow(betao2, 2));
  double tsi     = 1 / (s->aodp - sfour);
      s->eta     = s->aodp * s->eccentricity * tsi;
  double eta2    = pow(s->eta, 2);
  double eeta    = s->eccentricity * s->eta;
  double psi2    = fabs(1 - eta2);
  double coef    = qoms24 * pow(tsi, 4);
  double coef1   = coef / pow(psi2, 3.5);
  double C2      = coef1 * s->xnodp
                 * (s->aodp * (1 + 1.5 * eta2 + eeta * (4 + eta2))
                 + 0.375 * J2 * tsi / psi2 * x3th2m1
                 * (8 + 3 * eta2 * (8 + eta2)));
      s->C1      = s->Bstar * C2;
      s->C4      = 2 * s->xnodp * coef1 * s->aodp * betao2
                 * (s->eta * (2 + 0.5 * eta2) + s->eccentricity
                 * (0.5 + 2 * eta2) - J2 * tsi / (s->aodp * psi2)
                 * (-3 * x3th2m1 * (1 - 2 * eeta + eta2
                 * (1.5 - 0.5 * eeta)) + 0.75 * s->x1mth2
                 * (2 * eta2 - eeta * (1 + eta2))
                 * cos(2 * s->argument_perigee)));
  double theta4  = pow(cosio, 4);
  double temp1   = 1.5 * J2 * pinv2 * s->xnodp;
  double temp2   = 0.5 * temp1 * J2 * pinv2;
  double temp3   = -0.46875 * J4 * pow(pinv2, 2) * s->xnodp;
      s->xmdot   = s->xnodp + 0.5 * temp1 * betao * x3th2m1 + 0.0625
                 * temp2 * betao * (13 - 78 * theta2 + 137 * theta4);
      s->omgdot  = -0.5 * temp1 * s->x1m5th2
                 + 0.0625 * temp2 * (7 - 114 * theta2 + 395 * theta4)
                 + temp3 * (3 - 36 * theta2 + 49 * theta4);
  double xhdot1  = -temp1 * cosio;
      s->xnodot  = xhdot1 + (0.5 * temp2 * (4 - 19 * theta2)
                 + 2 * temp3 * (3 - 7 * theta2)) * cosio;
      s->xnodcf  = 3.5 * betao2 * xhdot1 * s->C1;
      s->t2cof   = 1.5 * s->C1;
  // Division by zero check then inclination = 180 deg
      s->xlcof   = 0.125 * A3OVK2 * sinio * (3.0 + 5.0 * cosio)
                 / ((fabs(cosio+1.0) > 1.5e-12) ? ((1.0 + cosio)) : (1.5e-12));
      s->aycof   = 0.25 * A3OVK2 * sinio;

#ifdef MATH_TRACE
  printf("---------------------------------------- i3\n");
  printf("qoms24 %+.15e\n", qoms24);
  printf("s4     %+.15e\n", sfour);
  printf("pinv2  %+.15e\n", pinv2);
  printf("tsi    %+.15e\n", tsi);
  printf("eta    %+.15e\n", s->eta);
  printf("eta2   %+.15e\n", eta2);
  printf("eeta   %+.15e\n", eeta);
  printf("psi2   %+.15e\n", psi2);
  printf("coef   %+.15e\n", coef);
  printf("coef1  %+.15e\n", coef1);
  printf("C2     %+.15e\n", C2);
  printf("C1     %+.15e\n", s->C1);
  printf("C4     %+.15e\n", s->C4);
  printf("theta4 %+.15e\n", theta4);
  printf("xmdot  %+.15e\n", s->xmdot);
  printf("omgdot %+.15e\n", s->omgdot);
  printf("xhdot1 %+.15e\n", xhdot1);
  printf("xnodot %+.15e\n", s->xnodot);
  printf("xnodcf %+.15e\n", s->xnodcf);
  printf("xlcof  %+.15e\n", s->xlcof);
  printf("aycof  %+.15e\n", s->aycof);
  printf("x7thm1 %+.15e\n", s->x7thm1);
  printf("x1mth2 %+.15e\n", s->x1mth2);
  printf("dspace %d\n", s->is_deep_space);
  printf("simple %d\n", s->use_simple_model);
#endif

  if (s->is_deep_space == true) // Deep space init here
  {

    s->GSTo  = jul2gst(s->julian_epoch);

    // Constants
    const double zes    =  0.01675;
    const double zel    =  0.05490;
    const double c1ss   =  2.9864797e-6;
    const double c1l    =  4.7968065e-7;
    const double zsinis =  0.39785416;
    const double zcosis =  0.91744867;
    const double zcosgs =  0.1945905;
    const double zsings = -0.98088458;

    double sinomm = sin(s->argument_perigee);
    double cosomm = cos(s->argument_perigee);
    double sinim  = sin(s->inclination);
    double cosim  = cos(s->inclination);

    // Initialize lunar solar terms
    s->peo    = 0;
    s->pinco  = 0;
    s->plo    = 0;
    s->pgho   = 0;
    s->pho    = 0;

    double day    = s->julian_epoch + 18261.5 - B1950;
    double xnodce = fmod(4.5236020 - 9.2422029e-4 * day, TAU);
    double stem   = sin(xnodce);
    double ctem   = cos(xnodce);
    double zcosil = 0.91375164 - 0.03568096 * ctem;
    double zsinil = sqrt(1.0 - zcosil * zcosil);
    double zsinhl = 0.089683511 * stem / zsinil;
    double zcoshl = sqrt(1.0 - zsinhl * zsinhl);
    double gam    = 5.8351514 + 0.0019443680 * day;
    double zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    double zx     = gam + atan2(0.39785416 * stem / zsinil, zy) - xnodce;
    double zcosgl = cos(zx);
    double zsingl = sin(zx);

#ifdef MATH_TRACE
    printf("======================================== id1\n");
    printf("[DS] julian %+.15e\n", s->julian_epoch);
    printf("[DS] day    %+.15e\n", day);
    printf("[DS] xnodce %+.15e\n", xnodce);
    printf("[DS] stem   %+.15e\n", stem);
    printf("[DS] ctem   %+.15e\n", ctem);
    printf("[DS] zcosil %+.15e\n", zcosil);
    printf("[DS] zsinil %+.15e\n", zsinil);
    printf("[DS] zsinhl %+.15e\n", zsinhl);
    printf("[DS] zcoshl %+.15e\n", zcoshl);
    printf("[DS] gam    %+.15e\n", gam);
    printf("[DS] zy     %+.15e\n", zy);
    printf("[DS] zx     %+.15e\n", zx);
    printf("[DS] zcosgl %+.15e\n", zcosgl);
    printf("[DS] zsingl %+.15e\n", zsingl);
#endif

    s->zmos = fmod(6.2565837 + 0.017201977 * day, TAU);

    // Do solar terms
    double cosq  = cos(s->right_asc_node);
    double sinq  = sin(s->right_asc_node);
    double zcosg = zcosgs;
    double zsing = zsings;
    double zcosi = zcosis;
    double zsini = zsinis;
    double zcosh = cosq;
    double zsinh = sinq;
    double cc    = c1ss;
    double xnoi  = 1.0 / s->xnodp;

    // Iterative terms
    double a1,  a2,  a3,  a4,   a5,   a6,   a7,   a8,   a9,   a10;
    double s1,  s2,  s3,  s4,   s5,   s6,   s7;
    double ss1, ss2, ss3, ss4,  ss5,  ss6,  ss7;
    double sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33;
    double x1,  x2,  x3,  x4,   x5,   x6,   x7,   x8;
    double z1,  z2,  z3,  z11,  z12,  z13,  z21,  z22,  z23,  z31,  z32,  z33;

    for (uint8_t ls_flag = 1; ls_flag <= 2; ls_flag++)
    {
      a1  =  zcosg * zcosh + zsing * zcosi * zsinh;
      a3  = -zsing * zcosh + zcosg * zcosi * zsinh;
      a7  = -zcosg * zsinh + zsing * zcosi * zcosh;
      a8  =  zsing * zsini;
      a9  =  zsing * zsinh + zcosg * zcosi * zcosh;
      a10 =  zcosg * zsini;
      a2  =  cosim * a7 + sinim * a8;
      a4  =  cosim * a9 + sinim * a10;
      a5  = -sinim * a7 + cosim * a8;
      a6  = -sinim * a9 + cosim * a10;

      x1  =  a1 * cosomm + a2 * sinomm;
      x2  =  a3 * cosomm + a4 * sinomm;
      x3  = -a1 * sinomm + a2 * cosomm;
      x4  = -a3 * sinomm + a4 * cosomm;
      x5  =  a5 * sinomm;
      x6  =  a6 * sinomm;
      x7  =  a5 * cosomm;
      x8  =  a6 * cosomm;

      z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
      z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
      z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
      z1  =  3.0 *  (a1 * a1 + a2 * a2) + z31 * eo2;
      z2  =  6.0 *  (a1 * a3 + a2 * a4) + z32 * eo2;
      z3  =  3.0 *  (a3 * a3 + a4 * a4) + z33 * eo2;
      z11 = -6.0 * a1 * a5 + eo2 *  (-24.0 * x1 * x7-6.0 * x3 * x5);
      z12 = -6.0 *  (a1 * a6 + a3 * a5) + eo2 *
          (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
      z13 = -6.0 * a3 * a6 + eo2 * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
      z21 =  6.0 * a2 * a5 + eo2 * (24.0 * x1 * x5 - 6.0 * x3 * x7);
      z22 =  6.0 *  (a4 * a5 + a2 * a6) + eo2 *
          (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
      z23 =  6.0 * a4 * a6 + eo2 * (24.0 * x2 * x6 - 6.0 * x4 * x8);
      z1  = z1 + z1 + betao2 * z31;
      z2  = z2 + z2 + betao2 * z32;
      z3  = z3 + z3 + betao2 * z33;
      s3  = cc * xnoi;
      s2  = -0.5 * s3 / betao;
      s4  = s3 * betao;
      s1  = -15.0 * s->eccentricity * s4;
      s5  = x1 * x3 + x2 * x4;
      s6  = x2 * x3 + x1 * x4;
      s7  = x2 * x4 - x1 * x3;

      // Do lunar terms
      if (ls_flag == 1)
      {
        ss1   = s1;
        ss2   = s2;
        ss3   = s3;
        ss4   = s4;
        ss5   = s5;
        ss6   = s6;
        ss7   = s7;
        sz1   = z1;
        sz2   = z2;
        sz3   = z3;
        sz11  = z11;
        sz12  = z12;
        sz13  = z13;
        sz21  = z21;
        sz22  = z22;
        sz23  = z23;
        sz31  = z31;
        sz32  = z32;
        sz33  = z33;
        zcosg = zcosgl;
        zsing = zsingl;
        zcosi = zcosil;
        zsini = zsinil;
        zcosh = zcoshl * cosq + zsinhl * sinq;
        zsinh = sinq * zcoshl - cosq * zsinhl;
        cc    = c1l;
      }
    }

    s->zmol = fmod(4.7199672 + 0.22997150  * day - gam, TAU);

#ifdef MATH_TRACE
    printf("======================================== id2\n");
    printf("[DS] xnoi  %+.15e\n", xnoi);
    printf("[DS] a1    %+.15e\n", a1);
    printf("[DS] a2    %+.15e\n", a2);
    printf("[DS] a3    %+.15e\n", a3);
    printf("[DS] a4    %+.15e\n", a4);
    printf("[DS] a5    %+.15e\n", a5);
    printf("[DS] a6    %+.15e\n", a6);
    printf("[DS] a7    %+.15e\n", a7);
    printf("[DS] a8    %+.15e\n", a8);
    printf("[DS] a9    %+.15e\n", a9);
    printf("[DS] a10   %+.15e\n", a10);
    printf("[DS] x1    %+.15e\n", x1);
    printf("[DS] x2    %+.15e\n", x2);
    printf("[DS] x3    %+.15e\n", x3);
    printf("[DS] x4    %+.15e\n", x4);
    printf("[DS] x5    %+.15e\n", x5);
    printf("[DS] x6    %+.15e\n", x6);
    printf("[DS] x7    %+.15e\n", x7);
    printf("[DS] x8    %+.15e\n", x8);
    printf("[DS] z1    %+.15e\n", z1);
    printf("[DS] z2    %+.15e\n", z2);
    printf("[DS] z3    %+.15e\n", z3);
    printf("[DS] z11   %+.15e\n", z11);
    printf("[DS] z12   %+.15e\n", z12);
    printf("[DS] z13   %+.15e\n", z13);
    printf("[DS] z21   %+.15e\n", z21);
    printf("[DS] z22   %+.15e\n", z22);
    printf("[DS] z23   %+.15e\n", z23);
    printf("[DS] z31   %+.15e\n", z31);
    printf("[DS] z32   %+.15e\n", z32);
    printf("[DS] z33   %+.15e\n", z33);
    printf("[DS] s1    %+.15e\n", s1);
    printf("[DS] s2    %+.15e\n", s2);
    printf("[DS] s3    %+.15e\n", s3);
    printf("[DS] s4    %+.15e\n", s4);
    printf("[DS] s5    %+.15e\n", s5);
    printf("[DS] s6    %+.15e\n", s6);
    printf("[DS] s7    %+.15e\n", s7);
    printf("[DS] ss1   %+.15e\n", ss1);
    printf("[DS] ss2   %+.15e\n", ss2);
    printf("[DS] ss3   %+.15e\n", ss3);
    printf("[DS] ss4   %+.15e\n", ss4);
    printf("[DS] ss5   %+.15e\n", ss5);
    printf("[DS] ss6   %+.15e\n", ss6);
    printf("[DS] ss7   %+.15e\n", ss7);
    printf("[DS] sz11  %+.15e\n", sz11);
    printf("[DS] sz12  %+.15e\n", sz12);
    printf("[DS] sz13  %+.15e\n", sz13);
    printf("[DS] sz21  %+.15e\n", sz21);
    printf("[DS] sz22  %+.15e\n", sz22);
    printf("[DS] sz23  %+.15e\n", sz23);
    printf("[DS] sz31  %+.15e\n", sz31);
    printf("[DS] sz32  %+.15e\n", sz32);
    printf("[DS] sz33  %+.15e\n", sz33);
    printf("[DS] zcosh %+.15e\n", zcosh);
    printf("[DS] zsinh %+.15e\n", zsinh);
    printf("[DS] zmos  %+.15e\n", s->zmos);
    printf("[DS] zmol  %+.15e\n", s->zmol);
#endif

    // Do final solar terms
    s->se2  =   2 * ss1 * ss6;
    s->se3  =   2 * ss1 * ss7;
    s->si2  =   2 * ss2 * sz12;
    s->si3  =   2 * ss2 * (sz13 - sz11);
    s->sl2  =  -2 * ss3 * sz2;
    s->sl3  =  -2 * ss3 * (sz3 - sz1);
    s->sl4  =  -2 * ss3 * (-21 - 9 * eo2) * zes;
    s->sgh2 =   2 * ss4 * sz32;
    s->sgh3 =   2 * ss4 * (sz33 - sz31);
    s->sgh4 = -18 * ss4 * zes;
    s->sh2  =  -2 * ss2 * sz22;
    s->sh3  =  -2 * ss2 * (sz23 - sz21);

#ifdef MATH_TRACE
    printf("======================================== id3\n");
    printf("[DS] se2  %+.15e\n", s->se2);
    printf("[DS] se3  %+.15e\n", s->se3);
    printf("[DS] si2  %+.15e\n", s->si2);
    printf("[DS] si3  %+.15e\n", s->si3);
    printf("[DS] sl2  %+.15e\n", s->sl2);
    printf("[DS] sl3  %+.15e\n", s->sl3);
    printf("[DS] sl4  %+.15e\n", s->sl4);
    printf("[DS] sgh2 %+.15e\n", s->sgh2);
    printf("[DS] sgh3 %+.15e\n", s->sgh3);
    printf("[DS] sgh4 %+.15e\n", s->sgh4);
    printf("[DS] sh2  %+.15e\n", s->sh2);
    printf("[DS] sh3  %+.15e\n", s->sh3);
#endif

    // Do final lunar terms
    s->ee2  =   2 * s1 * s6;
    s->e3   =   2 * s1 * s7;
    s->xi2  =   2 * s2 * z12;
    s->xi3  =   2 * s2 * (z13 - z11);
    s->xl2  =  -2 * s3 * z2;
    s->xl3  =  -2 * s3 * (z3 - z1);
    s->xl4  =  -2 * s3 * (-21 - 9 * eo2) * zel;
    s->xgh2 =   2 * s4 * z32;
    s->xgh3 =   2 * s4 * (z33 - z31);
    s->xgh4 = -18 * s4 * zel;
    s->xh2  =  -2 * s2 * z22;
    s->xh3  =  -2 * s2 * (z23 - z21);

#ifdef MATH_TRACE
    printf("======================================== id4\n");
    printf("[DS] ee2  %+.15e\n", s->ee2);
    printf("[DS] e3   %+.15e\n", s->e3);
    printf("[DS] xi2  %+.15e\n", s->xi2);
    printf("[DS] xi3  %+.15e\n", s->xi3);
    printf("[DS] xl2  %+.15e\n", s->xl2);
    printf("[DS] xl3  %+.15e\n", s->xl3);
    printf("[DS] xl4  %+.15e\n", s->xl4);
    printf("[DS] xgh2 %+.15e\n", s->xgh2);
    printf("[DS] xgh3 %+.15e\n", s->xgh3);
    printf("[DS] xgh4 %+.15e\n", s->xgh4);
    printf("[DS] xh2  %+.15e\n", s->xh2);
    printf("[DS] xh3  %+.15e\n", s->xh3);
    printf("[DS] GSTo %+.15e\n", s->GSTo);
#endif

    const double q22    = 1.7891679e-6;
    const double q31    = 2.1460748e-6;
    const double q33    = 2.2123015e-7;
    const double root22 = 1.7891679e-6;
    const double root44 = 7.3636953e-9;
    const double root54 = 2.1765803e-9;
    const double root32 = 3.7393792e-7;
    const double root52 = 1.1428639e-7;

    // Deep space resonance initialization
    s->is_12h_resonant = false;
    s->is_24h_resonant = false;

    if ((s->mean_motion < 0.0052359877) && (s->mean_motion > 0.0034906585))
    {
      s->is_24h_resonant = true; // 24h resonance
    }
    if ((s->mean_motion >= 8.26e-3)
        && (s->mean_motion <= 9.24e-3)
        && (s->eccentricity >= 0.5))
    {
      s->is_12h_resonant = true; // 12h resonance
    }

    // Do solar terms
    double ses  =  ss1 * ZNS * ss5;
    double sis  =  ss2 * ZNS * (sz11 + sz13);
    double sls  = -ZNS * ss3 * (sz1 + sz3 - 14 - 6 * eo2);
    double sghs =  ss4 * ZNS * (sz31 + sz33 - 6);
    double shs  = -ZNS * ss2 * (sz21 + sz23);

    // Fix for 180 deg inclination
    if ((s->inclination < 5.2359877e-2)
        || (s->inclination > PI - 5.2359877e-2))
    {
      shs = 0;
    }

    if (sinim != 0.0)
    {
      shs = shs / sinim;
    }

    double sgs  = sghs - cosim * shs;

    // Do lunar terms
    s->dedt = ses + s1 * ZNL * s5;
    s->didt = sis + s2 * ZNL * (z11 + z13);
    s->dmdt = sls - ZNL * s3 * (z1 + z3 - 14.0 - 6.0 * eo2);

    double sghl = s4 * ZNL * (z31 + z33 - 6.0);
    double shll = -ZNL * s2 * (z21 + z23);

    // sgp4fix for 180 deg incl
    if ((s->inclination < 5.2359877e-2)
        || (s->inclination > PI - 5.2359877e-2))
    {
      shll = 0.0;
    }

    s->domdt = sgs + sghl;
    s->dnodt = shs;

    if (sinim != 0.0)
    {
      s->domdt -= cosim / sinim * shll;
      s->dnodt += shll / sinim;
    }

#ifdef MATH_TRACE
    printf("======================================== id5\n");
    printf("[DS] ses   %+.15e\n", ses);
    printf("[DS] sis   %+.15e\n", sis);
    printf("[DS] sls   %+.15e\n", sls);
    printf("[DS] sghs  %+.15e\n", sghs);
    printf("[DS] shs   %+.15e\n", shs);
    printf("[DS] sgs   %+.15e\n", sgs);
    printf("[DS] dedt  %+.15e\n", s->dedt);
    printf("[DS] didt  %+.15e\n", s->didt);
    printf("[DS] dmdt  %+.15e\n", s->dmdt);
    printf("[DS] domdt %+.15e\n", s->domdt);
    printf("[DS] dnodt %+.15e\n", s->dnodt);
    printf("[DS] sghl  %+.15e\n", sghl);
    printf("[DS] shll  %+.15e\n", shll);
#endif

    // Calculate deep space resonance effects
    double aonv  = pow(s->xnodp / XKE, TWOTHIRD);
    double ainv2 = pow(aonv, 2);

#ifdef MATH_TRACE
    printf("======================================== id6\n");
    printf("[DS] 12h res %d\n", (int)s->is_12h_resonant);
    printf("[DS] 24h res %d\n", (int)s->is_24h_resonant);
#endif

    double sini2  =  sinim * sinim;
    double cosisq =  cosim * cosim;

    // Initialize the resonance terms
    if ((s->is_12h_resonant == true)
        || (s->is_24h_resonant == true))
    {
      // Geopotential resonance for 12 hour orbits
      if (s->is_12h_resonant == true)
      {
        double eocu   = s->eccentricity * s->eccentricity * s->eccentricity;
        double g201   = -0.306 - (s->eccentricity - 0.64) * 0.440;

        // 12h Resonant polynomials
        double g211, g310, g322, g410, g422, g520, g521, g523, g532, g533;

        if (s->eccentricity <= 0.65)
        {
          g211 =    3.616  -  13.2470 * s->eccentricity
               +   16.2900 * eo2;
          g310 =  -19.302  + 117.3900 * s->eccentricity
               -  228.4190 * eo2 +  156.5910 * eocu;
          g322 =  -18.9068 + 109.7927 * s->eccentricity
               -  214.6334 * eo2 +  146.5816 * eocu;
          g410 =  -41.122  + 242.6940 * s->eccentricity
               -  471.0940 * eo2 +  313.9530 * eocu;
          g422 = -146.407  + 841.8800 * s->eccentricity
               - 1629.014  * eo2 + 1083.4350 * eocu;
          g520 = -532.114  + 3017.977 * s->eccentricity
               - 5740.032  * eo2 + 3708.2760 * eocu;
        }
        else
        {
          g211 =   -72.099 +   331.819 * s->eccentricity
               -   508.738 * eo2 +   266.724 * eocu;
          g310 =  -346.844 +  1582.851 * s->eccentricity
               -  2415.925 * eo2 +  1246.113 * eocu;
          g322 =  -342.585 +  1554.908 * s->eccentricity
               -  2366.899 * eo2 +  1215.972 * eocu;
          g410 = -1052.797 +  4758.686 * s->eccentricity
               -  7193.992 * eo2 +  3651.957 * eocu;
          g422 = -3581.690 + 16178.110 * s->eccentricity
               - 24462.770 * eo2 + 12422.520 * eocu;

          if (s->eccentricity > 0.715)
          {
            g520 = -5149.66 + 29936.92 * s->eccentricity
                 - 54087.36 * eo2 + 31324.56 * eocu;
          }
          else
          {
            g520 = 1464.74 - 4664.75 * s->eccentricity + 3763.64 * eo2;
          }
        }
        if (s->eccentricity < 0.7)
        {
          g533 = -919.22770 + 4988.6100 * s->eccentricity
               - 9064.7700  * eo2 + 5542.21  * eocu;
          g521 = -822.71072 + 4568.6173 * s->eccentricity
               - 8491.4146  * eo2 + 5337.524 * eocu;
          g532 = -853.66600 + 4690.2500 * s->eccentricity
               - 8624.7700  * eo2 + 5341.4  * eocu;
        }
        else
        {
          g533 = -37995.780 + 161616.52 * s->eccentricity
               - 229838.20  * eo2 + 109377.94 * eocu;
          g521 = -51752.104 + 218913.95 * s->eccentricity
               - 309468.16  * eo2 + 146349.42 * eocu;
          g532 = -40023.880 + 170470.89 * s->eccentricity
               - 242699.48  * eo2 + 115605.82 * eocu;
        }

        double f220, f221, f321, f322, f441, f442, f522, f523, f542, f543;

        f220 =  0.75 * (1.0 + 2.0 * cosim+cosisq);
        f221 =  1.5 * sini2;
        f321 =  1.875 * sinim  *  (1.0 - 2.0 * cosim - 3.0 * cosisq);
        f322 = -1.875 * sinim  *  (1.0 + 2.0 * cosim - 3.0 * cosisq);
        f441 = 35.0 * sini2 * f220;
        f442 = 39.3750 * sini2 * sini2;
        f522 =  9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim- 5.0 * cosisq)
             + 0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq) );
        f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim
             + 10.0 * cosisq) + 6.56250012 * (1.0+2.0 * cosim - 3.0 * cosisq));
        f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim+cosisq
             * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
        f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim+cosisq
             * (12.0 + 8.0 * cosim - 10.0 * cosisq));

               temp1 = 3.0 * pow(s->xnodp, 2) * ainv2;
        double temp  = temp1 * root22;
            s->d2201 = temp * f220 * g201;
            s->d2211 = temp * f221 * g211;
               temp1 = temp1 * aonv;
               temp  = temp1 * root32;
            s->d3210 = temp * f321 * g310;
            s->d3222 = temp * f322 * g322;
               temp1 = temp1 * aonv;
               temp  = 2.0 * temp1 * root44;
            s->d4410 = temp * f441 * g410;
            s->d4422 = temp * f442 * g422;
               temp1 = temp1 * aonv;
               temp  = temp1 * root52;
            s->d5220 = temp * f522 * g520;
            s->d5232 = temp * f523 * g532;
               temp  = 2.0 * temp1 * root54;
            s->d5421 = temp * f542 * g521;
            s->d5433 = temp * f543 * g533;
            s->xlamo = fmod(s->mean_anomaly + 2 * s->right_asc_node
                     - 2 * s->GSTo, TAU);
            s->xfact = s->xmdot + s->dmdt
                     + 2 * (s->xnodot + s->dnodt - RPTIM) - s->xnodp;
      }

      // Synchronous resonance terms

      double g200, g310, g300, f220, f311, f330;

      if (s->is_24h_resonant == true)
      {
            g200 = 1 + eo2 * (-2.5 + 0.8125 * eo2);
            g310 = 1 + 2 * eo2;
            g300 = 1 + eo2 * (-6 + 6.60937 * eo2);
            f220 = 0.75 * (1 + cosim) * (1 + cosim);
            f311 = 0.9375 * sini2 * (1 + 3 * cosim) - 0.75 * (1 + cosim);
            f330 = 1 + cosim;
            f330 = 1.875 * f330 * f330 * f330;
         s->del1 = 3 * pow(s->xnodp, 2) * ainv2;
         s->del2 = 2 * s->del1 * f220 * g200 * q22;
         s->del3 = 3 * s->del1 * f330 * g300 * q33 * aonv;
         s->del1 = s->del1 * f311 * g310 * q31 * aonv;
        s->xfact = s->xmdot + (s->omgdot + s->xnodot) - RPTIM + s->dmdt
                 + s->domdt + s->dnodt - s->xnodp;
        s->xlamo = fmod(s->mean_anomaly + s->right_asc_node
                 + s->argument_perigee - s->GSTo, TAU);
      }

      // Initialize the integrator
      s->xli    = s->xlamo;
      s->xni    = s->xnodp;
      s->atime  = 0;
    }
  }

  // Compute near space constants
  double C3      = 0;
      s->xmcof   = 0;

  if (s->eccentricity > 1.0e-4)
  {
    C3           = coef * tsi * A3OVK2 * s->xnodp * sinio / s->eccentricity;
    s->xmcof     = -TWOTHIRD * coef * s->Bstar / eeta;
  }
      s->C5      = 2 * coef1 * s->aodp * betao2
                 * (1 + 2.75 * (eta2 + eeta) + eeta * eta2);
      s->omgcof  = s->Bstar * C3 * cos(s->argument_perigee);
      s->delmo   = pow(1 + s->eta * cos(s->mean_anomaly), 3);
      s->sinmo   = sin(s->mean_anomaly);

  if (s->use_simple_model == false)
  {
    double C12   = pow(s->C1, 2);
        s->D2    = 4 * s->aodp * tsi * C12;
    double temp  = s->D2 * tsi * s->C1 / 3;
        s->D3    = (17 * s->aodp + sfour) * temp;
        s->D4    = 0.5 * temp * s->aodp * tsi
                 * (221 * s->aodp + 31 * sfour) * s->C1;
        s->t3cof = s->D2 + 2 * C12;
        s->t4cof = 0.25 * (3 * s->D3 + s->C1 * (12 * s->D2 + 10 * C12));
        s->t5cof = 0.2  * (3 * s->D4 + 12 * s->C1 * s->D3 + 6 * pow(s->D2, 2)
                 + 15 * C12 * (2 * s->D2 + C12));
  }

#ifdef MATH_TRACE
  printf("---------------------------------------- i4\n");
  printf("C3     %+.15e\n", C3);
  printf("xmcof  %+.15e\n", s->xmcof);
  printf("C5     %+.15e\n", s->C5);
  printf("omgcof %+.15e\n", s->omgcof);
  printf("delmo  %+.15e\n", s->delmo);
  printf("sinmo  %+.15e\n", s->sinmo);
  printf("D2     %+.15e\n", s->D2);
  printf("D3     %+.15e\n", s->D3);
  printf("D4     %+.15e\n", s->D4);
  printf("t3cof  %+.15e\n", s->t3cof);
  printf("t4cof  %+.15e\n", s->t4cof);
  printf("t5cof  %+.15e\n", s->t5cof);

  vec3 p, v;
  return sat_propagate(s, 0, 10, 1.0e-12, &p, &v);
#else
  // Propagate at zero time since epoch
  return sat_propagate(s, 0, 10, 1.0e-12, NULL, NULL);
#endif
}

/*
 * Get position and velocity vectors in the TEME frame at given time since epoch
 *
 * Inputs:  s         - sat struct pointer with initialized orbital data
 *          tdelta    - Time since orbit TLE epoch, minutes
 *          maxiter   - Kepler's equation maximum iteration count
 *          tolerance - Kepler's equation desired precision tolerance
 * Outputs: p         - 3D position vector in TEME frame in km
 *          v         - 3D velocity vector in TEME frame in km/sec
 * Returns: 0         - Success
 *         -1         - Invalid inputs or parametres
 *         -2         - Negative mean motion
 *         -3         - Eccentricity out of range (e >= 1; e < -1.0e-12)
 *         -4         - Short period preliminary quantities error
 *         -5         - Decayed satellite
 */
int
sat_propagate
(
    sat* s,
    double tdelta,
    unsigned int maxiter,
    double tolerance,
    vec3* p,
    vec3* v
)
{
  if ((s == NULL) || (maxiter < 1) || (tolerance <= 0.0))
  {
    return -1;
  }

  // Update for secular gravity and atmospheric drag
  double xmdf     = s->mean_anomaly + s->xmdot * tdelta;
  double omgadf   = s->argument_perigee + s->omgdot * tdelta;
  double xnoddf   = s->right_asc_node + s->xnodot * tdelta;
  double t2       = pow(tdelta, 2);
  double xnode    = xnoddf + s->xnodcf * t2;
  double tempa    = 1 - s->C1 * tdelta;
  double tempe    = s->Bstar * s->C4 * tdelta;
  double templ    = s->t2cof * t2;
  double omega    = omgadf;
  double xmp      = xmdf;

  if (s->use_simple_model == false)
  {
    double delomg = s->omgcof * tdelta;
    double delm   = s->xmcof * (pow(1.0 + s->eta * cos(xmdf), 3) - s->delmo);
    xmp           = xmdf + delomg + delm;
    omega         = omgadf - delomg - delm;
    double t3     = t2 * tdelta;
    double t4     = t3 * tdelta;
    tempa         = tempa - s->D2 * t2 - s->D3 * t3 - s->D4 * t4;
    tempe         = tempe + s->Bstar * s->C5 * (sin(xmp) - s->sinmo);
    templ         = templ + s->t3cof * t3 + t4 * (s->t4cof + tdelta * s->t5cof);
  }

  double nm    = s->xnodp;
  double em    = s->eccentricity;
  double inclm = s->inclination;

#ifdef MATH_TRACE
  printf("---------------------------------------- p1\n");
  printf("xmdf   %+.15e\n", xmdf);
  printf("xmp    %+.15e\n", xmp);
  printf("omgadf %+.15e\n", omgadf);
  printf("xnoddf %+.15e\n", xnoddf);
  printf("t2     %+.15e\n", t2);
  printf("tempa  %+.15e\n", tempa);
  printf("tempe  %+.15e\n", tempe);
  printf("templ  %+.15e\n", templ);
  printf("omega  %+.15e\n", omega);
  printf("nm     %+.15e\n", nm);
  printf("em     %+.15e\n", em);
  printf("inclm  %+.15e\n", inclm);
#endif

  if (s->is_deep_space == true)
  {
    // Deep space contributions to mean elements for perturbing third body
    const double fasx2 = 0.13130908;
    const double fasx4 = 2.8843198;
    const double fasx6 = 0.37448087;
    const double g22   = 5.7686396;
    const double g32   = 0.95240898;
    const double g44   = 1.8014998;
    const double g52   = 1.0508330;
    const double g54   = 4.4108898;
    const double stepp = 720;
    const double stepn = -720;
    const double step2 = 259200;

    // Calculate deep space resonance effects
    s->dndt       = 0;
    double theta  = fmod(s->GSTo + tdelta * RPTIM, TAU);

    // Perturbed quantities
    em += s->dedt * tdelta;
    inclm += s->didt  * tdelta;
    omega += s->domdt * tdelta;
    xnode += s->dnodt * tdelta;
    xmp   += s->dmdt  * tdelta;

    // Euler-Maclaurin numerical integration
    double ft = 0;
    double delta;
    if ((s->is_12h_resonant == true) || (s->is_24h_resonant == true))
    {
      if ((s->atime == 0)
          || (tdelta * s->atime <= 0)
          || (fabs(tdelta) < fabs(s->atime)))
      {
        s->atime = 0;
        s->xni   = s->xnodp;
        s->xli   = s->xlamo;
      }

      if (tdelta > 0)
        delta = stepp;
      else
        delta = stepn;

      bool integrating = true;
      double xndt, xldot, xnddt, xomi, x2omi, x2li;

      while (integrating == true)
      {
        // Synchronous resonance dot terms
        if (s->is_24h_resonant == true)
        {
          xndt  = s->del1 * sin(s->xli - fasx2) + s->del2
                * sin(2 * (s->xli - fasx4))
                + s->del3 * sin(3 * (s->xli - fasx6));
          xldot = s->xni + s->xfact;
          xnddt = s->del1 * cos(s->xli - fasx2) +
              2 * s->del2 * cos(2 * (s->xli - fasx4)) +
              3 * s->del3 * cos(3 * (s->xli - fasx6));
          xnddt = xnddt * xldot;
        }
        // Geopotential resonance terms
        else
        {
          xomi  = s->argument_perigee + s->omgdot * s->atime;
          x2omi = xomi + xomi;
          x2li  = s->xli + s->xli;
          xndt  = s->d2201 * sin(x2omi + s->xli - g22) + s->d2211 * sin(s->xli - g22) +
              s->d3210 * sin(xomi + s->xli - g32)  + s->d3222 * sin(-xomi + s->xli - g32)+
              s->d4410 * sin(x2omi + x2li - g44)+ s->d4422 * sin(x2li - g44) +
              s->d5220 * sin(xomi + s->xli - g52)  + s->d5232 * sin(-xomi + s->xli - g52)+
              s->d5421 * sin(xomi + x2li - g54) + s->d5433 * sin(-xomi + x2li - g54);
          xldot = s->xni + s->xfact;
          xnddt = s->d2201 * cos(x2omi + s->xli - g22) + s->d2211 * cos(s->xli - g22) +
              s->d3210 * cos(xomi + s->xli - g32) + s->d3222 * cos(-xomi + s->xli - g32) +
              s->d5220 * cos(xomi + s->xli - g52) + s->d5232 * cos(-xomi + s->xli - g52) +
              2 * (s->d4410 * cos(x2omi + x2li - g44) +
                  s->d4422 * cos(x2li - g44) + s->d5421 * cos(xomi + x2li - g54) +
                  s->d5433 * cos(-xomi + x2li - g54));
          xnddt = xnddt * xldot;
        }

        // Integrator
        if (fabs(tdelta - s->atime) < stepp)
        {
          ft          = tdelta - s->atime;
          integrating = false;
        }

        if (integrating == true)
        {
          s->xli   = s->xli + xldot * delta + xndt * step2;
          s->xni   = s->xni + xndt * delta + xnddt * step2;
          s->atime = s->atime + delta;
        }
      }

      nm = s->xni + xndt * ft + xnddt * ft * ft * 0.5;
      double xl = s->xli + xldot * ft + xndt * ft * ft * 0.5;

      if (s->is_12h_resonant)
      {
        xmp     = xl - 2 * xnode + 2 * theta;
        s->dndt = nm - s->xnodp;
      }
      else
      {
        xmp     = xl - xnode - omega + theta;
        s->dndt = nm - s->xnodp;
      }
      nm = s->xnodp + s->dndt;

#ifdef MATH_TRACE
      printf("======================================== dp1\n");
      printf("tdelta %+.15e\n", tdelta);
      printf("atime  %+.15e\n", s->atime);
      printf("xndt   %+.15e\n", xndt);
      printf("xldot  %+.15e\n", xldot);
      printf("xnddt  %+.15e\n", xnddt);
      printf("xomi   %+.15e\n", xomi);
      printf("x2omi  %+.15e\n", x2omi);
      printf("x2li   %+.15e\n", x2li);
      printf("xli    %+.15e\n", s->xli);
      printf("xni    %+.15e\n", s->xni );
#endif
    }
  }

  if (nm <= 0)
  {
    return -2;
  }

  double am = pow((XKE / nm), TWOTHIRD) * tempa * tempa;
  nm = XKE / pow(am, 1.5);
  em = em - tempe;

  if ((em >= 1) || (em < -0.001))
  {
    return -3;
  }

  // Avoid division by zero
  if (em < 1.0e-6)
  {
    em  = 1.0e-6;
  }

         xmp += s->xnodp * templ;
  double xlm  = xmp + omega + xnode;
  double em2  = em * em;

  xnode  = fmod(xnode, TAU);
  omega  = fmod(omega, TAU);
  xlm    = fmod(xlm, TAU);
  xmp    = fmod(xlm - omega - xnode, TAU);

  s->inclination_lp      = inclm;
  s->eccentricity_lp     = em;
  s->right_asc_node_lp   = xnode;
  s->argument_perigee_lp = omega;
  s->mean_anomaly_lp     = xmp;
  double sinip = sin(s->inclination_lp);
  double cosip = cos(s->inclination_lp);

#ifdef MATH_TRACE
  printf("---------------------------------------- p2\n");
  printf("am     %+.15e\n", am);
  printf("nm     %+.15e\n", nm);
  printf("em     %+.15e\n", em);
  printf("xlm    %+.15e\n", xlm);
  printf("em2    %+.15e\n", em2);
  printf("omega  %+.15e\n", omega);
  printf("inclp  %+.15e\n", s->inclination_lp);
  printf("ep     %+.15e\n", s->eccentricity_lp);
  printf("nodep  %+.15e\n", s->right_asc_node_lp);
  printf("argpp  %+.15e\n", s->argument_perigee_lp);
  printf("mp     %+.15e\n", s->mean_anomaly_lp);
  printf("sinip  %+.15e\n", sinip);
  printf("cosip  %+.15e\n", cosip);
#endif

  // Add lunar-solar periodics
  if (s->is_deep_space == true)
  {
    // Constants
    const double zes = 0.01675;
    const double zel = 0.05490;

    // Calculate time varying periodics
    double zm    = s->zmos + ZNS * tdelta;
    double zf    = zm + 2 * zes * sin(zm);
    double sinzf = sin(zf);
    double f2    =  0.5 * sinzf * sinzf - 0.25;
    double f3    = -0.5 * sinzf * cos(zf);
    double ses   = s->se2* f2 + s->se3 * f3;
    double sis   = s->si2 * f2 + s->si3 * f3;
    double sls   = s->sl2 * f2 + s->sl3 * f3 + s->sl4 * sinzf;
    double sghs  = s->sgh2 * f2 + s->sgh3 * f3 + s->sgh4 * sinzf;
    double shs   = s->sh2 * f2 + s->sh3 * f3;
    zm    = s->zmol + ZNL * tdelta;

    zf    = zm + 2 * zel * sin(zm);
    sinzf = sin(zf);
    f2    =  0.5 * sinzf * sinzf - 0.25;
    f3    = -0.5 * sinzf * cos(zf);
    double sel   = s->ee2 * f2 + s->e3 * f3;
    double sil   = s->xi2 * f2 + s->xi3 * f3;
    double sll   = s->xl2 * f2 + s->xl3 * f3 + s->xl4 * sinzf;
    double sghl  = s->xgh2 * f2 + s->xgh3 * f3 + s->xgh4 * sinzf;
    double shll  = s->xh2 * f2 + s->xh3 * f3;
    double pe    = ses + sel;
    double pinc  = sis + sil;
    double pl    = sls + sll;
    double pgh   = sghs + sghl;
    double ph    = shs + shll;

    pe    = pe - s->peo;
    pinc  = pinc - s->pinco;
    pl    = pl - s->plo;
    pgh   = pgh - s->pgho;
    ph    = ph - s->pho;
    s->inclination_lp  += pinc;
    s->eccentricity_lp += pe;
    sinip = sin(s->inclination_lp);
    cosip = cos(s->inclination_lp);

    // Apply periodics directly
    if (s->inclination_lp >= 0.2) // Lyddane choice
    {
      ph  = ph / sinip;
      pgh = pgh - cosip * ph;

      s->argument_perigee_lp += pgh;
      s->right_asc_node_lp   += ph;
      s->mean_anomaly_lp     += pl;
    }
    else
    {
      // Apply periodics with Lyddane modifications
      double sinop  = sin(s->right_asc_node_lp);
      double cosop  = cos(s->right_asc_node_lp);
      double alfdp  = sinip * sinop;
      double betdp  = sinip * cosop;
      double dalf   =  ph * cosop + pinc * cosip * sinop;
      double dbet   = -ph * sinop + pinc * cosip * cosop;
      alfdp  = alfdp + dalf;
      betdp  = betdp + dbet;
      s->right_asc_node_lp  = fmod(s->right_asc_node_lp, TAU);

      // Wrap negative node for atan below
      if (s->right_asc_node_lp < 0)
      {
        s->right_asc_node_lp += TAU;
      }

      double xls    = s->mean_anomaly_lp + s->argument_perigee_lp + cosip * s->right_asc_node_lp;
      xls          += pl + pgh - pinc * s->right_asc_node_lp * sinip;
      double xnoh   = s->right_asc_node_lp;
      s->right_asc_node_lp  = atan2(alfdp, betdp);

      // Wrap negative node for fabs below
      if (s->right_asc_node_lp < 0)
      {
        s->right_asc_node_lp += TAU;
      }

      if (fabs(xnoh - s->right_asc_node_lp) > PI)
      {
        if (s->right_asc_node_lp < xnoh)
        {
          s->right_asc_node_lp = s->right_asc_node_lp + TAU;
        }
        else
        {
          s->right_asc_node_lp = s->right_asc_node_lp - TAU;
        }
      }

      s->mean_anomaly_lp    += pl;
      s->argument_perigee_lp = xls - s->mean_anomaly_lp - cosip * s->right_asc_node_lp;
    }

#ifdef MATH_TRACE
    printf("======================================== dp2\n");
    printf("xmp     %+.15e\n", xmp);
    printf("xlm     %+.15e\n", xlm);
    printf("em2     %+.15e\n", em2);
    printf("omega   %+.15e\n", omega);
    printf("incl_lp %+.15e\n", s->inclination_lp);
    printf("node_lp %+.15e\n", s->right_asc_node_lp);
    printf("argp_lp %+.15e\n", s->argument_perigee_lp);
    printf("ecc_lp  %+.15e\n", s->eccentricity_lp);
    printf("mo_lp   %+.15e\n", s->mean_anomaly_lp);
#endif

    if (s->inclination_lp < 0)
    {
      s->inclination_lp      *= -1;
      s->right_asc_node_lp   += PI;
      s->argument_perigee_lp -= PI;
    }

    if ((s->eccentricity_lp < 0)
     || (s->eccentricity_lp > 1))
    {
      return -3;
    }

    // Long period periodics
    sinip    =  sin(s->inclination_lp);
    cosip    =  cos(s->inclination_lp);
    s->aycof = -0.5 * J3DIVJ2 * sinip;

    // Avoid division by zero for s->inclination_lp = 180 deg
    if (fabs(cosip + 1) > 1.5e-12)
    {
      s->xlcof = -0.25 * J3DIVJ2 * sinip * (3 + 5 * cosip) / (1 + cosip);
    }
    else
    {
      s->xlcof = -0.25 * J3DIVJ2 * sinip * (3 + 5 * cosip) / 1.5e-12;
    }
  }

  double axnl    = s->eccentricity_lp * cos(s->argument_perigee_lp);
  double a1e2inv = 1 / (am * (1 - pow(s->eccentricity_lp, 2)));
  double aynl    = s->eccentricity_lp * sin(s->argument_perigee_lp)
                 + a1e2inv * s->aycof;
  double xl      = s->mean_anomaly_lp + s->argument_perigee_lp
                 + s->right_asc_node_lp + a1e2inv * s->xlcof * axnl;

  // Kepler's equation
  double  u       = fmod(xl - s->right_asc_node_lp, TAU);
  double  eo1     = u;
  double  kdelta  = 9999.9;
  uint8_t ktr     = 0;
  double  sineo1, coseo1;

#ifdef MATH_TRACE
  printf("---------------------------------------- p3\n");
  printf("xlcof   %+.15e\n", s->xlcof);
  printf("axnl    %+.15e\n", axnl);
  printf("a1e2inv %+.15e\n", a1e2inv);
  printf("aynl    %+.15e\n", aynl);
  printf("xl      %+.15e\n", xl);
  printf("u       %+.15e\n", u);
#endif

  while ((fabs(kdelta) >= tolerance) && (ktr < maxiter))
  {
    sineo1 = sin(eo1);
    coseo1 = cos(eo1);
    kdelta = (u - aynl * coseo1 + axnl * sineo1 - eo1)
           / (1 - coseo1 * axnl - sineo1 * aynl);

    if(fabs(kdelta) >= 0.95)
    {
      kdelta = kdelta > 0 ? 0.95 : -0.95;
    }

#ifdef MATH_TRACE
    printf("> >>>>>> Kepler\n");
    printf("> kdelta %+.15e\n", kdelta);
    printf("> sineo1 %+.15e\n", sineo1);
    printf("> coseo1 %+.15e\n", coseo1);
    printf("> eo1    %+.15e\n", eo1);
    printf("> ktr    %d\n", ktr);
#endif

    eo1 += kdelta;
    ktr++;
  }

  // Short period preliminary quantities
  double ecose = axnl * coseo1 + aynl * sineo1;
  double esine = axnl * sineo1 - aynl * coseo1;
  double el2   = axnl * axnl   + aynl * aynl;
  double pl    = am * (1 - el2);
  double mrt;

#ifdef MATH_TRACE
  printf("---------------------------------------- p4\n");
  printf("ecose %+.15e\n", ecose);
  printf("esine %+.15e\n", esine);
  printf("el2   %+.15e\n", el2);
  printf("pl    %+.15e\n", pl);
#endif

  if (pl < 0)
  {
    return -4;
  }
  else
  {
    double rl     = am * (1 - ecose);
    double rdotl  = sqrt(am) * esine/rl;
    double rvdotl = sqrt(pl) / rl;
    double betal  = sqrt(1 - el2);
    double sinu   = am / rl * (sineo1 - aynl - axnl * esine / (1 + betal));
    double cosu   = am / rl * (coseo1 - axnl + aynl * esine / (1 + betal));
    double su     = atan2(sinu, cosu);
    double sin2u  = (cosu + cosu) * sinu;
    double cos2u  = 1 - 2 * sinu * sinu;
    double temp1  = 0.5 * J2 * (1 / pl);
    double temp2  = temp1 * (1 / pl);

    // Update for short period periodics
    if (s->is_deep_space == true)
    {
      double cosip2 = cosip * cosip;
      s->con41      = 3 * cosip2 - 1;
      s->x1mth2     = 1 - cosip2;
      s->x7thm1     = 7 * cosip2 - 1;
    }

           mrt   = rl * (1 - 1.5 * temp2 * betal * s->con41)
                 + 0.5 * temp1 * s->x1mth2 * cos2u;
           su    = su - 0.25 * temp2 * s->x7thm1 * sin2u;
           xnode = s->right_asc_node_lp + 1.5 * temp2 * cosip * sin2u;
    double xinc  = s->inclination_lp    + 1.5 * temp2 * cosip * sinip * cos2u;
    double mvt   = rdotl  - nm * temp1 * s->x1mth2  * sin2u / XKE;
    double rvdot = rvdotl + nm * temp1 * (s->x1mth2 * cos2u
                 + 1.5 * s->con41) / XKE;

#ifdef MATH_TRACE
    printf("---------------------------------------- p5\n");
    printf("rl     %+.15e\n", rl);
    printf("rdotl  %+.15e\n", rdotl);
    printf("rvdotl %+.15e\n", rvdotl);
    printf("betal  %+.15e\n", betal);
    printf("sinu   %+.15e\n", sinu);
    printf("cosu   %+.15e\n", cosu);
    printf("su     %+.15e\n", su);
    printf("sin2u  %+.15e\n", sin2u);
    printf("cos2u  %+.15e\n", cos2u);
    printf("temp1  %+.15e\n", temp1);
    printf("con41  %+.15e\n", s->con41);
    printf("x1mth2 %+.15e\n", s->x1mth2);
    printf("temp2  %+.15e\n", temp2);
    printf("mrt    %+.15e\n", mrt);
    printf("xnode  %+.15e\n", xnode);
    printf("xinc   %+.15e\n", xinc);
    printf("mvt    %+.15e\n", mvt);
    printf("rvdot  %+.15e\n", rvdot);
#endif

    if ((p != NULL) && (v != NULL))
    {
      // Orientation vectors
      double sinsu =  sin(su);
      double cossu =  cos(su);
      double snod  =  sin(xnode);
      double cnod  =  cos(xnode);
      double sini  =  sin(xinc);
      double cosi  =  cos(xinc);
      double xmx   = -snod * cosi;
      double xmy   =  cnod * cosi;
      double ux    =  xmx * sinsu + cnod * cossu;
      double uy    =  xmy * sinsu + snod * cossu;
      double uz    =  sini * sinsu;
      double vx    =  xmx * cossu - cnod * sinsu;
      double vy    =  xmy * cossu - snod * sinsu;
      double vz    =  sini * cossu;

      // Position and velocity vectors
      p->x = (mrt * ux) * RE;
      p->y = (mrt * uy) * RE;
      p->z = (mrt * uz) * RE;
      v->x = (mvt * ux + rvdot * vx) * VKMPS;
      v->y = (mvt * uy + rvdot * vy) * VKMPS;
      v->z = (mvt * uz + rvdot * vz) * VKMPS;

#ifdef MATH_TRACE
      printf("---------------------------------------- p6\n");
      printf("t %lf\n", tdelta);
      printf("px %12.6lf\n", p->x);
      printf("py %12.6lf\n", p->y);
      printf("pz %12.6lf\n", p->z);
      printf("vx %12.6lf\n", v->x);
      printf("vy %12.6lf\n", v->y);
      printf("vz %12.6lf\n", v->z);
#endif
    }
  }

  // Satellite decayed?
  if (mrt < 1.0)
  {
    return -5;
  }

  return 0;
}

/*
 * Get observational data about the satellite from ground station
 *
 * Inputs:  s         - sat struct pointer with initialized orbital data
 *          time      - Unix timestamp of observation
 *          time_ms   - Milllisecond portion of the above
 *          obs_geo   - Geodetic coordinates of the ground station
 * Outputs: result    - Observational data
 * Returns: 0         - Success
 *         -1         - Invalid inputs or parametres
 *         -2         - Negative mean motion
 *         -3         - Eccentricity out of range (e >= 1; e < -1.0e-12)
 *         -4         - Short period preliminary quantities error
 *         -5         - Decayed satellite
 */
int
sat_observe
(
        sat*    s,
        time_t  time,
        float   time_ms,
  const vec3*   obs_geo,
        obs*    result
)
{
  if ((s        == NULL) ||
      (time_ms  >= 1000) ||
      (obs_geo  == NULL) ||
      (result   == NULL))
  {
    return -1;
  }

  double tdelta = difftime(time + time_ms / 1000,
                           s->epoch + s->epoch_ms / 1000) / 60;

  vec3 posteme, velteme;

  int retval = sat_propagate(s, tdelta, 10, 1.0e-12, &posteme, &velteme);

  if (retval != 0)
  {
    return retval;
  }
  // Switching to ECEF common frame to fix to Earth
  vec3 posecef, velecef, obsposecef, obsvelecef;

  teme2ecef(&posteme, &velteme, unix2jul(time, time_ms), &posecef, &velecef);

  result->latlonalt = ecef2geo(&posecef);

  // Vis-viva equation
  result->velocity  = sqrt(GM * (2 / (RE + result->latlonalt.alt)
                               - 1 / (RE * s->aodp)));

  obsposecef = geo2ecef(obs_geo);

  // Observer to satellite vector in ECEF frame
  vec3 posdiffecef  = vec3_add(1, &posecef, -1, &obsposecef);

  result->azelrng   = ecef2azelrng(&obsposecef, &posdiffecef);
  result->rng_rate  = vec3_dot(&posdiffecef, &velecef) / result->azelrng.rng;

  // TODO: Implementation pending
  result->is_illum  = false;

  return retval;
}

int
sat_passes
(
        sat*         s,
  const time_t*      start_time,
  const time_t*      stop_time,
  const vec3*        obs_geo,
        unsigned int delta_t,
        double       horizon
)
{
  obs           o            = {0};
  unsigned int  aos_count    = 0;
  unsigned int  los_count    = 0;
  double        prev_el      = -TAU;
  double        prev_prev_el = -TAU;

  unsigned int* AOS_t;
  unsigned int* LOS_t;

  unsigned int  max_samples  = (unsigned int)ceil((*stop_time - *start_time)
                                                  / delta_t);

  AOS_t = malloc(max_samples * sizeof(unsigned int));
  LOS_t = malloc(max_samples * sizeof(unsigned int));

  if ((AOS_t == NULL) && (LOS_t == NULL))
  {
    return -1;
  }

  horizon = fmax(0, horizon);

  // for plotting
  FILE* outfile = fopen("elevations.out", "w");

  // Find principle flares on coarse time pass and count them
  for (time_t t = *start_time; t <= *stop_time; t += delta_t)
  {
    sat_observe(s, t, 0, obs_geo, &o);

    // for plotting
    fprintf(outfile, "%ld,%8.3lf\n", t - *start_time, o.azelrng.el * RAD2DEG);

    if (o.azelrng.el > horizon)
    {
      if (t == *start_time)
      {
        AOS_t[aos_count] = *start_time;
        aos_count++;
      }
      else if (prev_el <= horizon)
      {
        AOS_t[aos_count] = find_aos_los(s, obs_geo, t - delta_t,
                                        delta_t, horizon, true);
        aos_count++;
      }
      else if (t + delta_t > *stop_time)
      {
        LOS_t[los_count] = *stop_time;
        los_count++;
      }
    }
    else if (prev_el > horizon)
    {
      LOS_t[los_count] = find_aos_los(s, obs_geo, t - delta_t,
                                      delta_t, horizon, false);
      los_count++;
    }

    prev_el = o.azelrng.el;
  }

  // for plotting
  fclose(outfile);

  printf("AOS count: %ld\n", aos_count);
  printf("LOS count: %ld\n", los_count);

  for (unsigned int i = 0; i < aos_count; i++)
  {
    printf("AOS[%ld] = %ld\n", i, AOS_t[i] - *start_time);
  }

  for (unsigned int i = 0; i < los_count; i++)
  {
    printf("LOS[%ld] = %ld\n", i, LOS_t[i] - *start_time);
  }

  if (aos_count != los_count)
  {
    return -2;
  }

  double tca_el = 0;

  printf("Flare count: %ld\n", aos_count);

  // TODO: Do more optimal local maximum search
  for (unsigned int i = 0; i < aos_count; i ++)
  {
    for (time_t t = AOS_t[i]; t <= LOS_t[i]; t += (LOS_t[i] - AOS_t[i]) / 13)
    {
      sat_observe(s, t, 0, obs_geo, &o);
      tca_el = fmax(tca_el, o.azelrng.el);
    }
    printf("TCA el coarse = %lf\n", tca_el * RAD2DEG);
    tca_el = 0;
  }

//  For speed comparison
  for (unsigned int i = 0; i < aos_count; i ++)
  {
    for (time_t t = AOS_t[i]; t <= LOS_t[i]; t++)
    {
      sat_observe(s, t, 0, obs_geo, &o);
      tca_el = fmax(tca_el, o.azelrng.el);
    }
    printf("TCA el = %lf\n", tca_el * RAD2DEG);
    tca_el = 0;
  }

  free(AOS_t);
  free(LOS_t);

  return 0;
}

/*
 * Recursively zero in on an AOS or an LOS event down to 1 sec resolution
 *
 * Inputs:  s          - sat struct pointer with initialized orbital data
 *          obs_geo    - Geodetic coordinates of the ground station
 *          start_time - Unix timestamp of observation
 *          delta_t    - Time step
 *          horizon    - Elevation of observer horizon, [0; pi/2] rad
 *          is_aos     - Are we looking for AOS? (false for LOS)
 * Outputs: None
 * Returns: zero crossing unix time with second precision
 */
time_t
find_aos_los
(
        sat*         s,
  const vec3*        obs_geo,
  const time_t       start_time,
        unsigned int delta_t,
        double       horizon,
        bool         is_aos
)
{
  obs    o = {0};
  double diff;

  unsigned int half_dt = ceil(delta_t / 2.0);

  if (delta_t > 1)
  {
    sat_observe(s, start_time + half_dt, 0, obs_geo, &o);

    if (is_aos == true)
    {
      diff = o.azelrng.el - horizon;
    }
    else
    {
      diff = horizon - o.azelrng.el;
    }

    if (diff > 0)
    {
      return find_aos_los(s, obs_geo, start_time,
                          half_dt, horizon, is_aos);
    }
    else
    {
      return find_aos_los(s, obs_geo, start_time + half_dt,
                          half_dt, horizon, is_aos);
    }
  }

  return start_time;
}
