/*
 * libsgp4ansi.c - an ANSI C-11 SGP4/SDP4 implementation library for sgp4ansid.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#include <ctype.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "libsgp4ansi.h"
#include "epoch.h"
#include "const.h"
#include "coord.h"
#include "vector.h"

// ************************************************************************* //
//                                VERSION                                    //
// ************************************************************************* //

const char library_version[] = VERSION;

// ************************************************************************* //
//                            PRIVATE PROTOTYPES                             //
// ************************************************************************* //

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
 *
 * Calls: orbit_init
 */
size_t str_trim(const char* str, size_t len, char* out)
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
 * Inputs:  name_str - String containing satellite name
 *          line1    - TLE line 1
 *          line2    - TLE line 2
 * Outputs: s        - sat struct pointer with full orbital data expanded
 * Returns: 0        - Success
 *         -1        - Corrup TLE
 *
 * Calls: orbit_init
 *
 * TODO: Add checksum check
 */
int
sat_load_tle(char* name_str, char* line1, char* line2, sat* s)
{
  // Pre-formatting the raw TLE input
  for (char j = 10; j <= 15; j++)
    if (line1[j] == ' ')
      line1[j] = '_';

  if (line1[44] != ' ')
    line1[43] = line1[44];

  line1[44] = '.';
  if (line1[7] == ' ')
    line1[7] = 'U';

  if (line1[9] == ' ')
    line1[9] = '.';

  for (char j = 45; j <= 49; j++)
    if (line1[j] == ' ')
      line1[j] = '0';

  if (line1[51] == ' ')
    line1[51] = '0';

  if (line1[53] != ' ')
    line1[52] = line1[53];

  line1[53] = '.';
  line2[25] = '.';

  for (char j = 26; j <= 32; j++)
    if (line2[j] == ' ')
      line2[j] = '0';

  if (line1[62] == ' ')
    line1[62] = '0';

  if (line1[68] == ' ')
    line1[68] = '0';

  str_trim(name_str, 24, s->tle.name);

  struct tm epoch_tm, prop_tm;

  unsigned int cardnum;
  unsigned int norad_number1;
  unsigned int norad_number2;
  char int_designator[9];
  int epochyr;
  double epochdays;
  double mean_motion_dt2;
  double mean_motion_ddt6;
  unsigned int nexp;
  double bstar;
  unsigned int bexp;
  int ephem_type;
  int elset_number;
  double inclination;
  double right_asc_node;
  double argument_perigee;
  double mean_anomaly;
  double mean_motion;

  int retval = sscanf(line1,"%2d %5ld %1c %9s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",
         &cardnum, &norad_number1, &s->tle.sec_class, &s->tle.int_designator, &epochyr,
         &epochdays,&mean_motion_dt2, &mean_motion_ddt6, &nexp, &bstar,
         &bexp, &ephem_type, &elset_number);

  if (retval != 14 || cardnum != 1)
  {
    return -1;
  }

  if (line2[52] == ' ') // check for minus sign
    retval = sscanf(line2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %10lf %6ld \n",
           &cardnum, &norad_number2, &inclination, &right_asc_node,
           &s->tle.eccentricity, &argument_perigee, &mean_anomaly,
           &mean_motion, &s->tle.orbit_number);
  else
    retval = sscanf(line2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld \n",
           &cardnum, &norad_number2, &inclination, &right_asc_node,
           &s->tle.eccentricity, &argument_perigee, &mean_anomaly,
           &mean_motion, &s->tle.orbit_number);

  if (retval != 9 || cardnum != 2 || norad_number1 != norad_number2)
  {
    return -1;
  }

  s->tle.norad_number = norad_number2;

  // Convert to SGP4 units
  if (fractday2unix(epochyr, epochdays, &s->tle.epoch, &s->tle.epoch_ms) != 0)
  {
    return -1;
  }
  s->tle.epoch_jul        = unix2jul(&s->tle.epoch, s->tle.epoch_ms);
  s->tle.mean_motion_dt2  = mean_motion_dt2  / (RPD2RADPM * 1440.0);
  s->tle.mean_motion_ddt6 = mean_motion_ddt6 * pow(10, nexp)
                            / (RPD2RADPM * 1440.0 * 1440);
  s->tle.bstar            = bstar * pow(10, bexp);
  s->tle.inclination      = inclination * DEG2RAD;
  s->tle.right_asc_node   = right_asc_node * DEG2RAD;
  s->tle.argument_perigee = argument_perigee * DEG2RAD;
  s->tle.mean_anomaly     = mean_anomaly * DEG2RAD;
  s->tle.mean_motion      = mean_motion * TWOPI / 1440;

  // Recover original mean motion and semimajor axis
  const double a1     = pow(XKE / s->tle.mean_motion, TWOTHIRD);
  const double cosio  = cos(s->tle.inclination);
  const double theta2 = cosio * cosio;
  const double x3thm1 = 3.0 * theta2 - 1.0;
  const double eosq   = s->tle.eccentricity * s->tle.eccentricity;
  const double betao2 = 1.0 - eosq;
  const double betao  = sqrt(betao2);
  const double temp   = (1.5 * J2 / 2.0) * x3thm1 / (betao * betao2);
  const double del1   = temp / (a1 * a1);
  const double a0     = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + del1
                      * 134.0 / 81.0)));
  const double del0   = temp / (a0 * a0);

  s->comm.xnodp       = s->tle.mean_motion / (1.0 + del0);
  s->comm.aodp        = a0 / (1.0 - del0);

  // Find perigee and period
  s->comm.perigee     = (s->comm.aodp * (1.0 - s->tle.eccentricity) - 1) * RE;
  s->comm.period      = TWOPI / s->comm.xnodp;

  // Expand constants
  return sat_init(s);
}

/*
 * Expand SGP4/SDP4 orbit elements from an orbit containing NORAD TLE portion
 *
 * Inputs:  s  - sat struct pointer to an unexpanded elements set
 * Outputs: s  - sat struct pointer with full orbital data
 * Returns: 0  - Success
 *         -1  - Eccenricity out of range [0, 1)
 *         -2  - Inclination out of range [0, pi)
 *
 * TODO: add return values
 */
int
sat_init(sat* s)
{

  /*
   * error checks
   */
  if (s->tle.eccentricity < 0.0 || s->tle.eccentricity > 0.999999)
  {
    return -1;
  }

  if (s->tle.inclination < 0.0 || s->tle.inclination > PI)
  {
    return -2;
  }

  s->comm.cosio       = cos(s->tle.inclination);
  s->comm.sinio       = sin(s->tle.inclination);
  const double theta2 = s->comm.cosio * s->comm.cosio;
  s->comm.x3thm1      = 3.0 * theta2 - 1.0;
  const double eosq   = s->tle.eccentricity * s->tle.eccentricity;
  const double betao2 = 1.0 - eosq;
  const double betao  = sqrt(betao2);

  if (s->comm.period >= 225.0)
  {
    s->comm.is_deep_space = true;
  }
  else
  {
    s->comm.is_deep_space = false;
    s->comm.use_simple_model = false;

    // Use simple model for perigee less than 220 kilometers
    if (s->comm.perigee < 220.0)
    {
      s->comm.use_simple_model = true;
    }
  }

  // find new s4 and qoms2t for perigees below 156km
  double s4 = 1.0 + 78.0 / RE;
  double qoms24 = pow(((120.0 - 78.0) / RE), 4.0);
  if (s->comm.perigee < 156.0)
  {
    s4 = s->comm.perigee - 78.0;
    if (s->comm.perigee < 98.0)
    {
      s4 = 20.0;
    }
    qoms24 = pow((120.0 - s4) * 1.0 / RE, 4.0);
    s4 = s4 / RE + 1;
  }

   // Expand constants
   const double pinvsq = 1.0 / (s->comm.aodp * s->comm.aodp * betao2 * betao2);
   const double tsi    = 1.0 / (s->comm.aodp - s4);
   s->comm.eta         = s->comm.aodp * s->tle.eccentricity * tsi;
   const double etasq  = s->comm.eta * s->comm.eta;
   const double eeta   = s->tle.eccentricity * s->comm.eta;
   const double psisq  = fabs(1.0 - etasq);
   const double coef   = qoms24 * pow(tsi, 4.0);
   const double coef1  = coef / pow(psisq, 3.5);
   const double c2     = coef1 * s->comm.xnodp * (s->comm.aodp * (1.0 + 1.5
                       * etasq + eeta * (4.0 + etasq)) + 0.75 * J2DIV2 * tsi
                       / psisq * s->comm.x3thm1 * (8.0 + 3.0 * etasq
                       * (8.0 + etasq)));
   s->comm.c1          = s->tle.bstar * c2;
   s->comm.a3ovk2      = -J3 / J2DIV2;
   s->comm.x1mth2      = 1.0 - theta2;
   s->comm.c4          = 2.0 * s->comm.xnodp * coef1 * s->comm.aodp * betao2
                       * (s->comm.eta * (2.0 + 0.5 * etasq)
                       + s->tle.eccentricity * (0.5 + 2.0 * etasq) - J2 * tsi
                       / (s->comm.aodp * psisq) * (-3.0 * s->comm.x3thm1 * (1.0
                       - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) + 0.75
                       * s->comm.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq))
                       * cos(2.0 * s->tle.argument_perigee)));
   const double theta4 = theta2 * theta2;
   const double temp1  = 3.0 * J2DIV2 * pinvsq * s->comm.xnodp;
   const double temp2  = temp1 * J2DIV2 * pinvsq;
   const double temp3  = 1.25 * J2DIV2 * pinvsq * pinvsq * s->comm.xnodp;
   s->comm.xmdot       = s->comm.xnodp + 0.5 * temp1 * betao * s->comm.x3thm1
                       + 0.0625 * temp2 * betao * (13.0 - 78.0 * theta2 + 137.0
                       * theta4);
   const double x1m5th = 1.0 - 5.0 * theta2;
   s->comm.omgdot      = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7.0 - 114.0
                       * theta2 + 395.0 * theta4) + temp3 * (3.0 - 36.0
                       * theta2 + 49.0 * theta4);
   const double xhdot1 = -temp1 * s->comm.cosio;
   s->comm.xnodot      = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0
                       * temp3 * (3.0 - 7.0 * theta2)) * s->comm.cosio;
   s->comm.xnodcf      = 3.5 * betao2 * xhdot1 * s->comm.c1;
   s->comm.t2cof       = 1.5 * s->comm.c1;

   if (fabs(s->comm.cosio + 1.0) > 1.5e-12)
   {
     s->comm.xlcof = 0.125 * s->comm.a3ovk2 * s->comm.sinio * (3.0 + 5.0 * s->comm.cosio) / (1.0 + s->comm.cosio);
   }
   else
   {
     s->comm.xlcof = 0.125 * s->comm.a3ovk2 * s->comm.sinio * (3.0 + 5.0 * s->comm.cosio) / 1.5e-12;
   }

   s->comm.aycof = 0.25 * s->comm.a3ovk2 * s->comm.sinio;
   s->comm.x7thm1 = 7.0 * theta2 - 1.0;

   if (s->comm.is_deep_space)
   {
     s->deep.gsto = jul2gst(s->tle.epoch_jul);

     //DeepSpaceInitialise(eosq, s->comm.sinio, s->comm.cosio, betao,
     //                    theta2, betao2,
     //                    s->comm.xmdot, s->comm.omgdot, s->comm.xnodot);
   }
   else
   {
     // Shortcut for round orbits
     double c3 = 0.0;
     if (s->tle.eccentricity > 1.0e-4)
     {
       c3 = coef * tsi * s->comm.a3ovk2 * s->comm.xnodp * s->comm.sinio
          / s->tle.eccentricity;
     }

     s->near.c5      = 2.0 * coef1 * s->comm.aodp * betao2 * (1.0 + 2.75
                     * (etasq + eeta) + eeta * etasq);
     s->near.omgcof  = s->tle.bstar * c3 * cos(s->tle.argument_perigee);

     // Shortcut for round orbits
     s->near.xmcof   = 0.0;
     if (s->tle.eccentricity > 1.0e-4)
     {
       s->near.xmcof = -TWOTHIRD * coef * s->tle.bstar / eeta;
     }

     s->near.delmo = pow(1.0 + s->comm.eta * (cos(s->tle.mean_anomaly)), 3.0);
     s->near.sinmo = sin(s->tle.mean_anomaly);

     if (!s->comm.use_simple_model)
     {
       const double c1sq = s->comm.c1 * s->comm.c1;
       s->near.d2        = 4.0 * s->comm.aodp * tsi * c1sq;
       const double temp = s->near.d2 * tsi * s->comm.c1 / 3.0;
       s->near.d3        = (17.0 * s->comm.aodp + s4) * temp;
       s->near.d4        = 0.5 * temp * s->comm.aodp * tsi * (221.0
                         * s->comm.aodp + 31.0 * s4) * s->comm.c1;
       s->near.t3cof     = s->near.d2 + 2.0 * c1sq;
       s->near.t4cof     = 0.25 * (3.0 * s->near.d3 + s->comm.c1
                         * (12.0 * s->near.d2 + 10.0 * c1sq));
       s->near.t5cof     = 0.2 * (3.0 * s->near.d4 + 12.0 * s->comm.c1 *
                         s->near.d3 + 6.0 * s->near.d2 * s->near.d2 + 15.0 *
                         c1sq * (2.0 * s->near.d2 + c1sq));
     }
   }
  return 0;
}
