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
#include <time.h>

#include "libsgp4ansi.h"
#include "epoch.h"
#include "const.h"
#include "coord.h"
#include "vector.h"

// ************************************************************************* //
//                                VERSION                                    //
// ************************************************************************* //

const char libsgp4ansi_version[] = VERSION;

// ************************************************************************* //
//                            PRIVATE PROTOTYPES                             //
// ************************************************************************* //

// Trim leading and trailing whitespaces from a string
size_t
str_trim(const char*, size_t, char*);

/*
AIAA call structure:
twoline2rv()
  - reading TLE
  - converting to SGP4 units
  sgp4init()
    - init with params
    initl()
    if deepspace
      dscom()
      dpper()
      dsinit()
    sgp4(0)

sgp4()
  if deepspace
    dspace()
    dpper()          <------
*/


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

void sat_print(sat* s, char* caption)
{
  FILE *f = fopen("ansi_sat_print.txt", "a");
  fprintf(f, "======================= %s ====================\n", caption);
  fprintf(f, "------------------------ NORAD TLE --------------------------\n");
  fprintf(f, "name             %s\n", s->name);
  fprintf(f, "sec_class        %c\n", s->sec_class);
  fprintf(f, "int_designator   %s\n", s->int_designator);
  fprintf(f, "epoch            %s%lfms\n", ctime(&s->epoch), s->epoch_ms);
  fprintf(f, "julian_epoch     %22.15lf\n", s->julian_epoch);
  fprintf(f, "mean_motion_dt2  %+.15e\n", s->mean_motion_dt2);
  fprintf(f, "mean_motion_ddt6 %+.15e\n", s->mean_motion_ddt6);
  fprintf(f, "Bstar            %+.15e\n", s->Bstar);
  fprintf(f, "inclination      %+.15e\n", s->inclination);
  fprintf(f, "right_asc_node   %+.15e\n", s->right_asc_node);
  fprintf(f, "eccentricity     %+.15e\n", s->eccentricity);
  fprintf(f, "argument_perigee %+.15e\n", s->argument_perigee);
  fprintf(f, "mean_anomaly     %+.15e\n", s->mean_anomaly);
  fprintf(f, "mean_motion      %+.15e\n", s->mean_motion);
  fprintf(f, "norad_number     %d\n", s->norad_number);
  fprintf(f, "orbit_number     %d\n", s->orbit_number);
  fprintf(f, "-------------------------- Flags ----------------------------\n");
  fprintf(f, "is_deep_space    %d\n", s->is_deep_space);
  fprintf(f, "use_simple_model %d\n", s->use_simple_model);
  fprintf(f, "is_24h_resonant  %d\n", s->is_24h_resonant);
  fprintf(f, "is_12h_resonant  %d\n", s->is_12h_resonant);
  fprintf(f, "---------------- Standard orbital elements ------------------\n");
  fprintf(f, "GSTo             %+.15e\n", s->GSTo);
  fprintf(f, "xnodp            %+.15e\n", s->xnodp);
  fprintf(f, "aodp             %+.15e\n", s->aodp);
  fprintf(f, "perigee          %+.15e\n", s->perigee);
  fprintf(f, "perigee_alt      %+.15e\n", s->perigee_alt);
  fprintf(f, "period           %+.15e\n", s->period);
  fprintf(f, "---------------------- Common constants ---------------------\n");
  fprintf(f, "aycof            %+.15e\n", s->aycof);
  fprintf(f, "C1               %+.15e\n", s->C1);
  fprintf(f, "C4               %+.15e\n", s->C4);
  fprintf(f, "eta              %+.15e\n", s->eta);
  fprintf(f, "omgdot           %+.15e\n", s->omgdot);
  fprintf(f, "t2cof            %+.15e\n", s->t2cof);
  fprintf(f, "x1mth2           %+.15e\n", s->x1mth2);
  fprintf(f, "x1m5th2          %+.15e\n", s->x1m5th2);
  fprintf(f, "x1m7th2          %+.15e\n", s->con41);
  fprintf(f, "x7thm1           %+.15e\n", s->x7thm1);
  fprintf(f, "xlcof            %+.15e\n", s->xlcof);
  fprintf(f, "xnodcf           %+.15e\n", s->xnodcf);
  fprintf(f, "xnodot           %+.15e\n", s->xnodot);
  fprintf(f, "xmdot            %+.15e\n", s->xmdot);
  fclose(f);
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
sat_load_tle(char* tlestr0, char* tlestr1, char* tlestr2, sat* s)
{
  // Pre-formatting the raw TLE input
  for (char j = 10; j <= 15; j++)
    if (tlestr1[j] == ' ')
      tlestr1[j] = '_';

  if (tlestr1[44] != ' ')
    tlestr1[43] = tlestr1[44];

  tlestr1[44] = '.';
  if (tlestr1[7] == ' ')
    tlestr1[7] = 'U';

  if (tlestr1[9] == ' ')
    tlestr1[9] = '.';

  for (char j = 45; j <= 49; j++)
    if (tlestr1[j] == ' ')
      tlestr1[j] = '0';

  if (tlestr1[51] == ' ')
    tlestr1[51] = '0';

  if (tlestr1[53] != ' ')
    tlestr1[52] = tlestr1[53];

  tlestr1[53] = '.';
  tlestr2[25] = '.';

  for (char j = 26; j <= 32; j++)
    if (tlestr2[j] == ' ')
      tlestr2[j] = '0';

  if (tlestr1[62] == ' ')
    tlestr1[62] = '0';

  if (tlestr1[68] == ' ')
    tlestr1[68] = '0';

  str_trim(tlestr0, 24, s->name);

  struct tm epoch_tm, prop_tm;

  // TODO: Do real checks

  int cardnum, epochyr, nexp, Bexp, checksum, ephem_type, elset_number;
  double nddot, Bstar, epochdays;

  int retval = sscanf(tlestr1,"%2d %5ld %1c %8s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",
         &cardnum, &s->norad_number, &s->sec_class, s->int_designator, &epochyr,
         &epochdays,&s->mean_motion_dt2, &nddot, &nexp, &Bstar,
         &Bexp, &ephem_type, &elset_number);

  if (retval != 13)
  {
    return -1;
  }

  if (tlestr2[52] == ' ') // check for minus sign
    retval = sscanf(tlestr2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %10lf %6ld \n",
           &cardnum,&s->norad_number, &s->inclination, &s->right_asc_node, &s->eccentricity, &s->argument_perigee,
           &s->mean_anomaly, &s->mean_motion, &s->orbit_number);
  else
    retval = sscanf(tlestr2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld \n",
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

  s->julian_epoch = unix2jul(&s->epoch, s->epoch_ms);

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
 * Expand SGP4/SDP4 orbit elements from a sat record containing NORAD TLE data
 *
 * Inputs:  s - sat struct pointer to a record containing NORAD TLE data
 * Outputs: s - sat struct pointer with full orbital data
 * Returns: 0 - Success
 *         -1 - Eccentricity out of range [0, 1)
 *         -2 - Inclination out of range [0, pi]
 */
int
sat_init(sat* s)
{

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
  //s->aodp        = a0 / (1 - delta0);
  s->aodp        = pow(XKE / s->xnodp, TWOTHIRD);
  s->perigee     = (s->aodp * (1 - s->eccentricity)) * RE; // TODO: Remove?
  s->perigee_alt = s->perigee - RE;
  s->period      = TWOPI / s->xnodp;

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
      s->t2cof   = 1.5 * s->C1; // TODO: Remove?
  // Division by zero check then inclination = 180 deg
      s->xlcof   = 0.125 * A3OVK2 * sinio * (3.0 + 5.0 * cosio)
                 / ((fabs(cosio+1.0) > 1.5e-12) ? ((1.0 + cosio)) : (1.5e-12));
      s->aycof   = 0.25 * A3OVK2 * sinio;

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

  if (s->is_deep_space == true) // Deep space init here
  {

    s->GSTo  = jul2gst(s->julian_epoch);

    // Constants
    const double zes    =  0.01675; // TODO: Move to macros?
    const double zel    =  0.05490;
    const double c1ss   =  2.9864797e-6;
    const double c1l    =  4.7968065e-7;
    const double zsinis =  0.39785416;
    const double zcosis =  0.91744867;
    const double zcosgs =  0.1945905;
    const double zsings = -0.98088458;

    double snodm  = sin(s->right_asc_node); // TODO: duplicates?
    double cnodm  = cos(s->right_asc_node);
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

    double day    = s->julian_epoch + 18261.5 - DEC31_1949_0000H;
    double xnodce = fmod(4.5236020 - 9.2422029e-4 * day, TWOPI);
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

    printf("======================================== id1\n");
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

    s->zmos = fmod(6.2565837 + 0.017201977 * day, TWOPI);

    // Do solar terms
    double cosq  = cos(s->right_asc_node);
    double sinq  = sin(s->right_asc_node);
    double zcosg = zcosgs;
    double zsing = zsings;
    double zcosi = zcosis;
    double zsini = zsinis;
    double zcosh = cosq; // TODO: Move down?
    double zsinh = sinq; // TODO: Move down?
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
        zcosh = zcoshl * cnodm + zsinhl * snodm;
        zsinh = snodm * zcoshl - cnodm * zsinhl;
        cc    = c1l;
      }
    }

    s->zmol = fmod(4.7199672 + 0.22997150  * day - gam, TWOPI);

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

    //dpper(s, 0); // TODO: Investigate further

    const double q22    = 1.7891679e-6;  // TODO: Move to macros?
    const double q31    = 2.1460748e-6;
    const double q33    = 2.2123015e-7;
    const double root22 = 1.7891679e-6;
    const double root44 = 7.3636953e-9;
    const double root54 = 2.1765803e-9;
    const double rptim  = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
    const double root32 = 3.7393792e-7;
    const double root52 = 1.1428639e-7;
    const double znl    = 1.5835218e-4;
    const double zns    = 1.19459e-5;

    // Deep space resonance initialization
    s->is_12h_resonant = false;
    s->is_24h_resonant = false;
    // TODO: Check for period instead for readability?
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
    double ses  =  ss1 * zns * ss5;
    double sis  =  ss2 * zns * (sz11 + sz13);
    double sls  = -zns * ss3 * (sz1 + sz3 - 14 - 6 * eo2);
    double sghs =  ss4 * zns * (sz31 + sz33 - 6);
    double shs  = -zns * ss2 * (sz21 + sz23);

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

//    printf("ses:    %+.15e\n", ses);
//    printf("sis:    %+.15e\n", sis);
//    printf("sls:    %+.15e\n", sls);
//    printf("sghs:   %+.15e\n", sghs);
//    printf("shs:    %+.15e\n", shs);
//    printf("sgs:    %+.15e\n", sgs);

    // Do lunar terms
    s->dedt = ses + s1 * znl * s5;
    s->didt = sis + s2 * znl * (z11 + z13);
    s->dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * eo2);

    double sghl = s4 * znl * (z31 + z33 - 6.0);
    double shll = -znl * s2 * (z21 + z23);

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

    // Calculate deep space resonance effects

    /* TODO: t is zero at init, right?
     *em     = *em + *dedt * t;
     *inclm  = *inclm + *didt * t;
     *argpm  = *argpm + *domdt * t;
     *nodem  = *nodem + *dnodt * t;
     *mm     = *mm + *dmdt * t;
     */

    double aonv  = pow(s->xnodp / XKE, TWOTHIRD);
    double ainv2 = pow(aonv, 2);

    printf("======================================== id6\n");
    printf("[DS] 12h res %d\n", (int)s->is_12h_resonant);
    printf("[DS] 24h res %d\n", (int)s->is_24h_resonant);

    // Initialize the resonance terms
    if ((s->is_12h_resonant == true)
        || (s->is_24h_resonant == true))
    {
      // Geopotential resonance for 12 hour orbits
      if (s->is_12h_resonant == true)
      {
        //cosisq = cosim * cosim;
        //emo    = *em;
        //*em     = ecco;
        //emsqo  = emsq;
        //emsq   = eccsq;
        double eocu   = pow(s->eccentricity, 3); // TODO: Unroll?
        double g201   = -0.306 - (s->eccentricity - 0.64) * 0.440;

//        printf("----------------------------\n");
//        printf("em:     %70.60lf\n", s->eccentricity);
//        printf("emsq:   %70.60lf\n", eosq);
//        printf("eoc:    %70.60lf\n", eocu);

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

//        printf("----------------------------\n");
//        printf("g201:   %+.15e\n", g201);
//        printf("g211:   %+.15e\n", g211);
//        printf("g310:   %+.15e\n", g310);
//        printf("g322:   %+.15e\n", g322);
//        printf("g410:   %+.15e\n", g410);
//        printf("g422:   %+.15e\n", g422);
//        printf("g520:   %+.15e\n", g520);
//        printf("g533:   %+.15e\n", g533);
//        printf("g521:   %+.15e\n", g521);
//        printf("g532:   %+.15e\n", g532);

        double sini2  =  sinim * sinim; // TODO: Copy of a copy of a copy?
        double cosisq =  cosim * cosim; // TODO: Copy of a copy of a copy?

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

//        printf("----------------------------\n");
//        printf("f220:   %+.15e\n", f220);
//        printf("f221:   %+.15e\n", f221);
//        printf("f321:   %+.15e\n", f321);
//        printf("f322:   %+.15e\n", f322);
//        printf("f441:   %+.15e\n", f441);
//        printf("f442:   %+.15e\n", f442);
//        printf("f522:   %+.15e\n", f522);
//        printf("f523:   %+.15e\n", f523);
//        printf("f542:   %+.15e\n", f542);
//        printf("f543:   %+.15e\n", f543);

               temp1 = 3.0 * pow(s->xnodp, 2) * ainv2;
        double temp  = temp1 * root22;
            s->d2201 = temp * f220 * g201;
            s->d2211 = temp * f221 * g211;
               temp1 = temp1 * aonv; // TODO: Rename/Remove?
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
                     - 2 * s->GSTo, TWOPI);
            s->xfact = s->xmdot + s->dmdt
                     + 2 * (s->xnodot + s->dnodt - rptim) - s->xnodp;
        //double em    = emo;
        //double emsq  = emsqo;

//        printf("----------------------------\n");
//        printf("ainv2:  %+.15e\n", ainv2);
//        printf("d2201:  %+.15e\n", d2201);
//        printf("d2211:  %+.15e\n", d2211);
//        printf("d3210:  %+.15e\n", d3210);
//        printf("d3222:  %+.15e\n", d3222);
//        printf("d4410:  %+.15e\n", d4410);
//        printf("d4422:  %+.15e\n", d4422);
//        printf("d5220:  %+.15e\n", d5220);
//        printf("d5232:  %+.15e\n", d5232);
//        printf("d5232:  %+.15e\n", d5232);
//        printf("d5433:  %+.15e\n", d5433);
      }

      // Synchronous resonance terms

      double g200, g310, g300, f220, f311, f330;

      if (s->is_24h_resonant == true)
      {
            g200 = 1 + eo2 * (-2.5 + 0.8125 * eo2);
            g310 = 1 + 2 * eo2;
            g300 = 1 + eo2 * (-6 + 6.60937 * eo2);
            f220 = 0.75 * (1 + cosim) * (1 + cosim);
            f311 = 0.9375 * sinim * sinim * (1 + 3 * cosim) - 0.75 * (1 + cosim);
            f330 = 1 + cosim;
            f330 = 1.875 * f330 * f330 * f330;
         s->del1 = 3 * pow(s->xnodp, 2) * ainv2;
         s->del2 = 2 * s->del1 * f220 * g200 * q22;
         s->del3 = 3 * s->del1 * f330 * g300 * q33 * aonv;
         s->del1 = s->del1 * f311 * g310 * q31 * aonv;
        s->xfact = s->xmdot + (s->omgdot + s->xnodot) - rptim + s->dmdt
                 + s->domdt + s->dnodt - s->xnodp;
        s->xlamo = fmod(s->mean_anomaly + s->right_asc_node
                 + s->argument_perigee - s->GSTo, TWOPI);

//        printf("----------------------------\n");
//        printf("g200:   %+.15e\n", g200);
//        printf("g310:   %+.15e\n", g310);
//        printf("g300:   %+.15e\n", g300);
//        printf("f220:   %+.15e\n", f220);
//        printf("f311:   %+.15e\n", f311);
//        printf("f330:   %+.15e\n", f330);
//        printf("del1:   %+.15e\n", del1);
//        printf("del2:   %+.15e\n", del2);
//        printf("del3:   %+.15e\n", del3);
      }

      // Initialize the integrator
      s->xli    = s->xlamo;
      s->xni    = s->xnodp;
      s->atime  = 0.0;

//      printf("xfact:  %+.15e\n", s->xfact);
//      printf("xlamo:  %+.15e\n", s->xlamo);
//      printf("xli:    %+.15e\n", s->xli);
//      printf("xni:    %+.15e\n", s->xni);
//      printf("atime:  %+.15e\n", s->atime);
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
 // printf("xmcof:     %20.15lf\n", s->xmcof);
      s->C5      = 2 * coef1 * s->aodp * betao2
                 * (1 + 2.75 * (eta2 + eeta) + eeta * eta2);
      s->omgcof  = s->Bstar * C3 * cos(s->argument_perigee);
     // printf("omgcof:     %20.15lf\n", s->omgcof);
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

  vec3 p, v; // TODO: Test only
  return sat_propagate(s, 1440.0, 4, 1.0e-12, &p, &v);

  // Propagate at zero time since epoch
  //return sat_propagate(s, 0.0, 4, 1.0e-12, NULL, NULL);
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
 *          3         - Long periodics result error
 *          4         - Short period preliminary quantities error
 *          5         - Decaying satellite
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

  double nm    = s->xnodp; // TODO: Rename? Optimize?
  double em    = s->eccentricity; // TODO: Optimize?
  double inclm = s->inclination; // TODO: Optimize?

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

  if (s->is_deep_space == true)
  {
    //tc = tdelta;
    // Deep space contributions to mean elements for perturbing third body
    const double fasx2 = 0.13130908;
    const double fasx4 = 2.8843198;
    const double fasx6 = 0.37448087;
    const double g22   = 5.7686396;
    const double g32   = 0.95240898;
    const double g44   = 1.8014998;
    const double g52   = 1.0508330;
    const double g54   = 4.4108898;
    const double rptim = 4.37526908801129966e-3; // TODO: Move to macros?
    const double stepp =    720.0;
    const double stepn =   -720.0;
    const double step2 = 259200.0;

    // Calculate deep space resonance effects
    s->dndt       = 0;
    double theta  = fmod(s->GSTo + tdelta * rptim, TWOPI); // TODO: Move to struct?

    // Perturbed quantities
    em += s->dedt * tdelta;
    inclm += s->didt  * tdelta;
    omega += s->domdt * tdelta;
    xnode += s->dnodt * tdelta;
    xmp   += s->dmdt  * tdelta;

    // Euler-Maclaurin numerical integration
    double ft = 0; // TODO:Remove?
    double delt; // TODO: Rename
    if ((s->is_12h_resonant == true) || (s->is_24h_resonant == true))
    {
      if ((s->atime == 0)
          || (tdelta * s->atime <= 0.0)
          || (fabs(tdelta) < fabs(s->atime)))
      {
        s->atime = 0.0;
        s->xni   = s->xnodp;
        s->xli   = s->xlamo;
      }

      if (tdelta > 0.0)
        delt = stepp;
      else
        delt = stepn;

      int iretn = 381; // added for do loop TODO: Alternatives?
      int iret  =   0; // added for loop
      double xndt, xldot, xnddt, xomi, x2omi, x2li;

      while (iretn == 381)
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
        else
        {
          // Geopotential resonance terms
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
              2.0 * (s->d4410 * cos(x2omi + x2li - g44) +
                  s->d4422 * cos(x2li - g44) + s->d5421 * cos(xomi + x2li - g54) +
                  s->d5433 * cos(-xomi + x2li - g54));
          xnddt = xnddt * xldot;
        }

        // Integrator
        if (fabs(tdelta - s->atime) >= stepp)
        {
          iret  = 0;
          iretn = 381;
        }
        else // exit here
        {
          ft    = tdelta - s->atime;
          iretn = 0;
        }

        if (iretn == 381)
        {
          s->xli   = s->xli + xldot * delt + xndt * step2;
          s->xni   = s->xni + xndt * delt + xnddt * step2;
          s->atime = s->atime + delt;
        }
//        printf("tdelta: %+.15e\n", tdelta);
//        printf("atime:  %+.15e\n", s->atime);
//        printf("xndt:   %+.15e\n", xndt);
//        printf("xldot:  %+.15e\n", xldot);
//        printf("xnddt:  %+.15e\n", xnddt);
//        printf("xomi:   %+.15e\n", xomi);
//        printf("x2omi:  %+.15e\n", x2omi);
//        printf("x2li:   %+.15e\n", x2li);
//        printf("xli:    %+.15e\n", s->xli);
//        printf("xni:    %+.15e\n", s->xni );
//        printf("---------------------------------\n");
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
        xmp   = xl - xnode - omega + theta;
        s->dndt = nm - s->xnodp;
      }
      nm = s->xnodp + s->dndt;
    }

  }

  if (nm <= 0)
  {
    return -2;
  }

  double am = pow((XKE / nm), TWOTHIRD) * pow(tempa, 2); // TODO: Unroll
  nm = XKE / pow(am, 1.5);
  em = em - tempe;

  if ((em >= 1) || (em < -1.0e-12))
  {
    return -3;
  }

  // Avoid division by zero
  if (em < 1.0e-12) // TODO: Move tolerance to a constant?
  {
    em  = 1.0e-12;
  }

         xmp += s->xnodp * templ;
  double xlm  = xmp + omega + xnode;
  double em2  = pow(em, 2); // TODO: Unroll?
  //double temp = 1 - em2; // TODO: Remove?

  xnode  = fmod(xnode, TWOPI);
  omega  = fmod(omega, TWOPI);
  xlm    = fmod(xlm, TWOPI);
  xmp    = fmod(xlm - omega - xnode, TWOPI);

  s->inclination_lp      = inclm;
  s->eccentricity_lp     = em;
  s->right_asc_node_lp   = xnode;
  s->argument_perigee_lp = omega;
  s->mean_anomaly_lp     = xmp;
  double sinip = sin(s->inclination_lp);
  double cosip = cos(s->inclination_lp);

  printf("---------------------------------------- p2\n");
  printf("am     %+.15e\n", am);
  printf("nm     %+.15e\n", nm);
  printf("em     %+.15e\n", em);
  printf("xlm    %+.15e\n", xlm);
  printf("em2    %+.15e\n", em2);
  printf("em     %+.15e\n", em);
  printf("omega  %+.15e\n", omega);
  printf("inclp  %+.15e\n", s->inclination_lp);
  printf("ep     %+.15e\n", s->eccentricity_lp);
  printf("nodep  %+.15e\n", s->right_asc_node_lp);
  printf("argpp  %+.15e\n", s->argument_perigee_lp);
  printf("mp     %+.15e\n", s->mean_anomaly_lp);
  printf("sinip  %+.15e\n", sinip);
  printf("cosip  %+.15e\n", cosip);

  // Add lunar-solar periodics
  if (s->is_deep_space == true)
  {
    dpper(s, tdelta);

    printf("======================================== dp2\n");
    printf("xmp     %+.15e\n", xmp);
    printf("xlm     %+.15e\n", xlm);
    printf("em2     %+.15e\n", em2);
    printf("omega   %+.15e\n", omega);
    printf("incl_lp %+.15e\n", s->inclination_lp);
    printf("node_lp %+.15e\n", s->right_asc_node_lp);
    printf("argplp  %+.15e\n", s->argument_perigee_lp);
    printf("ecc_lp  %+.15e\n", s->eccentricity_lp);
    printf("mo_lp   %+.15e\n", s->mean_anomaly_lp);

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
  }

/*
  // Long period periodics
  if (s->method == 'd')
  {
    sinip =  sin(xincp);
    cosip =  cos(xincp);
    s->aycof = -0.5*j3oj2*sinip;
    // sgp4fix for divide by zero for xincp = 180 deg
    if ((fabs(cosip+1.0) > 1.5e-12) && (s->operationmode != 'a'))
      s->xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
    else
      s->xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
  }
*/

  double axnl    = s->eccentricity_lp * cos(s->argument_perigee_lp);
  double a1e2inv = 1 / (am * (1 - pow(s->eccentricity_lp, 2)));
  double aynl    = s->eccentricity_lp * sin(s->argument_perigee_lp)
                 + a1e2inv * s->aycof;
  double xl      = s->mean_anomaly_lp + s->argument_perigee_lp
                 + s->right_asc_node_lp + a1e2inv * s->xlcof * axnl;

  // Kepler's equation
  double  u       = fmod(xl - s->right_asc_node_lp, TWOPI);
  double  eo1     = u;
  double  kdelta  = 9999.9;
  uint8_t ktr     = 0;
  double  sineo1, coseo1;

  printf("---------------------------------------- p3\n");
  printf("xlcof   %+.15e\n", s->xlcof);
  printf("axnl    %+.15e\n", axnl);
  printf("a1e2inv %+.15e\n", a1e2inv);
  printf("aynl    %+.15e\n", aynl);
  printf("xl      %+.15e\n", xl);
  printf("u       %+.15e\n", u);

  while ((fabs(kdelta) >= tolerance) && (ktr < maxiter) )
  {
    sineo1 = sin(eo1);
    coseo1 = cos(eo1);
    kdelta = 1 - coseo1 * axnl - sineo1 * aynl;
    kdelta = (u - aynl * coseo1 + axnl * sineo1 - eo1) / kdelta;

    if(fabs(kdelta) >= 0.95)
    {
      kdelta = kdelta > 0 ? 0.95 : -0.95;
    }

    eo1 += kdelta;
    ktr++;

    printf("> >>>>>> Kepler\n");
    printf("> sineo1 %+.15e\n", sineo1);
    printf("> coseo1 %+.15e\n", coseo1);
    printf("> eo1    %+.15e\n", eo1);
    printf("> ktr    %d\n", ktr);
  }

  // Short period preliminary quantities
  double ecose = axnl * coseo1 + aynl * sineo1;
  double esine = axnl * sineo1 - aynl * coseo1;
  double el2   = axnl * axnl   + aynl * aynl;
  double pl    = am * (1 - el2);
  double mrt;

  printf("---------------------------------------- p4\n");
  printf("ecose %+.15e\n", ecose);
  printf("esine %+.15e\n", esine);
  printf("el2   %+.15e\n", el2);
  printf("pl    %+.15e\n", pl);

  if (pl < 0)
  {
    return 4;
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
/*
    // Update for short period periodics
    if (s->method == 'd')
    {
      cosisq                 = cosip * cosip;
      s->con41  = 3.0*cosisq - 1.0;
      s->x1mth2 = 1.0 - cosisq;
      s->x7thm1 = 7.0*cosisq - 1.0;
    }
*/
           mrt   = rl * (1 - 1.5 * temp2 * betal * s->con41)
                 + 0.5 * temp1 * s->x1mth2 * cos2u;
           su    = su - 0.25 * temp2 * s->x7thm1 * sin2u;
           xnode = s->right_asc_node_lp + 1.5 * temp2 * cosip * sin2u;
    double xinc  = s->inclination_lp    + 1.5 * temp2 * cosip * sinip * cos2u;
    double mvt   = rdotl  - nm * temp1 * s->x1mth2  * sin2u / XKE;
    double rvdot = rvdotl + nm * temp1 * (s->x1mth2 * cos2u
                 + 1.5 * s->con41) / XKE;

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

      printf("---------------------------------------- p8\n");
      printf("px %12.6lf\n", p->x);
      printf("py %12.6lf\n", p->y);
      printf("pz %12.6lf\n", p->z);
      printf("vx %12.6lf\n", v->x);
      printf("vy %12.6lf\n", v->y);
      printf("vz %12.6lf\n", v->z);
    }
  }

  // Satellite decayed?
  if (mrt < 1.0)
  {
    return 6;
  }

  return 0;
}

/*
 * Get position and velocity vectors in the TEME frame at given unix time
 *
 * Inputs:  sat       - orbit struct pointer with initialized orbital data
 *          time      - UTC time to find satellite data at
 *          msec      - Millisecond portion of time
 *          maxiter   - Kepler's equation maximum iteration count
 *          tolerance - Kepler's equation desired precision tolerance
 * Outputs: pos       - 3D position vector in TEME frame in km
 *          vel       - 3D velocity vector in TEME frame in km/sec
 * Returns: 0         - Success
 *         -1         - Invalid inputs or parametres
 *          1         - Mean motion zero or less
 *          2         - Nonsensical orbital eccentricity (e >= 1; e < -0.00001)
 *          3         - Long periodics result error
 *          4         - Short period preliminary quantities error
 *          5         - Decaying satellite
 */
int
sat_get_teme_at
(
    sat* sat,
    time_t* time,
    unsigned int msec,
    unsigned int maxiter,
    double tolerance,
    vec3* pos,
    vec3* vel
)
{
  if ((sat == NULL) ||
      (time == NULL) ||
      (msec >= 1000) ||
      (maxiter < 1) ||
      (tolerance <= 0.0) ||
      (pos == NULL) ||
      (vel == NULL))
  {
    return -1;
  }

  double tdelta = difftime(*time + msec / 1000,
                           sat->epoch + sat->epoch_ms / 1000) / 60.0L;

  return sat_propagate(sat, tdelta, maxiter, tolerance, pos, vel);
}


//ong period periodic contributions to the mean elements
void
dpper(sat* s, double tdelta) // TODO: Rename
{
  double alfdp, betdp, cosip, cosop, dalf, dbet, dls,
  f2,    f3,    pe,    pgh,   ph,   pinc, pl ,
  sel,   ses,   sghl,  sghs,  shll, shs,  sil,
  sinip, sinop, sinzf, sis,   sll,  sls,  xls,
  xnoh,  zf,    zm,    zel,   zes,  znl,  zns;

  // Constants
  zns   = 1.19459e-5; // TODO: const keyword?
  zes   = 0.01675;
  znl   = 1.5835218e-4;
  zel   = 0.05490;

  // Calculate time varying periodics
  zm    = s->zmos + zns * tdelta;
  zf    = zm + 2 * zes * sin(zm);
  sinzf = sin(zf);
  f2    =  0.5 * sinzf * sinzf - 0.25;
  f3    = -0.5 * sinzf * cos(zf);
  ses   = s->se2* f2 + s->se3 * f3;
  sis   = s->si2 * f2 + s->si3 * f3;
  sls   = s->sl2 * f2 + s->sl3 * f3 + s->sl4 * sinzf;
  sghs  = s->sgh2 * f2 + s->sgh3 * f3 + s->sgh4 * sinzf;
  shs   = s->sh2 * f2 + s->sh3 * f3;
  zm    = s->zmol + znl * tdelta;

  zf    = zm + 2 * zel * sin(zm);
  sinzf = sin(zf);
  f2    =  0.5 * sinzf * sinzf - 0.25;
  f3    = -0.5 * sinzf * cos(zf);
  sel   = s->ee2 * f2 + s->e3 * f3;
  sil   = s->xi2 * f2 + s->xi3 * f3;
  sll   = s->xl2 * f2 + s->xl3 * f3 + s->xl4 * sinzf;
  sghl  = s->xgh2 * f2 + s->xgh3 * f3 + s->xgh4 * sinzf;
  shll  = s->xh2 * f2 + s->xh3 * f3;
  pe    = ses + sel;
  pinc  = sis + sil;
  pl    = sls + sll;
  pgh   = sghs + sghl;
  ph    = shs + shll;


  pe    = pe - s->peo;
  pinc  = pinc - s->pinco;
  pl    = pl - s->plo;
  pgh   = pgh - s->pgho;
  ph    = ph - s->pho;
  s->inclination_lp  = s->inclination + pinc;
  s->eccentricity_lp = s->eccentricity + pe;
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
    sinop  = sin(s->right_asc_node_lp);
    cosop  = cos(s->right_asc_node_lp);
    alfdp  = sinip * sinop;
    betdp  = sinip * cosop;
    dalf   =  ph * cosop + pinc * cosip * sinop;
    dbet   = -ph * sinop + pinc * cosip * cosop;
    alfdp  = alfdp + dalf;
    betdp  = betdp + dbet;
    s->right_asc_node_lp  = fmod(s->right_asc_node_lp, TWOPI);
    xls    = s->mean_anomaly_lp + s->argument_perigee_lp + cosip * s->right_asc_node_lp;
    dls    = pl + pgh - pinc * s->right_asc_node_lp * sinip;
    xls    = xls + dls;
    xnoh   = s->right_asc_node_lp;
    s->right_asc_node_lp  = atan2(alfdp, betdp);

    if (fabs(xnoh - s->right_asc_node_lp) > PI)
    {
      if (s->right_asc_node_lp < xnoh)
      {
        s->right_asc_node_lp = s->right_asc_node_lp + TWOPI;
      }
      else
      {
        s->right_asc_node_lp = s->right_asc_node_lp - TWOPI;
      }
    }

    s->mean_anomaly_lp    = s->mean_anomaly_lp + pl;
    s->argument_perigee_lp = xls - s->mean_anomaly_lp - cosip * s->right_asc_node_lp;
  }
}

double  sgn
(
    double x
)
{
  if (x < 0.0)
  {
    return -1.0;
  }
  else
  {
    return 1.0;
  }

}

double  mag2
(
    double x[3]
)
{
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

void    cross2
(
    double vec1[3], double vec2[3], double outvec[3]
)
{
  outvec[0]= vec1[1]*vec2[2] - vec1[2]*vec2[1];
  outvec[1]= vec1[2]*vec2[0] - vec1[0]*vec2[2];
  outvec[2]= vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

double  dot2
(
    double x[3], double y[3]
)
{
  return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);
}  // end dot

double  angle
(
    double vec1[3],
    double vec2[3]
)
{
  double small, undefined, magv1, magv2, temp;
  small     = 0.00000001;
  undefined = 999999.1;

  magv1 = mag2(vec1);
  magv2 = mag2(vec2);

  if (magv1*magv2 > small*small)
  {
    temp= dot2(vec1,vec2) / (magv1*magv2);
    if (fabs( temp ) > 1.0)
      temp= sgn(temp) * 1.0;
    return acos( temp );
  }
  else
    return undefined;
}  // end angle

double  asinh
(
    double xval
)
{
  return log( xval + sqrt( xval*xval + 1.0 ) );
}  // end asinh

/* -----------------------------------------------------------------------------
*
*                           function newtonnu
*
*  this function solves keplers equation when the true anomaly is known.
*    the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
*    the parabolic limit at 168ø is arbitrary. the hyperbolic anomaly is also
*    limited. the hyperbolic sine is used because it's not double valued.
*
*  author        : david vallado                  719-573-2600   27 may 2002
*
*  revisions
*    vallado     - fix small                                     24 sep 2002
*
*  inputs          description                    range / units
*    ecc         - eccentricity                   0.0  to
*    nu          - true anomaly                   -2pi to 2pi rad
*
*  outputs       :
*    e0          - eccentric anomaly              0.0  to 2pi rad       153.02 ø
*    m           - mean anomaly                   0.0  to 2pi rad       151.7425 ø
*
*  locals        :
*    e1          - eccentric anomaly, next value  rad
*    sine        - sine of e
*    cose        - cosine of e
*    ktr         - index
*
*  coupling      :
*    asinh       - arc hyperbolic sine
*
*  references    :
*    vallado       2007, 85, alg 5
* --------------------------------------------------------------------------- */

void newtonnu
(
    double ecc, double nu,
    double* e0, double* m
)
{
  double small, sine, cose;

  // ---------------------  implementation   ---------------------
  *e0= 999999.9;
  *m = 999999.9;
  small = 0.00000001;

  // --------------------------- circular ------------------------
  if ( fabs( ecc ) < small  )
  {
    *m = nu;
    *e0= nu;
  }
  else
    // ---------------------- elliptical -----------------------
    if ( ecc < 1.0-small  )
    {
      sine= ( sqrt( 1.0 -ecc*ecc ) * sin(nu) ) / ( 1.0 +ecc*cos(nu) );
      cose= ( ecc + cos(nu) ) / ( 1.0  + ecc*cos(nu) );
      *e0  = atan2( sine,cose );
      *m   = *e0 - ecc*sin(*e0);
    }
    else
      // -------------------- hyperbolic  --------------------
      if ( ecc > 1.0 + small  )
      {
        if ((ecc > 1.0 ) && (fabs(nu)+0.00001 < PI-acos(1.0 /ecc)))
        {
          sine= ( sqrt( ecc*ecc-1.0  ) * sin(nu) ) / ( 1.0  + ecc*cos(nu) );
          *e0  = asinh( sine );
          *m   = ecc*sinh(*e0) - *e0;
        }
      }
      else
        // ----------------- parabolic ---------------------
        if ( fabs(nu) < 168.0*PI/180.0  )
        {
          *e0= tan( nu*0.5  );
          *m = *e0 + (*e0**e0**e0)/3.0;
        }

  if ( ecc < 1.0  )
  {
    *m = fmod( *m,2.0 *PI );
    if ( *m < 0.0  )
      *m = *m + 2.0 *PI;
    *e0 = fmod( *e0,2.0 *PI );
  }
}  // end newtonnu

/* -----------------------------------------------------------------------------
*
*                           function rv2coe
*
*  this function finds the classical orbital elements given the geocentric
*    equatorial position and velocity vectors.
*
*  author        : david vallado                  719-573-2600   21 jun 2002
*
*  revisions
*    vallado     - fix special cases                              5 sep 2002
*    vallado     - delete extra check in inclination code        16 oct 2002
*    vallado     - add constant file use                         29 jun 2003
*    vallado     - add mu                                         2 apr 2007
*
*  inputs          description                    range / units
*    r           - ijk position vector            km
*    v           - ijk velocity vector            km / s
*    mu          - gravitational parameter        km3 / s2
*
*  outputs       :
*    p           - semilatus rectum               km
*    a           - semimajor axis                 km
*    ecc         - eccentricity
*    incl        - inclination                    0.0  to pi rad
*    omega       - longitude of ascending node    0.0  to 2pi rad
*    argp        - argument of perigee            0.0  to 2pi rad
*    nu          - true anomaly                   0.0  to 2pi rad
*    m           - mean anomaly                   0.0  to 2pi rad
*    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
*    truelon     - true longitude            (ce) 0.0  to 2pi rad
*    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
*
*  locals        :
*    hbar        - angular momentum h vector      km2 / s
*    ebar        - eccentricity     e vector
*    nbar        - line of nodes    n vector
*    c1          - v**2 - u/r
*    rdotv       - r dot v
*    hk          - hk unit vector
*    sme         - specfic mechanical energy      km2 / s2
*    i           - index
*    e           - eccentric, parabolic,
*                  hyperbolic anomaly             rad
*    temp        - temporary variable
*    typeorbit   - type of orbit                  ee, ei, ce, ci
*
*  coupling      :
*    mag         - magnitude of a vector
*    cross       - cross product of two vectors
*    angle       - find the angle between two vectors
*    newtonnu    - find the mean anomaly
*
*  references    :
*    vallado       2007, 126, alg 9, ex 2-5
* --------------------------------------------------------------------------- */

void rv2coe
     (
       double r[3], double v[3],
       double* p, double* a, double* ecc, double* incl, double* omega, double* argp,
       double* nu, double* m, double* arglat, double* truelon, double* lonper
     )
     {
       double undefined, small, hbar[3], nbar[3], magr, magv, magn, ebar[3], sme,
              rdotv, infinite, temp, c1, hk, twopi, magh, halfpi, e;

       int i;
       char typeorbit[3];

     twopi  = 2.0 * PI;
     halfpi = 0.5 * PI;
     small  = 0.00000001;
     undefined = 999999.1;
     infinite  = 999999.9;

     // -------------------------  implementation   -----------------
     magr = mag2( r );
     magv = mag2( v );

     // ------------------  find h n and e vectors   ----------------
     cross2( r,v, hbar );
     magh = mag2( hbar );
     if ( magh > small )
       {
         nbar[0]= -hbar[1];
         nbar[1]=  hbar[0];
         nbar[2]=   0.0;
         magn = mag2( nbar );
         c1 = magv*magv - GM  / magr;
         rdotv = dot2( r,v );
         for (i= 0; i <= 2; i++)
             ebar[i]= (c1*r[i] - rdotv*v[i]) / GM;
         *ecc = mag2( ebar );

         // ------------  find *a e and semi-latus rectum   ----------
         sme= ( magv*magv*0.5  ) - ( GM  / magr );
         if ( fabs( sme ) > small )
             *a= -GM  / (2.0 *sme);
           else
             *a= infinite;
         *p = magh*magh / GM;

         // -----------------  find inclination   -------------------
         hk= hbar[2] / magh;
         *incl= acos( hk );

         // --------  determine type of orbit for later use  --------
         // ------ elliptical, parabolic, hyperbolic inclined -------
         strcpy(typeorbit,"ei");
         if ( *ecc < small )
           {
             // ----------------  circular equatorial ---------------
             if  ((*incl<small) | (fabs(*incl-PI)<small))
                 strcpy(typeorbit,"ce");
               else
                 // --------------  circular inclined ---------------
                 strcpy(typeorbit,"ci");
           }
           else
           {
             // - elliptical, parabolic, hyperbolic equatorial --
             if  ((*incl<small) | (fabs(*incl-PI)<small))
                 strcpy(typeorbit,"ee");
           }

         // ----------  find longitude of ascending node ------------
         if ( magn > small )
           {
             temp= nbar[0] / magn;
             if ( fabs(temp) > 1.0  )
                 temp= sgn(temp);
             *omega= acos( temp );
             if ( nbar[1] < 0.0  )
                 *omega= twopi - *omega;
           }
           else
             *omega= undefined;

         // ---------------- find argument of perigee ---------------
         if ( strcmp(typeorbit,"ei") == 0 )
           {
             *argp = angle( nbar,ebar);
             if ( ebar[2] < 0.0  )
                 *argp= twopi - *argp;
           }
           else
             *argp= undefined;

         // ------------  find true anomaly at epoch    -------------
         if ( typeorbit[0] == 'e' )
           {
             *nu =  angle( ebar,r);
             if ( rdotv < 0.0  )
                 *nu= twopi - *nu;
           }
           else
             *nu= undefined;

         // ----  find argument of latitude - circular inclined -----
         if ( strcmp(typeorbit,"ci") == 0 )
           {
             *arglat = angle( nbar,r );
             if ( r[2] < 0.0  )
                 *arglat= twopi - *arglat;
             *m = *arglat;
           }
           else
             *arglat= undefined;

         // -- find longitude of perigee - elliptical equatorial ----
         if  (( *ecc>small ) && (strcmp(typeorbit,"ee") == 0))
           {
             temp= ebar[0] / *ecc;
             if ( fabs(temp) > 1.0  )
                 temp= sgn(temp);
             *lonper= acos( temp );
             if ( ebar[1] < 0.0  )
                 *lonper= twopi - *lonper;
             if ( *incl > halfpi )
                 *lonper= twopi - *lonper;
           }
           else
             *lonper= undefined;

         // -------- find true longitude - circular equatorial ------
         if  (( magr>small ) && ( strcmp(typeorbit,"ce") == 0 ))
           {
             temp= r[0] / magr;
             if ( fabs(temp) > 1.0  )
                 temp= sgn(temp);
             *truelon= acos( temp );
             if ( r[1] < 0.0  )
                 *truelon= twopi - *truelon;
             if ( *incl > halfpi )
                 *truelon= twopi - *truelon;
             *m = *truelon;
           }
           else
             *truelon= undefined;

         // ------------ find mean anomaly for all orbits -----------
         if ( typeorbit[0] == 'e' )
             newtonnu(*ecc,*nu,  &e, m);
     }
      else
     {
        *p    = undefined;
        *a    = undefined;
        *ecc  = undefined;
        *incl = undefined;
        *omega= undefined;
        *argp = undefined;
        *nu   = undefined;
        *m    = undefined;
        *arglat = undefined;
        *truelon= undefined;
        *lonper = undefined;
     }
   }  // end rv2coe
