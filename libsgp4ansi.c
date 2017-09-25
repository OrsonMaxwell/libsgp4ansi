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

const char library_version[] = VERSION;

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

void sat_print(sat* sat, char* caption)
{
  FILE *f = fopen("sat.txt", "a");
  fprintf(f, "===== %s =====\n", caption);
  fprintf(f, "----- TLE portion -----\n");
  fprintf(f, "norad_number    : %d\n", sat->norad_number);
  fprintf(f, "mean_motion_dt2 : %32.25lf\n", sat->mean_motion_dt2);
  fprintf(f, "mean_motion_ddt6: %32.25lf\n", sat->mean_motion_ddt6);
  fprintf(f, "Bstar           : %32.25lf\n", sat->Bstar);
  fprintf(f, "inclination     : %32.25lf\n", sat->inclination);
  fprintf(f, "right_asc_node  : %32.25lf\n", sat->right_asc_node);
  fprintf(f, "eccentricity    : %32.25lf\n", sat->eccentricity);
  fprintf(f, "argument_perigee: %32.25lf\n", sat->argument_perigee);
  fprintf(f, "mean_anomaly    : %32.25lf\n", sat->mean_anomaly);
  fprintf(f, "mean_motion     : %32.25lf\n", sat->mean_motion);
  fprintf(f, "----- Time -----\n");
  fprintf(f, "epoch:\t\t%s%lfms\n", ctime(&sat->epoch), sat->epoch_ms);
  fprintf(f, "----- Standard orbital elements -----\n");
  fprintf(f, "xnodp           : %32.25lf\n", sat->xnodp);
  fprintf(f, "aodp            : %32.25lf\n", sat->aodp);
  fprintf(f, "perigee         : %32.25lf\n", sat->perigee);
  fprintf(f, "perigee_alt     : %32.25lf\n", sat->perigee_alt);
  fprintf(f, "period          : %32.25lf\n", sat->period);
  /*
  fprintf(f, "julepoch:\t%32.25lf\n", sat->julepoch);
  fprintf(f, "GSTo:\t\t%32.25lf\n", sat->GSTo);
  fprintf(f, "----- Flags -----\n");
  fprintf(f, "Deepsp:\t\t%d\n", sat->is_deep_space);
  fprintf(f, "Loworb:\t\t%d\n", sat->use_simple_model);
  fprintf(f, "----- Standard terms -----\n");
  fprintf(f, "a:\t\t%32.25lf\n", sat->a);
  fprintf(f, "altapoR:\t\t%32.25lf\n", sat->altapoR);
  fprintf(f, "altperR:\t\t%32.25lf\n", sat->altperR);
  fprintf(f, "aycof:\t\t%32.25lf\n", sat->aycof);
  fprintf(f, "C1:\t\t%32.25lf\n", sat->C1);
  fprintf(f, "C4:\t\t%32.25lf\n", sat->C4);
  fprintf(f, "C5:\t\t%32.25lf\n", sat->C5);
  fprintf(f, "con41:\t\t%32.25lf\n", sat->x3thm1);
  fprintf(f, "d2:\t\t%32.25lf\n", sat->d2);
  fprintf(f, "d3:\t\t%32.25lf\n", sat->d3);
  fprintf(f, "d4:\t\t%32.25lf\n", sat->d4);
  fprintf(f, "delMo:\t\t%32.25lf\n", sat->delMo);
  fprintf(f, "eta:\t\t%32.25lf\n", sat->eta);
  fprintf(f, "mdot:\t\t%32.25lf\n", sat->mdot);
  fprintf(f, "nodecf:\t\t%32.25lf\n", sat->nodecf);
  fprintf(f, "nodedot:\t%32.25lf\n", sat->nodedot);
  fprintf(f, "omegaprime:\t%32.25lf\n", sat->omegaprime);
  fprintf(f, "omgcof:\t\t%32.25lf\n", sat->omgcof);
  fprintf(f, "sinMo:\t\t%32.25lf\n", sat->sinMo);
  fprintf(f, "t2cof:\t\t%32.25lf\n", sat->t2cof);
  fprintf(f, "t3cof:\t\t%32.25lf\n", sat->t3cof);
  fprintf(f, "t4cof:\t\t%32.25lf\n", sat->t4cof);
  fprintf(f, "t5cof:\t\t%32.25lf\n", sat->t5cof);
  fprintf(f, "x1mth2:\t\t%32.25lf\n", sat->x1mth2);
  fprintf(f, "x7thm1:\t\t%32.25lf\n", sat->x7thm1);
  fprintf(f, "xlcof:\t\t%32.25lf\n", sat->xlcof);
  fprintf(f, "xmcof:\t\t%32.25lf\n", sat->xmcof);
  fprintf(f, "----- Deepspace terms -----\n");
  fprintf(f, "e3:\t\t%32.25lf\n", sat->e3);
  fprintf(f, "ee2:\t\t%32.25lf\n", sat->ee2);
  fprintf(f, "peo:\t\t%32.25lf\n", sat->peo);
  fprintf(f, "pgho:\t\t%32.25lf\n", sat->pgho);
  fprintf(f, "pho:\t\t%32.25lf\n", sat->pho);
  fprintf(f, "pinco:\t\t%32.25lf\n", sat->pinco);
  fprintf(f, "plo:\t\t%32.25lf\n", sat->plo);
  fprintf(f, "se2:\t\t%32.25lf\n", sat->se2);
  fprintf(f, "se3:\t\t%32.25lf\n", sat->se3);
  fprintf(f, "sgh2:\t\t%32.25lf\n", sat->sgh2);
  fprintf(f, "sgh3:\t\t%32.25lf\n", sat->sgh3);
  fprintf(f, "sgh4:\t\t%32.25lf\n", sat->sgh4);
  fprintf(f, "sh2:\t\t%32.25lf\n", sat->sh2);
  fprintf(f, "sh3:\t\t%32.25lf\n", sat->sh3);
  fprintf(f, "si2:\t\t%32.25lf\n", sat->si2);
  fprintf(f, "si3:\t\t%32.25lf\n", sat->si3);
  fprintf(f, "sl2:\t\t%32.25lf\n", sat->sl2);
  fprintf(f, "sl3:\t\t%32.25lf\n", sat->sl3);
  fprintf(f, "sl4:\t\t%32.25lf\n", sat->sl4);
  fprintf(f, "xgh2:\t\t%32.25lf\n", sat->xgh2);
  fprintf(f, "xgh3:\t\t%32.25lf\n", sat->xgh3);
  fprintf(f, "xgh4:\t\t%32.25lf\n", sat->xgh4);
  fprintf(f, "xh2:\t\t%32.25lf\n", sat->xh2);
  fprintf(f, "xh3:\t\t%32.25lf\n", sat->xh3);
  fprintf(f, "xi2:\t\t%32.25lf\n", sat->xi2);
  fprintf(f, "xi3:\t\t%32.25lf\n", sat->xi3);
  fprintf(f, "xl2:\t\t%32.25lf\n", sat->xl2);
  fprintf(f, "xl3:\t\t%32.25lf\n", sat->xl3);
  fprintf(f, "xl4:\t\t%32.25lf\n", sat->xl4);
  fprintf(f, "zmol:\t\t%32.25lf\n", sat->zmol);
  fprintf(f, "zmos:\t\t%32.25lf\n", sat->zmos);
  fprintf(f, "----- Resonant terms -----\n");
  fprintf(f, "d2201:\t\t%32.25lf\n", sat->d2201);
  fprintf(f, "d2211:\t\t%32.25lf\n", sat->d2211);
  fprintf(f, "d3210:\t\t%32.25lf\n", sat->d3210);
  fprintf(f, "d3222:\t\t%32.25lf\n", sat->d3222);
  fprintf(f, "d4410:\t\t%32.25lf\n", sat->d4410);
  fprintf(f, "d4422:\t\t%32.25lf\n", sat->d4422);
  fprintf(f, "d5220:\t\t%32.25lf\n", sat->d5220);
  fprintf(f, "d5232:\t\t%32.25lf\n", sat->d5232);
  fprintf(f, "d5421:\t\t%32.25lf\n", sat->d5421);
  fprintf(f, "d5433:\t\t%32.25lf\n", sat->d5433);
  fprintf(f, "dedt:\t\t%32.25lf\n", sat->dedt);
  fprintf(f, "didt:\t\t%32.25lf\n", sat->didt);
  fprintf(f, "dmdt:\t\t%32.25lf\n", sat->dmdt);
  fprintf(f, "dnodt:\t\t%32.25lf\n", sat->dnodt);
  fprintf(f, "domdt:\t\t%32.25lf\n", sat->domdt);
  fprintf(f, "del1:\t\t%32.25lf\n", sat->del1);
  fprintf(f, "del2:\t\t%32.25lf\n", sat->del2);
  fprintf(f, "del3:\t\t%32.25lf\n", sat->del3);
  fprintf(f, "xfact:\t\t%32.25lf\n", sat->xfact);
  fprintf(f, "xlamo:\t\t%32.25lf\n", sat->xlamo);
  fprintf(f, "xli:\t\t%32.25lf\n", sat->xli);
  fprintf(f, "xni:\t\t%32.25lf\n", sat->xni);
  fprintf(f, "----- Legacy -----\n");
  fprintf(f, "isimp:\t\t%d\n", sat->isimp);
  fprintf(f, "error:\t\t%d\n", sat->error);
  fprintf(f, "method:\t\t%c\n", sat->method);
  fprintf(f, "operationmode:\t\t%c\n", sat->operationmode);
  fprintf(f, "cc1:\t\t%32.25lf\n", sat->cc1);
  fprintf(f, "cc4:\t\t%32.25lf\n", sat->cc4);
  fprintf(f, "cc5:\t\t%32.25lf\n", sat->cc5);
  fprintf(f, "delmo:\t\t%32.25lf\n", sat->delmo);
  fprintf(f, "argpdot:\t\t%32.25lf\n", sat->argpdot);
  fprintf(f, "sinmao:\t\t%32.25lf\n", sat->sinmao);
  fprintf(f, "t:\t\t%32.25lf\n", sat->t);
  fprintf(f, "irez:\t\t%d\n", sat->irez);
  fprintf(f, "gsto:\t\t%32.25lf\n", sat->gsto);
  fprintf(f, "atime:\t\t%32.25lf\n", sat->atime);
  fprintf(f, "bstar:\t\t%32.25lf\n", sat->bstar);
  fprintf(f, "ecco:\t\t%32.25lf\n", sat->ecco);
  fprintf(f, "argpo:\t\t%32.25lf\n", sat->argpo);
  fprintf(f, "inclo:\t\t%32.25lf\n", sat->inclo);
  fprintf(f, "mo:\t\t%32.25lf\n", sat->mo);
  fprintf(f, "nodeo:\t\t%32.25lf\n", sat->nodeo);
  fprintf(f, "alta:\t\t%32.25lf\n", sat->alta);
  fprintf(f, "altp:\t\t%32.25lf\n", sat->altp);*/
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
  double eosq    = pow(s->eccentricity, 2);
  double theta2  = pow(cosio, 2);
  double x3thm1  = 3 * theta2 - 1;
  double betao2  = 1 - eosq;
  double betao   = sqrt(betao2);

  // Un-Kozai mean motion
  double delta1  = 0.75 * J2 * (3.0 * theta2 - 1.0) / (betao * betao2);
  double a0      = aa * (1.0 - pow(delta1 / pow(aa, 2), 2)
                 - delta1 / pow(aa, 2) * (1.0 / 3.0 + 134.0
                 * pow(delta1 / pow(aa, 2), 2) / 81.0));
  double delta0  = delta1 / pow(a0, 2);
  s->xnodp       = s->mean_motion / (1 + delta0);
  s->aodp        = a0 / (1 - delta0);
  s->perigee     = (s->aodp * (1 - s->eccentricity)) * RE; // TODO: Remove?
  s->perigee_alt = s->perigee - RE;
  s->period      = TWOPI / s->xnodp;

  //sat_print(s, "INIT");

  if (s->period >= 225)
  {
    s->is_deep_space      = true;
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
                 + 0.375 * J2 * tsi / psi2 * x3thm1
                 * (8 + 3 * eta2 * (8 + eta2)));
  double C1      = s->Bstar * C2;
  double x1mth2  = 1 - theta2;
  double C4      = 2 * s->xnodp * coef1 * s->aodp * betao2
                 * (s->eta * (2 + 0.5 * eta2) + s->eccentricity
                 * (0.5 + 2 * eta2) - J2 * tsi / (s->aodp * psi2)
                 * (-3 * x3thm1 * (1 - 2 * eeta + eta2
                 * (1.5 - 0.5 * eeta)) + 0.75 * x1mth2
                 * (2 * eta2 - eeta * (1 + eta2))
                 * cos(2 * s->argument_perigee)));
  double theta4  = pow(cosio, 4);
  double temp1   = 1.5 * J2 * pinv2 * s->xnodp;
  double temp2   = 0.5 * temp1 * J2 * pinv2;
  double temp3   = -0.46875 * J4 * pow(pinv2, 2) * s->xnodp;
      s->xmdot   = s->xnodp + 0.5 * temp1 * betao * x3thm1 + 0.0625
                 * temp2 * betao * (13 - 78 * theta2 + 137 * theta4);
  double x1m5th  = 1 - 5 * theta2;
      s->omgdot  = -0.5 * temp1 * x1m5th
                 + 0.0625 * temp2 * (7 - 114 * theta2 + 395 * theta4)
                 + temp3 * (3 - 36 * theta2 + 49 * theta4);
  double xhdot1  = -temp1 * cosio;
      s->xnodot  = xhdot1 + (0.5 * temp2 * (4 - 19 * theta2)
                 + 2 * temp3 * (3 - 7 * theta2)) * cosio;
      s->xnodcf  = 3.5 * betao * xhdot1 * C1;
  double t2cof   = 1.5 * C1; // TODO: Remove?
  // Division by zero check then inclination = 180 deg
  double xlcof = 0.125 * A3OVK2 * sinio * (3.0 + 5.0 * cosio)
            / ((fabs(cosio+1.0) > 1.5e-12)?((1.0 + cosio)):(1.5e-12));
  double aycof   = 0.25 * A3OVK2 * sinio;
  double x7thm1  = 7 * theta2 - 1;

  if (s->is_deep_space == true) // Deep space init here
  {
    //double tc    =  0.0;
    //double inclm = s->inclination;

    double GSTo  = jul2gst(s->julian_epoch);

    // Constants
    const double zes     =  0.01675;
    const double zel     =  0.05490;
    const double c1ss    =  2.9864797e-6;
    const double c1l     =  4.7968065e-7;
    const double zsinis  =  0.39785416;
    const double zcosis  =  0.91744867;
    const double zcosgs  =  0.1945905;
    const double zsings  = -0.98088458;

    double snodm  = sin(s->right_asc_node); // TODO: duplicates?
    double cnodm  = cos(s->right_asc_node); // TODO: duplicates?
    double sinomm = sin(s->argument_perigee); // TODO: duplicates?
    double cosomm = cos(s->argument_perigee); // TODO: duplicates?
    double sinim  = sin(s->inclination); // TODO: duplicates?
    double cosim  = cos(s->inclination); // TODO: duplicates?

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
    // printf("day:    %22.15lf\n", day);
    // printf("xnodce: %22.15lf\n", xnodce);
    // printf("stem:   %22.15lf\n", stem);
    // printf("ctem:   %22.15lf\n", ctem);
    // printf("zcosil: %22.15lf\n", zcosil);
    // printf("zsinil: %22.15lf\n", zsinil);
    // printf("zsinhl: %22.15lf\n", zsinhl);
    // printf("zcoshl: %22.15lf\n", zcoshl);
    // printf("gam:    %22.15lf\n", gam);
    // printf("zy:     %22.15lf\n", zy);
    // printf("zx:     %22.15lf\n", zx);
    // printf("zcosgl: %22.15lf\n", zcosgl);
    // printf("zsingl: %22.15lf\n", zsingl);

    s->zmos = fmod(6.2565837 + 0.017201977 * day, TWOPI);
//    printf("zmos:   %22.15lf\n", zmos);

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
   // printf("zcosh:  %22.15lf\n", zcosh);
   // printf("zsinh:  %22.15lf\n", zsinh);
   // printf("xnoi:   %22.15lf\n", xnoi);

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
//      printf("a1:     %22.15lf\n", a1);
//      printf("a3:     %22.15lf\n", a3);
//      printf("a7:     %22.15lf\n", a7);
//      printf("a8:     %22.15lf\n", a8);
//      printf("a9:     %22.15lf\n", a9);
//      printf("a10:    %22.15lf\n", a10);
//      printf("a2:     %22.15lf\n", a2);
//      printf("a4:     %22.15lf\n", a4);
//      printf("a5:     %22.15lf\n", a5);
//      printf("a6:     %22.15lf\n", a6);

      x1  =  a1 * cosomm + a2 * sinomm;
      x2  =  a3 * cosomm + a4 * sinomm;
      x3  = -a1 * sinomm + a2 * cosomm;
      x4  = -a3 * sinomm + a4 * cosomm;
      x5  =  a5 * sinomm;
      x6  =  a6 * sinomm;
      x7  =  a5 * cosomm;
      x8  =  a6 * cosomm;
//      printf("x1:     %22.15lf\n", x1);
//      printf("x2:     %22.15lf\n", x2);
//      printf("x3:     %22.15lf\n", x3);
//      printf("x4:     %22.15lf\n", x4);
//      printf("x5:     %22.15lf\n", x5);
//      printf("x6:     %22.15lf\n", x6);
//      printf("x7:     %22.15lf\n", x7);
//      printf("x8:     %22.15lf\n", x8);

      z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
      z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
      z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
      z1  =  3.0 *  (a1 * a1 + a2 * a2) + z31 * eosq;
      z2  =  6.0 *  (a1 * a3 + a2 * a4) + z32 * eosq;
      z3  =  3.0 *  (a3 * a3 + a4 * a4) + z33 * eosq;
      z11 = -6.0 * a1 * a5 + eosq *  (-24.0 * x1 * x7-6.0 * x3 * x5);
      z12 = -6.0 *  (a1 * a6 + a3 * a5) + eosq *
          (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
      z13 = -6.0 * a3 * a6 + eosq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
      z21 =  6.0 * a2 * a5 + eosq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
      z22 =  6.0 *  (a4 * a5 + a2 * a6) + eosq *
          (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
      z23 =  6.0 * a4 * a6 + eosq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
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
//      printf("z31:    %22.15lf\n", z31);
//      printf("z32:    %22.15lf\n", z32);
//      printf("z33:    %22.15lf\n", z33);
//      printf("z1:     %22.15lf\n", z1);
//      printf("z2:     %22.15lf\n", z2);
//      printf("z3:     %22.15lf\n", z3);
//      printf("z11:    %22.15lf\n", z11);
//      printf("z12:    %22.15lf\n", z12);
//      printf("z13:    %22.15lf\n", z13);
//      printf("z21:    %22.15lf\n", z21);
//      printf("z22:    %22.15lf\n", z22);
//      printf("z23:    %22.15lf\n", z23);
//      printf("s3:     %22.15lf\n", s3);
//      printf("s2:     %22.15lf\n", s2);
//      printf("s4:     %22.15lf\n", s4);
//      printf("s1:     %22.15lf\n", s1);
//      printf("s5:     %22.15lf\n", s5);
//      printf("s6:     %22.15lf\n", s6);
//      printf("s7:     %22.15lf\n", s7);

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
//    printf("zmol:   %22.15lf\n", zmol);
//    printf("------------------------------\n");
    // Do final solar terms
    s->se2  =   2 * ss1 * ss6;
    s->se3  =   2 * ss1 * ss7;
    s->si2  =   2 * ss2 * sz12;
    s->si3  =   2 * ss2 * (sz13 - sz11);
    s->sl2  =  -2 * ss3 * sz2;
    s->sl3  =  -2 * ss3 * (sz3 - sz1);
    s->sl4  =  -2 * ss3 * (-21 - 9 * eosq) * zes;
    s->sgh2 =   2 * ss4 * sz32;
    s->sgh3 =   2 * ss4 * (sz33 - sz31);
    s->sgh4 = -18 * ss4 * zes;
    s->sh2  =  -2 * ss2 * sz22;
    s->sh3  =  -2 * ss2 * (sz23 - sz21);
//    printf("se2:    %22.15lf\n", s->se2);
//    printf("se3:    %22.15lf\n", s->se3);
//    printf("si2:    %22.15lf\n", s->si2);
//    printf("si3:    %22.15lf\n", s->si3);
//    printf("sl2:    %22.15lf\n", s->sl2);
//    printf("sl3:    %22.15lf\n", s->sl3);
//    printf("sl4:    %22.15lf\n", s->sl4);
//    printf("sgh2:   %22.15lf\n", s->sgh2);
//    printf("sgh3:   %22.15lf\n", s->sgh3);
//    printf("sgh4:   %22.15lf\n", s->sgh4);
//    printf("sh2:    %22.15lf\n", s->sh2);
//    printf("sh3:    %22.15lf\n", s->sh3);

    // Do final lunar terms
    s->ee2  =   2 * s1 * s6;
    s->e3   =   2 * s1 * s7;
    s->xi2  =   2 * s2 * z12;
    s->xi3  =   2 * s2 * (z13 - z11);
    s->xl2  =  -2 * s3 * z2;
    s->xl3  =  -2 * s3 * (z3 - z1);
    s->xl4  =  -2 * s3 * (-21 - 9 * eosq) * zel;
    s->xgh2 =   2 * s4 * z32;
    s->xgh3 =   2 * s4 * (z33 - z31);
    s->xgh4 = -18 * s4 * zel;
    s->xh2  =  -2 * s2 * z22;
    s->xh3  =  -2 * s2 * (z23 - z21);
//    printf("ee2:    %22.15lf\n", s->ee2);
//    printf("e3:     %22.15lf\n", s->e3);
//    printf("xi2:    %22.15lf\n", s->xi2);
//    printf("xi3:    %22.15lf\n", s->xi3);
//    printf("xl2:    %22.15lf\n", s->xl2);
//    printf("xl3:    %22.15lf\n", s->xl3);
//    printf("xl4:    %22.15lf\n", s->xl4);
//    printf("xgh2:   %22.15lf\n", s->xgh2);
//    printf("xgh3:   %22.15lf\n", s->xgh3);
//    printf("xgh4:   %22.15lf\n", s->xgh4);
//    printf("xh2:    %22.15lf\n", s->xh2);
//    printf("xh3:    %22.15lf\n", s->xh3);
//
//    printf("------------------------------\n");
    // dpper(s, 0, true); TODO: Investigate further


    /*dsinit
    (
        whichconst,
        cosim, emsq, satrec.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
        ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, satrec.t, tc,
        satrec.gsto, satrec.mo, satrec.mdot, satrec.no, satrec.nodeo,
        satrec.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
        satrec.ecco, eccsq, em, argpm, inclm, mm, nm, nodem,
        satrec.irez,  satrec.atime,
        satrec.d2201, satrec.d2211, satrec.d3210, satrec.d3222 ,
        satrec.d4410, satrec.d4422, satrec.d5220, satrec.d5232,
        satrec.d5421, satrec.d5433, satrec.dedt,  satrec.didt,
        satrec.dmdt,  dndt,         satrec.dnodt, satrec.domdt ,
        satrec.del1,  satrec.del2,  satrec.del3,  satrec.xfact,
        satrec.xlamo, satrec.xli,   satrec.xni, opsmode
    );*/

    /*void dsinit
         (
           int whichconst,
           double cosim,  double emsq,   double argpo,   double s1,     double s2,
           double s3,     double s4,     double s5,      double sinim,  double ss1,
           double ss2,    double ss3,    double ss4,     double ss5,    double sz1,
           double sz3,    double sz11,   double sz13,    double sz21,   double sz23,
           double sz31,   double sz33,   double t,       double tc,     double gsto,
           double mo,     double mdot,   double no,      double nodeo,  double nodedot,
           double xpidot, double z1,     double z3,      double z11,    double z13,
           double z21,    double z23,    double z31,     double z33,    double ecco,
           double eccsq,  double* em,    double* argpm,  double* inclm, double* mm,
           double* nm,    double* nodem,
           int* irez,
           double* atime, double* d2201, double* d2211,  double* d3210, double* d3222,
           double* d4410, double* d4422, double* d5220,  double* d5232, double* d5421,
           double* d5433, double* dedt,  double* didt,   double* dmdt,  double* dndt,
           double* dnodt, double* domdt, double* del1,   double* del2,  double* del3,
           double* xfact, double* xlamo, double* xli,    double* xni
         )*/


    // Replace these
    const double q22    = 1.7891679e-6;
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
    double sls  = -zns * ss3 * (sz1 + sz3 - 14 - 6 * eosq);
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

//    printf("ses:    %22.15lf\n", ses);
//    printf("sis:    %22.15lf\n", sis);
//    printf("sls:    %22.15lf\n", sls);
//    printf("sghs:   %22.15lf\n", sghs);
//    printf("shs:    %22.15lf\n", shs);
//    printf("sgs:    %22.15lf\n", sgs);

    // Do lunar terms
    s->dedt = ses + s1 * znl * s5;
    s->didt = sis + s2 * znl * (z11 + z13);
    s->dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * eosq);

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

//    printf("dedt:   %22.15lf\n", s->dedt);
//    printf("didt:   %22.15lf\n", s->didt);
//    printf("dmdt:   %22.15lf\n", s->dmdt);
//    printf("dndt:   %22.15lf\n", s->dndt);
//    printf("sghl:   %22.15lf\n", sghl);
//    printf("shll:   %22.15lf\n", shll);
//    printf("domdt:  %22.15lf\n", s->domdt);
//    printf("dnodt:  %22.15lf\n", s->dnodt);

    // Calculate deep space resonance effects

    /* TODO: t is zero at init, right?
     *em     = *em + *dedt * t;
     *inclm  = *inclm + *didt * t;
     *argpm  = *argpm + *domdt * t;
     *nodem  = *nodem + *dnodt * t;
     *mm     = *mm + *dmdt * t;
     */

//    printf("----------------------------\n");
//    printf("theta:  %22.15lf\n", GSTo);
//    printf("em:     %22.15lf\n", s->eccentricity);
//    printf("inclm:  %22.15lf\n", s->inclination);

    double aonv  = pow(s->xnodp / XKE, TWOTHIRD); // TODO: Duplicate?
    double ainv2 = pow(aonv, 2);

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

        // TODO: Starting from here, the results of the polynomial calculations
        // begin to differ slightly due to numerical garbage in insignificant
        // decimal places of the eccentricity parameter. Requires testing to
        // fully understand the net impact of this on the resulting ephemerides.

        if (s->eccentricity <= 0.65)
        {
          g211 =    3.616  -  13.2470 * s->eccentricity
               +   16.2900 * eosq;
          g310 =  -19.302  + 117.3900 * s->eccentricity
               -  228.4190 * eosq +  156.5910 * eocu;
          g322 =  -18.9068 + 109.7927 * s->eccentricity
               -  214.6334 * eosq +  146.5816 * eocu;
          g410 =  -41.122  + 242.6940 * s->eccentricity
               -  471.0940 * eosq +  313.9530 * eocu;
          g422 = -146.407  + 841.8800 * s->eccentricity
               - 1629.014  * eosq + 1083.4350 * eocu;
          g520 = -532.114  + 3017.977 * s->eccentricity
               - 5740.032  * eosq + 3708.2760 * eocu;
        }
        else
        {
          g211 =   -72.099 +   331.819 * s->eccentricity
               -   508.738 * eosq +   266.724 * eocu;
          g310 =  -346.844 +  1582.851 * s->eccentricity
               -  2415.925 * eosq +  1246.113 * eocu;
          g322 =  -342.585 +  1554.908 * s->eccentricity
               -  2366.899 * eosq +  1215.972 * eocu;
          g410 = -1052.797 +  4758.686 * s->eccentricity
               -  7193.992 * eosq +  3651.957 * eocu;
          g422 = -3581.690 + 16178.110 * s->eccentricity
               - 24462.770 * eosq + 12422.520 * eocu;

          if (s->eccentricity > 0.715)
          {
            g520      = -5149.66 + 29936.92 * s->eccentricity
                      - 54087.36 * eosq + 31324.56 * eocu;
          }
          else
          {
            g520      = 1464.74 - 4664.75 * s->eccentricity + 3763.64 * eosq;
          }
        }
        if (s->eccentricity < 0.7)
        {
          g533 = -919.22770 + 4988.6100 * s->eccentricity
               - 9064.7700  * eosq + 5542.21  * eocu;
          g521 = -822.71072 + 4568.6173 * s->eccentricity
               - 8491.4146  * eosq + 5337.524 * eocu;
          g532 = -853.66600 + 4690.2500 * s->eccentricity
               - 8624.7700  * eosq + 5341.4  * eocu;
        }
        else
        {
          g533 = -37995.780 + 161616.52 * s->eccentricity
               - 229838.20  * eosq + 109377.94 * eocu;
          g521 = -51752.104 + 218913.95 * s->eccentricity
               - 309468.16  * eosq + 146349.42 * eocu;
          g532 = -40023.880 + 170470.89 * s->eccentricity
               - 242699.48  * eosq + 115605.82 * eocu;
        }

//        printf("----------------------------\n");
//        printf("g201:   %22.15lf\n", g201);
//        printf("g211:   %22.15lf\n", g211);
//        printf("g310:   %22.15lf\n", g310);
//        printf("g322:   %22.15lf\n", g322);
//        printf("g410:   %22.15lf\n", g410);
//        printf("g422:   %22.15lf\n", g422);
//        printf("g520:   %22.15lf\n", g520);
//        printf("g533:   %22.15lf\n", g533);
//        printf("g521:   %22.15lf\n", g521);
//        printf("g532:   %22.15lf\n", g532);

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
//        printf("f220:   %22.15lf\n", f220);
//        printf("f221:   %22.15lf\n", f221);
//        printf("f321:   %22.15lf\n", f321);
//        printf("f322:   %22.15lf\n", f322);
//        printf("f441:   %22.15lf\n", f441);
//        printf("f442:   %22.15lf\n", f442);
//        printf("f522:   %22.15lf\n", f522);
//        printf("f523:   %22.15lf\n", f523);
//        printf("f542:   %22.15lf\n", f542);
//        printf("f543:   %22.15lf\n", f543);

               temp1 = 3.0 * pow(s->xnodp, 2) * ainv2;
        double temp  = temp1 * root22;
        double d2201 = temp * f220 * g201;
        double d2211 = temp * f221 * g211;
               temp1 = temp1 * aonv;
               temp  = temp1 * root32;
        double d3210 = temp * f321 * g310;
        double d3222 = temp * f322 * g322;
               temp1 = temp1 * aonv;
               temp  = 2.0 * temp1 * root44;
        double d4410 = temp * f441 * g410;
        double d4422 = temp * f442 * g422;
               temp1 = temp1 * aonv;
               temp  = temp1 * root52;
        double d5220 = temp * f522 * g520;
        double d5232 = temp * f523 * g532;
               temp  = 2.0 * temp1 * root54;
        double d5421 = temp * f542 * g521;
        double d5433 = temp * f543 * g533;
            s->xlamo = fmod(s->mean_anomaly + 2 * s->right_asc_node
                     - 2 * GSTo, TWOPI);
            s->xfact = s->xmdot + s->dmdt
                     + 2 * (s->xnodot + s->dnodt - rptim) - s->xnodp;
        //double em    = emo;
        //double emsq  = emsqo;

//        printf("----------------------------\n");
//        printf("ainv2:  %22.15lf\n", ainv2);
//        printf("d2201:  %22.15lf\n", d2201);
//        printf("d2211:  %22.15lf\n", d2211);
//        printf("d3210:  %22.15lf\n", d3210);
//        printf("d3222:  %22.15lf\n", d3222);
//        printf("d4410:  %22.15lf\n", d4410);
//        printf("d4422:  %22.15lf\n", d4422);
//        printf("d5220:  %22.15lf\n", d5220);
//        printf("d5232:  %22.15lf\n", d5232);
//        printf("d5232:  %22.15lf\n", d5232);
//        printf("d5433:  %22.15lf\n", d5433);
      }

      // Synchronous resonance terms

      double g200, g310, g300, f220, f311, f330, del1, del2, del3;

      if (s->is_24h_resonant == true)
      {
            g200 = 1.0 + eosq * (-2.5 + 0.8125 * eosq);
            g310 = 1.0 + 2.0 * eosq;
            g300 = 1.0 + eosq * (-6.0 + 6.60937 * eosq);
            f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
            f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
            f330 = 1.0 + cosim;
            f330 = 1.875 * f330 * f330 * f330;
            del1 = 3.0 * pow(s->xnodp, 2) * ainv2;
            del2 = 2.0 * del1 * f220 * g200 * q22;
            del3 = 3.0 * del1 * f330 * g300 * q33 * aonv;
            del1 = del1 * f311 * g310 * q31 * aonv;
        s->xfact = s->xmdot + (s->omgdot + s->xnodot) - rptim + s->dmdt
                 + s->domdt + s->dnodt - s->xnodp;
        s->xlamo = fmod(s->mean_anomaly + s->right_asc_node
                 + s->argument_perigee - GSTo, TWOPI);

//        printf("----------------------------\n");
//        printf("g200:   %22.15lf\n", g200);
//        printf("g310:   %22.15lf\n", g310);
//        printf("g300:   %22.15lf\n", g300);
//        printf("f220:   %22.15lf\n", f220);
//        printf("f311:   %22.15lf\n", f311);
//        printf("f330:   %22.15lf\n", f330);
//        printf("del1:   %22.15lf\n", del1);
//        printf("del2:   %22.15lf\n", del2);
//        printf("del3:   %22.15lf\n", del3);
      }

      // Initialize the integrator
      s->xli    = s->xlamo;
      s->xni    = s->xnodp;
      s->atime  = 0.0;

//      printf("xfact:  %22.15lf\n", s->xfact);
//      printf("xlamo:  %22.15lf\n", s->xlamo);
//      printf("xli:    %22.15lf\n", s->xli);
//      printf("xni:    %22.15lf\n", s->xni);
//      printf("atime:  %22.15lf\n", s->atime);
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
    double C12   = pow(C1, 2);
        s->D2    = 4 * s->aodp * tsi * C12;
    double temp  = s->D2 * tsi * C1 / 3;
        s->D3    = (17 * s->aodp + sfour) * temp;
        s->D4    = 0.5 * temp * s->aodp * tsi * (221 * s->aodp + 31 * sfour) * C1;
        s->t3cof = s->D2 + 2 * C12;
        s->t4cof = 0.25 * (3 * s->D3 + C1 * (12 * s->D2 + 10 * C12));
        s->t5cof = 0.2  * (3 * s->D4 + 12 * C1 * s->D3 + 6 * pow(s->D2, 2)
                 + 15 * C12 * (2 * s->D2 + C12));
  }

  // Propagate at zero time since epoch here
  return sat_propagate(s, 0.0, 4, 1.0e-12, NULL, NULL);
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

  double vkmpersec     = RE * XKE / 60; // TODO: Remove?

  // Update for secular gravity and atmospheric drag
  double xmdf     = s->mean_anomaly + s->xmdot * tdelta;
  double omgadf   = s->argument_perigee + s->omgdot * tdelta;
  double xnoddf   = s->right_asc_node + s->xnodot * tdelta;
  double t2       = pow(tdelta, 2);
  double xnode    = xnoddf + s->xnodcf * t2;
  double tempa    = 1.0 - s->C1 * tdelta;
  double tempe    = s->Bstar * s->C4 * tdelta;
  double templ    = s->t2cof * t2;

  double omega    = omgadf;
  double xmp      = xmdf;

  if (s->use_simple_model == false)
  {
    double delomg = s->omgcof * tdelta;
    double delm   = s->xmcof * (pow(1.0 + s->eta * cos(xmdf), 3) - s->delmo);
    double temp   = delomg + delm;
    xmp           = xmdf + temp;
    omega         = omgadf - temp;
    double t3     = t2 * tdelta;
    double t4     = t3 * tdelta;
    tempa         = tempa - s->D2 * t2 - s->D3 * t3 - s->D4 * t4;
    tempe         = tempe + s->Bstar * s->C5 * (sin(xmp) - s->sinmo);
    templ         = templ + s->t3cof * t3 + t4 * (s->t4cof + tdelta * s->t5cof);
  }

  double nm    = s->mean_motion; // TODO: Optimize?
  double em    = s->eccentricity; // TODO: Optimize?
  double inclm = s->inclination; // TODO: Optimize?

  if (s->is_deep_space == true)
  {
    /*tc = tdelta;
*-----------------------------------------------------------------------------
*
*                           procedure dspace
*
*  this procedure provides deep space contributions to mean elements for
*    perturbing third body.  these effects have been averaged over one
*    revolution of the sun and moon.  for earth resonance effects, the
*    effects have been averaged over no revolutions of the satellite.
*    (mean motion)

  ----------------------------------------------------------------------------
    dspace
    (
        s->irez,
        s->d2201, s->d2211, s->d3210,
        s->d3222, s->d4410, s->d4422,
        s->d5220, s->d5232, s->d5421,
        s->d5433, s->dedt,  s->del1,
        s->del2,  s->del3,  s->didt,
        s->dmdt,  s->dnodt, s->domdt,
        s->argpo, s->argpdot, tdelta, tc,
        s->gsto, s->xfact, s->xlamo,
        s->no, s->atime,
        em, argpm, inclm, s->xli, mm, s->xni,
        nodem, dndt, nm, s->operationmode
    );


void dspace
     (
       int irez,
       double d2201,  double d2211,  double d3210,   double d3222,  double d4410,
       double d4422,  double d5220,  double d5232,   double d5421,  double d5433,
       double dedt,   double del1,   double del2,    double del3,   double didt,
       double dmdt,   double dnodt,  double domdt,   double argpo,  double argpdot,
       double t,      double tc,     double gsto,    double xfact,  double xlamo,
       double no,
       double* atime, double* em,    double* argpm,  double* inclm, double* xli,
       double* mm,    double* xni,   double* nodem,  double* dndt,  double* nm
     )
*/

     int iretn , iret;
     double delt, ft, theta, x2li, x2omi, xl, xldot , xnddt, xndt, xomi, g22, g32,
          g44, g52, g54, fasx2, fasx4, fasx6, rptim , step2, stepn , stepp;

     fasx2 = 0.13130908;
     fasx4 = 2.8843198;
     fasx6 = 0.37448087;
     g22   = 5.7686396;
     g32   = 0.95240898;
     g44   = 1.8014998;
     g52   = 1.0508330;
     g54   = 4.4108898;
     rptim = 4.37526908801129966e-3;
     stepp =    720.0;
     stepn =   -720.0;
     step2 = 259200.0;

     // ----------- calculate deep space resonance effects -----------
     *dndt   = 0.0;
     theta  = fmod(gsto + tc * rptim, twopi);
     *em     = *em + dedt * t;

     *inclm  = *inclm + didt * t;
     *argpm  = *argpm + domdt * t;
     *nodem  = *nodem + dnodt * t;
     *mm     = *mm + dmdt * t;

     //   sgp4fix for negative inclinations
     //   the following if statement should be commented out
     //  if (*inclm < 0.0)
     // {
     //    *inclm = -*inclm;
     //    *argpm = *argpm - pi;
     //    *nodem = *nodem + pi;
     //  }

     // - update resonances : numerical (euler-maclaurin) integration -
     // ------------------------- epoch restart ----------------------
     //   sgp4fix for propagator problems
     //   the following integration works for negative time steps and periods
     //   the specific changes are unknown because the original code was so convoluted

     // sgp4fix take out *atime = 0.0 and fix for faster operation
     ft    = 0.0;
     if (irez != 0)
       {
         // sgp4fix streamline check
         if ((*atime == 0.0) || (t * *atime <= 0.0) || (fabs(t) < fabs(*atime)) )
           {
             *atime  = 0.0;
             *xni    = no;
             *xli    = xlamo;
           }
           // sgp4fix move check outside loop
           if (t > 0.0)
               delt = stepp;
             else
               delt = stepn;

         iretn = 381; // added for do loop
         iret  =   0; // added for loop
         while (iretn == 381)
           {
             // ------------------- dot terms calculated -------------
             // ----------- near - synchronous resonance terms -------
             if (irez != 2)
               {
                 xndt  = del1 * sin(*xli - fasx2) + del2 * sin(2.0 * (*xli - fasx4)) +
                         del3 * sin(3.0 * (*xli - fasx6));
                 xldot = *xni + xfact;
                 xnddt = del1 * cos(*xli - fasx2) +
                         2.0 * del2 * cos(2.0 * (*xli - fasx4)) +
                         3.0 * del3 * cos(3.0 * (*xli - fasx6));
                 xnddt = xnddt * xldot;
               }
               else
               {
                 // --------- near - half-day resonance terms
                 xomi  = argpo + argpdot * *atime;
                 x2omi = xomi + xomi;
                 x2li  = *xli + *xli;
                 xndt  = d2201 * sin(x2omi + *xli - g22) + d2211 * sin(*xli - g22) +
                       d3210 * sin(xomi + *xli - g32)  + d3222 * sin(-xomi + *xli - g32)+
                       d4410 * sin(x2omi + x2li - g44)+ d4422 * sin(x2li - g44) +
                       d5220 * sin(xomi + *xli - g52)  + d5232 * sin(-xomi + *xli - g52)+
                       d5421 * sin(xomi + x2li - g54) + d5433 * sin(-xomi + x2li - g54);
                 xldot = *xni + xfact;
                 xnddt = d2201 * cos(x2omi + *xli - g22) + d2211 * cos(*xli - g22) +
                       d3210 * cos(xomi + *xli - g32) + d3222 * cos(-xomi + *xli - g32) +
                       d5220 * cos(xomi + *xli - g52) + d5232 * cos(-xomi + *xli - g52) +
                       2.0 * (d4410 * cos(x2omi + x2li - g44) +
                       d4422 * cos(x2li - g44) + d5421 * cos(xomi + x2li - g54) +
                       d5433 * cos(-xomi + x2li - g54));
                 xnddt = xnddt * xldot;
               }

             // ----------------------- integrator
             // sgp4fix move end checks to end of routine
             if (fabs(t - *atime) >= stepp)
               {
                 iret  = 0;
                 iretn = 381;
               }
               else // exit here
               {
                 ft    = t - *atime;
                 iretn = 0;
               }

             if (iretn == 381)
               {
                 *xli   = *xli + xldot * delt + xndt * step2;
                 *xni   = *xni + xndt * delt + xnddt * step2;
                 *atime = *atime + delt;
               }
           }  // while iretn = 381

         *nm = *xni + xndt * ft + xnddt * ft * ft * 0.5;
         xl = *xli + xldot * ft + xndt * ft * ft * 0.5;
         if (irez != 1)
           {
             *mm   = xl - 2.0 * *nodem + 2.0 * theta;
             *dndt = *nm - no;
           }
           else
           {
             *mm   = xl - *nodem - *argpm + theta;
             *dndt = *nm - no;
           }
         *nm = no + *dndt;
       }
    */
  }
/*
  if (nm <= 0.0)
  {
    return -2;
  }

  am = pow((xke / nm),x2o3) * tempa * tempa;
  nm = xke / pow(am, 1.5);
  em = em - tempe;

  if ((em >= 1.0) || (em < -1.0e-6)) // TODO: check tolerance
  {
    return -3;
  }

  // Avoid division by zero
  if ((em < 1.0e-12) && (s->operationmode != 'a'))
  {
    em  = 1.0e-6;
  }

  mm     = mm + s->no * templ;
  xlm    = mm + argpm + nodem;
  emsq   = em * em;
  temp   = 1.0 - emsq;

  nodem  = fmod(nodem, twopi);
  argpm  = fmod(argpm, twopi);
  xlm    = fmod(xlm, twopi);
  mm     = fmod(xlm - argpm - nodem, twopi);

  // ----------------- compute extra mean quantities -------------
  sinim = sin(inclm);
  cosim = cos(inclm);

  // -------------------- add lunar-solar periodics --------------
  ep     = em;
  xincp  = inclm;
  argpp  = argpm;
  nodep  = nodem;
  mp     = mm;
  sinip  = sinim;
  cosip  = cosim;

  /*
  if (s->method == 'd')
  {
    dpper
    (
        s->e3,   s->ee2,  s->peo,
        s->pgho, s->pho,  s->pinco,
        s->plo,  s->se2,  s->se3,
        s->sgh2, s->sgh3, s->sgh4,
        s->sh2,  s->sh3,  s->si2,
        s->si3,  s->sl2,  s->sl3,
        s->sl4,  tdelta,    s->xgh2,
        s->xgh3, s->xgh4, s->xh2,
        s->xh3,  s->xi2,  s->xi3,
        s->xl2,  s->xl3,  s->xl4,
        s->zmol, s->zmos, s->inclo,
        'n', ep, xincp, nodep, argpp, mp, s->operationmode
    );
    if (xincp < 0.0)
    {
      xincp  = -xincp;
      nodep = nodep + pi;
      argpp  = argpp - pi;
    }
    if ((ep < 0.0 ) || ( ep > 1.0))
    {
      //           // printf("# error ep %f\n", ep);
      s->error = 3;
      // sgp4fix add return
      return false;
    }
  } // if method = d

  // -------------------- long period periodics ------------------
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
  axnl = ep * cos(argpp);
  temp = 1.0 / (am * (1.0 - ep * ep));
  aynl = ep* sin(argpp) + temp * s->aycof;
  xl   = mp + argpp + nodep + temp * s->xlcof * axnl;

  // --------------------- solve kepler's equation ---------------
  u    = fmod(xl - nodep, twopi);
  eo1  = u;
  tem5 = 9999.9;
  ktr = 1;
  //   sgp4fix for kepler iteration
  //   the following iteration needs better limits on corrections
  while (( fabs(tem5) >= 1.0e-12) && (ktr <= 10) )
  {
    sineo1 = sin(eo1);
    coseo1 = cos(eo1);
    tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
    tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
    if(fabs(tem5) >= 0.95)
      tem5 = tem5 > 0.0 ? 0.95 : -0.95;
    eo1    = eo1 + tem5;
    ktr = ktr + 1;
  }

  // ------------- short period preliminary quantities -----------
  ecose = axnl*coseo1 + aynl*sineo1;
  esine = axnl*sineo1 - aynl*coseo1;
  el2   = axnl*axnl + aynl*aynl;
  pl    = am*(1.0-el2);
  if (pl < 0.0)
  {
    //        // printf("# error pl %f\n", pl);
    s->error = 4;
    // sgp4fix add return
    return false;
  }
  else
  {
    rl     = am * (1.0 - ecose);
    rdotl  = sqrt(am) * esine/rl;
    rvdotl = sqrt(pl) / rl;
    betal  = sqrt(1.0 - el2);
    temp   = esine / (1.0 + betal);
    sinu   = am / rl * (sineo1 - aynl - axnl * temp);
    cosu   = am / rl * (coseo1 - axnl + aynl * temp);
    su     = atan2(sinu, cosu);
    sin2u  = (cosu + cosu) * sinu;
    cos2u  = 1.0 - 2.0 * sinu * sinu;
    temp   = 1.0 / pl;
    temp1  = 0.5 * j2 * temp;
    temp2  = temp1 * temp;

    // -------------- update for short period periodics ------------
    if (s->method == 'd')
    {
      cosisq                 = cosip * cosip;
      s->con41  = 3.0*cosisq - 1.0;
      s->x1mth2 = 1.0 - cosisq;
      s->x7thm1 = 7.0*cosisq - 1.0;
    }
    mrt   = rl * (1.0 - 1.5 * temp2 * betal * s->con41) +
        0.5 * temp1 * s->x1mth2 * cos2u;
    su    = su - 0.25 * temp2 * s->x7thm1 * sin2u;
    xnode = nodep + 1.5 * temp2 * cosip * sin2u;
    xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
    mvt   = rdotl - nm * temp1 * s->x1mth2 * sin2u / xke;
    rvdot = rvdotl + nm * temp1 * (s->x1mth2 * cos2u +
        1.5 * s->con41) / xke;

    // --------------------- orientation vectors -------------------
    sinsu =  sin(su);
    cossu =  cos(su);
    snod  =  sin(xnode);
    cnod  =  cos(xnode);
    sini  =  sin(xinc);
    cosi  =  cos(xinc);
    xmx   = -snod * cosi;
    xmy   =  cnod * cosi;
    ux    =  xmx * sinsu + cnod * cossu;
    uy    =  xmy * sinsu + snod * cossu;
    uz    =  sini * sinsu;
    vx    =  xmx * cossu - cnod * sinsu;
    vy    =  xmy * cossu - snod * sinsu;
    vz    =  sini * cossu;

    // --------- position and velocity (in km and km/sec) ----------
    r[0] = (mrt * ux)* radiusearthkm;
    r[1] = (mrt * uy)* radiusearthkm;
    r[2] = (mrt * uz)* radiusearthkm;
    v[0] = (mvt * ux + rvdot * vx) * vkmpersec;
    v[1] = (mvt * uy + rvdot * vy) * vkmpersec;
    v[2] = (mvt * uz + rvdot * vz) * vkmpersec;
  }  // if pl > 0

  // sgp4fix for decaying satellites
  if (mrt < 1.0)
  {
    //        // printf("# decay condition %11.6f \n",mrt);
    s->error = 6;
    return false;
  }
*/
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


/* -----------------------------------------------------------------------------
*
*                           procedure dpper
*
*  this procedure provides deep space long period periodic contributions
*    to the mean elements.  by design, these periodics are zero at epoch.
*    this used to be dscom which included initialization, but it's really a
*    recurring function.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    e3          -
*    ee2         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    se2 , se3 , sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4 -
*    t           -
*    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
*    zmol        -
*    zmos        -
*    ep          - eccentricity                           0.0 - 1.0
*    inclo       - inclination - needed for lyddane modification
*    nodep       - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  outputs       :
*    ep          - eccentricity                           0.0 - 1.0
*    inclp       - inclination
*    nodep        - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  locals        :
*    alfdp       -
*    betdp       -
*    cosip  , sinip  , cosop  , sinop  ,
*    dalf        -
*    dbet        -
*    dls         -
*    f2, f3      -
*    pe          -
*    pgh         -
*    ph          -
*    pinc        -
*    pl          -
*    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
*    sll   , sls
*    xls         -
*    xnoh        -
*    zf          -
*    zm          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

void
dpper(sat* s, double tdelta) // TODO: Rename
{
  double alfdp, betdp, cosip, cosop, dalf, dbet, dls,
  f2,    f3,    pe,    pgh,   ph,   pinc, pl ,
  sel,   ses,   sghl,  sghs,  shll, shs,  sil,
  sinip, sinop, sinzf, sis,   sll,  sls,  xls,
  xnoh,  zf,    zm,    zel,   zes,  znl,  zns;

  // Constants
  zns   = 1.19459e-5;
  zes   = 0.01675;
  znl   = 1.5835218e-4;
  zel   = 0.05490;

  // Calculate time varying periodics
  zm    = s->zmos + zns * tdelta;

  zf    = zm + 2.0 * zes * sin(zm);
  sinzf = sin(zf);
  f2    =  0.5 * sinzf * sinzf - 0.25;
  f3    = -0.5 * sinzf * cos(zf);
  ses   = s->se2* f2 + s->se3 * f3;
  sis   = s->si2 * f2 + s->si3 * f3;
  sls   = s->sl2 * f2 + s->sl3 * f3 + s->sl4 * sinzf;
  sghs  = s->sgh2 * f2 + s->sgh3 * f3 + s->sgh4 * sinzf;
  shs   = s->sh2 * f2 + s->sh3 * f3;
  zm    = s->zmol + znl * tdelta;

  zf    = zm + 2.0 * zel * sin(zm);
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

  // Apply periodics directly using Lyddane choice
  if (s->inclination_lp >= 0.2)
  {
    ph  = ph / sinip;
    pgh = pgh - cosip * ph;

    s->argument_perigee_lp += pgh;
    s->right_asc_node_lp   += ph;
    s->mean_anomaly_lp     += pl;
  }
  else
  {
    /* ---- apply periodics with lyddane modification ---- */
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
    //  sgp4fix for afspc written intrinsic functions
    // nodep used without a trigonometric function ahead
    if (fabs(xnoh - s->right_asc_node_lp) > PI)
      if (s->right_asc_node_lp < xnoh)
        s->right_asc_node_lp = s->right_asc_node_lp + TWOPI;
      else
        s->right_asc_node_lp = s->right_asc_node_lp - TWOPI;
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
