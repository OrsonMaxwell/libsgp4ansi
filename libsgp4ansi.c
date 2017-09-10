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
      dsinit
    sgp4(0)

sgp4()
  if deepspace
    dspace()
    dpper()          <------
*/

// SGP4/SDP4 propagation function implementation
int
orbit_sgp4(sat*, double, unsigned int, double, vec3*, vec3*);

// Deep space long period periodics contributions to the mean elements
void
orbit_dslongper(sat*, double, double*, double*, double*, double*, double*);


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

/*
 * SGP4 propagation function implementation
 *
 * Inputs:  sat       - orbit struct pointer with initialized orbital data
 *          tdelta    - Time since epoch, minutes
 *          maxiter   - Kepler's equation maximum iteration count
 *          tolerance - Kepler's equation desired precision tolerance
 * Outputs: pos       - 3D position vector in TEME frame in km
 *          vel       - 3D velocity vector in TEME frame in km/sec
 * Returns: 0         - Success
 *          1         - Mean motion zero or less
 *          2         - Nonsensical orbital eccentricity (e >= 1; e < -0.00001)
 *          3         - Long periodics result error
 *          4         - Short period preliminary quantities error
 *          5         - Decaying satellite
 */
int
orbit_sgp4
(
    sat* sat,
    double tdelta,
    unsigned int maxiter,
    double tolerance,
    vec3* pos,
    vec3* vel
)
{
  double r[3], v[3];
  sgp4(1, sat, tdelta, r, v);
  pos->x = r[0];
  pos->y = r[1];
  pos->z = r[2];
  vel->x = v[0];
  vel->y = v[1];
  vel->z = v[2];
  return 0;
}

/*
 * Deep space long period periodics contributions to the mean elements
 *
 * Inputs:  sat    - orbit struct pointer with full orbital data
 *          tdelta - Time since epoch, minutes
 * Outputs: sat    - orbit struct pointer with long period periodics applied
 *          e      - Eccentricity with lpp applied
 *          i      - Inclination with lpp applied
 *          alpha  - Right ascension of ascension node with lpp applied
 *          omega  - Argument of perigee with lpp applied
 *          mp     - Mean motion term with lpp applied
 * Returns: None
 *
 * If any of the output pointers are NULL, assume initial run and apply lpp to
 * satellite data directly
 */
void
orbit_dslongper
(
    sat* sat,
    double tdelta,
    double* e,
    double* i,
    double* alpha,
    double* omega,
    double* mp
)
{

}

void print_orbit(sat* sat, char* caption)
{
  FILE *f = fopen("orbit.txt", "a");
  fprintf(f, "===== %s =====\n", caption);
  fprintf(f, "----- TLE portion -----\n");
  fprintf(f, "number:\t\t%d\n", sat->norad_number);
  fprintf(f, "nprimediv2:\t%20.10lf\n", sat->mean_motion_dt2);
  fprintf(f, "ndprimediv6:\t%20.10lf\n", sat->mean_motion_ddt6);
  fprintf(f, "Bstar:\t\t%20.10lf\n", sat->Bstar);
  fprintf(f, "i:\t\t%20.10lf\n", sat->inclination);
  fprintf(f, "alpha:\t\t%20.10lf\n", sat->right_asc_node);
  fprintf(f, "e:\t\t%20.10lf\n", sat->eccentricity);
  fprintf(f, "omega:\t\t%20.10lf\n", sat->argument_perigee);
  fprintf(f, "Mo:\t\t%20.10lf\n", sat->mean_anomaly);
  fprintf(f, "no:\t\t%20.10lf\n", sat->mean_motion);
  fprintf(f, "----- Time -----\n");
  fprintf(f, "epoch:\t\t%s.%lf\n", ctime(&sat->epoch), sat->epoch_ms);
  fprintf(f, "julepoch:\t%20.10lf\n", sat->julepoch);
  fprintf(f, "GSTo:\t\t%20.10lf\n", sat->GSTo);
  fprintf(f, "----- Flags -----\n");
  fprintf(f, "Deepsp:\t\t%d\n", sat->is_deep_space);
  fprintf(f, "Loworb:\t\t%d\n", sat->use_simple_model);
  fprintf(f, "----- Standard terms -----\n");
  fprintf(f, "a:\t\t%20.10lf\n", sat->a);
  fprintf(f, "altapoR:\t\t%20.10lf\n", sat->altapoR);
  fprintf(f, "altperR:\t\t%20.10lf\n", sat->altperR);
  fprintf(f, "aycof:\t\t%20.10lf\n", sat->aycof);
  fprintf(f, "C1:\t\t%20.10lf\n", sat->C1);
  fprintf(f, "C4:\t\t%20.10lf\n", sat->C4);
  fprintf(f, "C5:\t\t%20.10lf\n", sat->C5);
  fprintf(f, "con41:\t\t%20.10lf\n", sat->x3thm1);
  fprintf(f, "d2:\t\t%20.10lf\n", sat->d2);
  fprintf(f, "d3:\t\t%20.10lf\n", sat->d3);
  fprintf(f, "d4:\t\t%20.10lf\n", sat->d4);
  fprintf(f, "delMo:\t\t%20.10lf\n", sat->delMo);
  fprintf(f, "eta:\t\t%20.10lf\n", sat->eta);
  fprintf(f, "mdot:\t\t%20.10lf\n", sat->mdot);
  fprintf(f, "nodecf:\t\t%20.10lf\n", sat->nodecf);
  fprintf(f, "nodedot:\t%20.10lf\n", sat->nodedot);
  fprintf(f, "omegaprime:\t%20.10lf\n", sat->omegaprime);
  fprintf(f, "omgcof:\t\t%20.10lf\n", sat->omgcof);
  fprintf(f, "sinMo:\t\t%20.10lf\n", sat->sinMo);
  fprintf(f, "t2cof:\t\t%20.10lf\n", sat->t2cof);
  fprintf(f, "t3cof:\t\t%20.10lf\n", sat->t3cof);
  fprintf(f, "t4cof:\t\t%20.10lf\n", sat->t4cof);
  fprintf(f, "t5cof:\t\t%20.10lf\n", sat->t5cof);
  fprintf(f, "x1mth2:\t\t%20.10lf\n", sat->x1mth2);
  fprintf(f, "x7thm1:\t\t%20.10lf\n", sat->x7thm1);
  fprintf(f, "xlcof:\t\t%20.10lf\n", sat->xlcof);
  fprintf(f, "xmcof:\t\t%20.10lf\n", sat->xmcof);
  fprintf(f, "----- Deepspace terms -----\n");
  fprintf(f, "e3:\t\t%20.10lf\n", sat->e3);
  fprintf(f, "ee2:\t\t%20.10lf\n", sat->ee2);
  fprintf(f, "peo:\t\t%20.10lf\n", sat->peo);
  fprintf(f, "pgho:\t\t%20.10lf\n", sat->pgho);
  fprintf(f, "pho:\t\t%20.10lf\n", sat->pho);
  fprintf(f, "pinco:\t\t%20.10lf\n", sat->pinco);
  fprintf(f, "plo:\t\t%20.10lf\n", sat->plo);
  fprintf(f, "se2:\t\t%20.10lf\n", sat->se2);
  fprintf(f, "se3:\t\t%20.10lf\n", sat->se3);
  fprintf(f, "sgh2:\t\t%20.10lf\n", sat->sgh2);
  fprintf(f, "sgh3:\t\t%20.10lf\n", sat->sgh3);
  fprintf(f, "sgh4:\t\t%20.10lf\n", sat->sgh4);
  fprintf(f, "sh2:\t\t%20.10lf\n", sat->sh2);
  fprintf(f, "sh3:\t\t%20.10lf\n", sat->sh3);
  fprintf(f, "si2:\t\t%20.10lf\n", sat->si2);
  fprintf(f, "si3:\t\t%20.10lf\n", sat->si3);
  fprintf(f, "sl2:\t\t%20.10lf\n", sat->sl2);
  fprintf(f, "sl3:\t\t%20.10lf\n", sat->sl3);
  fprintf(f, "sl4:\t\t%20.10lf\n", sat->sl4);
  fprintf(f, "xgh2:\t\t%20.10lf\n", sat->xgh2);
  fprintf(f, "xgh3:\t\t%20.10lf\n", sat->xgh3);
  fprintf(f, "xgh4:\t\t%20.10lf\n", sat->xgh4);
  fprintf(f, "xh2:\t\t%20.10lf\n", sat->xh2);
  fprintf(f, "xh3:\t\t%20.10lf\n", sat->xh3);
  fprintf(f, "xi2:\t\t%20.10lf\n", sat->xi2);
  fprintf(f, "xi3:\t\t%20.10lf\n", sat->xi3);
  fprintf(f, "xl2:\t\t%20.10lf\n", sat->xl2);
  fprintf(f, "xl3:\t\t%20.10lf\n", sat->xl3);
  fprintf(f, "xl4:\t\t%20.10lf\n", sat->xl4);
  fprintf(f, "zmol:\t\t%20.10lf\n", sat->zmol);
  fprintf(f, "zmos:\t\t%20.10lf\n", sat->zmos);
  fprintf(f, "----- Resonant terms -----\n");
  fprintf(f, "d2201:\t\t%20.10lf\n", sat->d2201);
  fprintf(f, "d2211:\t\t%20.10lf\n", sat->d2211);
  fprintf(f, "d3210:\t\t%20.10lf\n", sat->d3210);
  fprintf(f, "d3222:\t\t%20.10lf\n", sat->d3222);
  fprintf(f, "d4410:\t\t%20.10lf\n", sat->d4410);
  fprintf(f, "d4422:\t\t%20.10lf\n", sat->d4422);
  fprintf(f, "d5220:\t\t%20.10lf\n", sat->d5220);
  fprintf(f, "d5232:\t\t%20.10lf\n", sat->d5232);
  fprintf(f, "d5421:\t\t%20.10lf\n", sat->d5421);
  fprintf(f, "d5433:\t\t%20.10lf\n", sat->d5433);
  fprintf(f, "dedt:\t\t%20.10lf\n", sat->dedt);
  fprintf(f, "didt:\t\t%20.10lf\n", sat->didt);
  fprintf(f, "dmdt:\t\t%20.10lf\n", sat->dmdt);
  fprintf(f, "dnodt:\t\t%20.10lf\n", sat->dnodt);
  fprintf(f, "domdt:\t\t%20.10lf\n", sat->domdt);
  fprintf(f, "del1:\t\t%20.10lf\n", sat->del1);
  fprintf(f, "del2:\t\t%20.10lf\n", sat->del2);
  fprintf(f, "del3:\t\t%20.10lf\n", sat->del3);
  fprintf(f, "xfact:\t\t%20.10lf\n", sat->xfact);
  fprintf(f, "xlamo:\t\t%20.10lf\n", sat->xlamo);
  fprintf(f, "xli:\t\t%20.10lf\n", sat->xli);
  fprintf(f, "xni:\t\t%20.10lf\n", sat->xni);
  fprintf(f, "----- Legacy -----\n");
  fprintf(f, "isimp:\t\t%d\n", sat->isimp);
  fprintf(f, "error:\t\t%d\n", sat->error);
  fprintf(f, "method:\t\t%c\n", sat->method);
  fprintf(f, "operationmode:\t\t%c\n", sat->operationmode);
  fprintf(f, "cc1:\t\t%20.10lf\n", sat->cc1);
  fprintf(f, "cc4:\t\t%20.10lf\n", sat->cc4);
  fprintf(f, "cc5:\t\t%20.10lf\n", sat->cc5);
  fprintf(f, "delmo:\t\t%20.10lf\n", sat->delmo);
  fprintf(f, "argpdot:\t\t%20.10lf\n", sat->argpdot);
  fprintf(f, "sinmao:\t\t%20.10lf\n", sat->sinmao);
  fprintf(f, "t:\t\t%20.10lf\n", sat->t);
  fprintf(f, "irez:\t\t%d\n", sat->irez);
  fprintf(f, "gsto:\t\t%20.10lf\n", sat->gsto);
  fprintf(f, "atime:\t\t%20.10lf\n", sat->atime);
  fprintf(f, "bstar:\t\t%20.10lf\n", sat->bstar);
  fprintf(f, "ecco:\t\t%20.10lf\n", sat->ecco);
  fprintf(f, "argpo:\t\t%20.10lf\n", sat->argpo);
  fprintf(f, "inclo:\t\t%20.10lf\n", sat->inclo);
  fprintf(f, "mo:\t\t%20.10lf\n", sat->mo);
  fprintf(f, "nodeo:\t\t%20.10lf\n", sat->nodeo);
  fprintf(f, "alta:\t\t%20.10lf\n", sat->alta);
  fprintf(f, "altp:\t\t%20.10lf\n", sat->altp);
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
tle2orbit(char* tlestr0, char* tlestr1, char* tlestr2, sat* s)
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

  // Make this return milliseconds
  if (fractday2unix(epochyr, epochdays, &s->epoch, &s->epoch_ms) != 0)
  {
    return -1;
  }

  s->julepoch = unix2jul(&s->epoch, s->epoch_ms);

  s->mean_motion_ddt6 = nddot * pow(10, nexp);
  s->Bstar = Bstar * pow(10, Bexp);

  // TODO: Convert to SGP4 units and expand

  return orbit_init(s);
}

/*
 * Expand SGP4/SDP4 orbit elements from an orbit containing NORAD TLE portion
 *
 * Inputs:  sat  - orbit struct pointer to an empty orbit
 * Outputs: sat  - orbit struct pointer with full orbital data
 * Returns: None TODO: add return values
 */
int
orbit_init(sat* sat)
{

  sgp4init(0, 'i', 0101, sat->julepoch-2433281.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, sat);
  return 0;
}

/*
 * Get position and velocity vectors in the TEME frame at given time since epoch
 *
 * Inputs:  sat       - orbit struct pointer with initialized orbital data
 *          tdelta    - Time since orbit TLE epoch, minutes
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
orbit_prop
(
    sat* sat,
    double tdelta,
    unsigned int maxiter,
    double tolerance,
    vec3* pos,
    vec3* vel
)
{
  if ((sat == NULL) ||
      (maxiter < 1) ||
      (tolerance <= 0.0) ||
      (pos == NULL) ||
      (vel == NULL))
  {
    return -1;
  }

//  double tdelta = difftime(*time + msec / 1000,
//                           sat->epoch + sat->epoch_ms / 1000) / 60.0L;

  return orbit_sgp4(sat, tdelta, maxiter, tolerance, pos, vel);
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
orbit_at
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

  return orbit_prop(sat, tdelta, maxiter, tolerance, pos, vel);
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

void dpper
     (
       double e3,     double ee2,    double peo,     double pgho,   double pho,
       double pinco,  double plo,    double se2,     double se3,    double sgh2,
       double sgh3,   double sgh4,   double sh2,     double sh3,    double si2,
       double si3,    double sl2,    double sl3,     double sl4,    double t,
       double xgh2,   double xgh3,   double xgh4,    double xh2,    double xh3,
       double xi2,    double xi3,    double xl2,     double xl3,    double xl4,
       double zmol,   double zmos,   double inclo,
       char init,
       double* ep,    double* inclp, double* nodep,  double* argpp, double* mp,
       char opsmode
     )
{
     /* --------------------- local variables ------------------------ */
     const double twopi = 2.0 * PI;
     double alfdp, betdp, cosip, cosop, dalf, dbet, dls,
          f2,    f3,    pe,    pgh,   ph,   pinc, pl ,
          sel,   ses,   sghl,  sghs,  shll, shs,  sil,
          sinip, sinop, sinzf, sis,   sll,  sls,  xls,
          xnoh,  zf,    zm,    zel,   zes,  znl,  zns;

     /* ---------------------- constants ----------------------------- */
     zns   = 1.19459e-5;
     zes   = 0.01675;
     znl   = 1.5835218e-4;
     zel   = 0.05490;

     /* --------------- calculate time varying periodics ----------- */
     zm    = zmos + zns * t;
     // be sure that the initial call has time set to zero
     if (init == 'y')
         zm = zmos;
     zf    = zm + 2.0 * zes * sin(zm);
     sinzf = sin(zf);
     f2    =  0.5 * sinzf * sinzf - 0.25;
     f3    = -0.5 * sinzf * cos(zf);
     ses   = se2* f2 + se3 * f3;
     sis   = si2 * f2 + si3 * f3;
     sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
     sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
     shs   = sh2 * f2 + sh3 * f3;
     zm    = zmol + znl * t;
     if (init == 'y')
         zm = zmol;
     zf    = zm + 2.0 * zel * sin(zm);
     sinzf = sin(zf);
     f2    =  0.5 * sinzf * sinzf - 0.25;
     f3    = -0.5 * sinzf * cos(zf);
     sel   = ee2 * f2 + e3 * f3;
     sil   = xi2 * f2 + xi3 * f3;
     sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
     sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
     shll  = xh2 * f2 + xh3 * f3;
     pe    = ses + sel;
     pinc  = sis + sil;
     pl    = sls + sll;
     pgh   = sghs + sghl;
     ph    = shs + shll;

     if (init == 'n')
       {
       pe    = pe - peo;
       pinc  = pinc - pinco;
       pl    = pl - plo;
       pgh   = pgh - pgho;
       ph    = ph - pho;
       *inclp = *inclp + pinc;
       *ep    = *ep + pe;
       sinip = sin(*inclp);
       cosip = cos(*inclp);

       /* ----------------- apply periodics directly ------------ */
       //  sgp4fix for lyddane choice
       //  strn3 used original inclination - this is technically feasible
       //  gsfc used perturbed inclination - also technically feasible
       //  probably best to readjust the 0.2 limit value and limit discontinuity
       //  0.2 rad = 11.45916 deg
       //  use next line for original strn3 approach and original inclination
       //  if (inclo >= 0.2)
       //  use next line for gsfc version and perturbed inclination
       if (*inclp >= 0.2)
         {
           ph     = ph / sinip;
           pgh    = pgh - cosip * ph;
           *argpp  = *argpp + pgh;
           *nodep  = *nodep + ph;
           *mp     = *mp + pl;
         }
         else
         {
           /* ---- apply periodics with lyddane modification ---- */
           sinop  = sin(*nodep);
           cosop  = cos(*nodep);
           alfdp  = sinip * sinop;
           betdp  = sinip * cosop;
           dalf   =  ph * cosop + pinc * cosip * sinop;
           dbet   = -ph * sinop + pinc * cosip * cosop;
           alfdp  = alfdp + dalf;
           betdp  = betdp + dbet;
           *nodep  = fmod(*nodep, twopi);
           //  sgp4fix for afspc written intrinsic functions
           // nodep used without a trigonometric function ahead
           if ((*nodep < 0.0) && (opsmode == 'a'))
               *nodep = *nodep + twopi;
           xls    = *mp + *argpp + cosip * *nodep;
           dls    = pl + pgh - pinc * *nodep * sinip;
           xls    = xls + dls;
           xnoh   = *nodep;
           *nodep  = atan2(alfdp, betdp);
           //  sgp4fix for afspc written intrinsic functions
           // nodep used without a trigonometric function ahead
           if ((*nodep < 0.0) && (opsmode == 'a'))
               *nodep = *nodep + twopi;
           if (fabs(xnoh - *nodep) > PI)
             if (*nodep < xnoh)
                *nodep = *nodep + twopi;
               else
                *nodep = *nodep - twopi;
           *mp    = *mp + pl;
           *argpp = xls - *mp - cosip * *nodep;
         }
       }   // if init == 'n'

//#include "debug1.cpp"
}  // end dpper

/*-----------------------------------------------------------------------------
*
*                           procedure dscom
*
*  this procedure provides deep space common items used by both the secular
*    and periodics subroutines.  input is provided as shown. this routine
*    used to be called dpper, but the functions inside weren't well organized.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    epoch       -
*    ep          - eccentricity
*    argpp       - argument of perigee
*    tc          -
*    inclp       - inclination
*    nodep       - right ascension of ascending node
*    np          - mean motion
*
*  outputs       :
*    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
*    day         -
*    e3          -
*    ee2         -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    gam         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    rtemsq      -
*    se2, se3         -
*    sgh2, sgh3, sgh4        -
*    sh2, sh3, si2, si3, sl2, sl3, sl4         -
*    s1, s2, s3, s4, s5, s6, s7          -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
*    nm          - mean motion
*    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*    zmol        -
*    zmos        -
*
*  locals        :
*    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
*    betasq      -
*    cc          -
*    ctem, stem        -
*    x1, x2, x3, x4, x5, x6, x7, x8          -
*    xnodce      -
*    xnoi        -
*    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
*    zcosi  , zsini  , zcosil , zsinil ,
*    zx          -
*    zy          -
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

void dscom
     (
       double epoch,  double ep,     double argpp,   double tc,     double inclp,
       double nodep,  double np,
       double* snodm, double* cnodm, double* sinim,  double* cosim, double* sinomm,
       double* cosomm,double* day,   double* e3,     double* ee2,   double* em,
       double* emsq,  double* gam,   double* peo,    double* pgho,  double* pho,
       double* pinco, double* plo,   double* rtemsq, double* se2,   double* se3,
       double* sgh2,  double* sgh3,  double* sgh4,   double* sh2,   double* sh3,
       double* si2,   double* si3,   double* sl2,    double* sl3,   double* sl4,
       double* s1,    double* s2,    double* s3,     double* s4,    double* s5,
       double* s6,    double* s7,    double* ss1,    double* ss2,   double* ss3,
       double* ss4,   double* ss5,   double* ss6,    double* ss7,   double* sz1,
       double* sz2,   double* sz3,   double* sz11,   double* sz12,  double* sz13,
       double* sz21,  double* sz22,  double* sz23,   double* sz31,  double* sz32,
       double* sz33,  double* xgh2,  double* xgh3,   double* xgh4,  double* xh2,
       double* xh3,   double* xi2,   double* xi3,    double* xl2,   double* xl3,
       double* xl4,   double* nm,    double* z1,     double* z2,    double* z3,
       double* z11,   double* z12,   double* z13,    double* z21,   double* z22,
       double* z23,   double* z31,   double* z32,    double* z33,   double* zmol,
       double* zmos
     )
{
     /* -------------------------- constants ------------------------- */
     const double zes     =  0.01675;
     const double zel     =  0.05490;
     const double c1ss    =  2.9864797e-6;
     const double c1l     =  4.7968065e-7;
     const double zsinis  =  0.39785416;
     const double zcosis  =  0.91744867;
     const double zcosgs  =  0.1945905;
     const double zsings  = -0.98088458;
     const double twopi   =  2.0 * PI;

     /* --------------------- local variables ------------------------ */
     int lsflg;
     double a1    , a2    , a3    , a4    , a5    , a6    , a7    ,
        a8    , a9    , a10   , betasq, cc    , ctem  , stem  ,
        x1    , x2    , x3    , x4    , x5    , x6    , x7    ,
        x8    , xnodce, xnoi  , zcosg , zcosgl, zcosh , zcoshl,
        zcosi , zcosil, zsing , zsingl, zsinh , zsinhl, zsini ,
        zsinil, zx    , zy;

     *nm     = np;
     *em     = ep;
     *snodm  = sin(nodep);
     *cnodm  = cos(nodep);
     *sinomm = sin(argpp);
     *cosomm = cos(argpp);
     *sinim  = sin(inclp);
     *cosim  = cos(inclp);
     *emsq   = *em * *em;
     betasq = 1.0 - *emsq;
     *rtemsq = sqrt(betasq);

     /* ----------------- initialize lunar solar terms --------------- */
     *peo    = 0.0;
     *pinco  = 0.0;
     *plo    = 0.0;
     *pgho   = 0.0;
     *pho    = 0.0;
     *day    = epoch + 18261.5 + tc / 1440.0;
     xnodce = fmod(4.5236020 - 9.2422029e-4 * *day, twopi);
     stem   = sin(xnodce);
     ctem   = cos(xnodce);
     zcosil = 0.91375164 - 0.03568096 * ctem;
     zsinil = sqrt(1.0 - zcosil * zcosil);
     zsinhl = 0.089683511 * stem / zsinil;
     zcoshl = sqrt(1.0 - zsinhl * zsinhl);
     *gam    = 5.8351514 + 0.0019443680 * *day;
     zx     = 0.39785416 * stem / zsinil;
     zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
     zx     = atan2(zx, zy);
     zx     = *gam + zx - xnodce;
     zcosgl = cos(zx);
     zsingl = sin(zx);

     /* ------------------------- do solar terms --------------------- */
     zcosg = zcosgs;
     zsing = zsings;
     zcosi = zcosis;
     zsini = zsinis;
     zcosh = *cnodm;
     zsinh = *snodm;
     cc    = c1ss;
     xnoi  = 1.0 / *nm;

     for (lsflg = 1; lsflg <= 2; lsflg++)
       {
         a1  =   zcosg * zcosh + zsing * zcosi * zsinh;
         a3  =  -zsing * zcosh + zcosg * zcosi * zsinh;
         a7  =  -zcosg * zsinh + zsing * zcosi * zcosh;
         a8  =   zsing * zsini;
         a9  =   zsing * zsinh + zcosg * zcosi * zcosh;
         a10 =   zcosg * zsini;
         a2  =   *cosim * a7 + *sinim * a8;
         a4  =   *cosim * a9 + *sinim * a10;
         a5  =  -*sinim * a7 + *cosim * a8;
         a6  =  -*sinim * a9 + *cosim * a10;

         x1  =  a1 * *cosomm + a2 * *sinomm;
         x2  =  a3 * *cosomm + a4 * *sinomm;
         x3  = -a1 * *sinomm + a2 * *cosomm;
         x4  = -a3 * *sinomm + a4 * *cosomm;
         x5  =  a5 * *sinomm;
         x6  =  a6 * *sinomm;
         x7  =  a5 * *cosomm;
         x8  =  a6 * *cosomm;

         *z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
         *z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
         *z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
         *z1  =  3.0 *  (a1 * a1 + a2 * a2) + *z31 * *emsq;
         *z2  =  6.0 *  (a1 * a3 + a2 * a4) + *z32 * *emsq;
         *z3  =  3.0 *  (a3 * a3 + a4 * a4) + *z33 * *emsq;
         *z11 = -6.0 * a1 * a5 + *emsq *  (-24.0 * x1 * x7-6.0 * x3 * x5);
         *z12 = -6.0 *  (a1 * a6 + a3 * a5) + *emsq *
                (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
         *z13 = -6.0 * a3 * a6 + *emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
         *z21 =  6.0 * a2 * a5 + *emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
         *z22 =  6.0 *  (a4 * a5 + a2 * a6) + *emsq *
                (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
         *z23 =  6.0 * a4 * a6 + *emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
         *z1  = *z1 + *z1 + betasq * *z31;
         *z2  = *z2 + *z2 + betasq * *z32;
         *z3  = *z3 + *z3 + betasq * *z33;
         *s3  = cc * xnoi;
         *s2  = -0.5 * *s3 / *rtemsq;
         *s4  = *s3 * *rtemsq;
         *s1  = -15.0 * *em * *s4;
         *s5  = x1 * x3 + x2 * x4;
         *s6  = x2 * x3 + x1 * x4;
         *s7  = x2 * x4 - x1 * x3;

         /* ----------------------- do lunar terms ------------------- */
         if (lsflg == 1)
           {
             *ss1   = *s1;
             *ss2   = *s2;
             *ss3   = *s3;
             *ss4   = *s4;
             *ss5   = *s5;
             *ss6   = *s6;
             *ss7   = *s7;
             *sz1   = *z1;
             *sz2   = *z2;
             *sz3   = *z3;
             *sz11  = *z11;
             *sz12  = *z12;
             *sz13  = *z13;
             *sz21  = *z21;
             *sz22  = *z22;
             *sz23  = *z23;
             *sz31  = *z31;
             *sz32  = *z32;
             *sz33  = *z33;
             zcosg = zcosgl;
             zsing = zsingl;
             zcosi = zcosil;
             zsini = zsinil;
             zcosh = zcoshl * *cnodm + zsinhl * *snodm;
             zsinh = *snodm * zcoshl - *cnodm * zsinhl;
             cc    = c1l;
          }
       }

     *zmol = fmod(4.7199672 + 0.22997150  * *day - *gam, twopi);
     *zmos = fmod(6.2565837 + 0.017201977 * *day, twopi);

     /* ------------------------ do solar terms ---------------------- */
     *se2  =   2.0 * *ss1 * *ss6;
     *se3  =   2.0 * *ss1 * *ss7;
     *si2  =   2.0 * *ss2 * *sz12;
     *si3  =   2.0 * *ss2 * (*sz13 - *sz11);
     *sl2  =  -2.0 * *ss3 * *sz2;
     *sl3  =  -2.0 * *ss3 * (*sz3 - *sz1);
     *sl4  =  -2.0 * *ss3 * (-21.0 - 9.0 * *emsq) * zes;
     *sgh2 =   2.0 * *ss4 * *sz32;
     *sgh3 =   2.0 * *ss4 * (*sz33 - *sz31);
     *sgh4 = -18.0 * *ss4 * zes;
     *sh2  =  -2.0 * *ss2 * *sz22;
     *sh3  =  -2.0 * *ss2 * (*sz23 - *sz21);

     /* ------------------------ do lunar terms ---------------------- */
     *ee2  =   2.0 * *s1 * *s6;
     *e3   =   2.0 * *s1 * *s7;
     *xi2  =   2.0 * *s2 * *z12;
     *xi3  =   2.0 * *s2 * (*z13 - *z11);
     *xl2  =  -2.0 * *s3 * *z2;
     *xl3  =  -2.0 * *s3 * (*z3 - *z1);
     *xl4  =  -2.0 * *s3 * (-21.0 - 9.0 * *emsq) * zel;
     *xgh2 =   2.0 * *s4 * *z32;
     *xgh3 =   2.0 * *s4 * (*z33 - *z31);
     *xgh4 = -18.0 * *s4 * zel;
     *xh2  =  -2.0 * *s2 * *z22;
     *xh3  =  -2.0 * *s2 * (*z23 - *z21);

//#include "debug2.cpp"
}  // end dscom

/*-----------------------------------------------------------------------------
*
*                           procedure dsinit
*
*  this procedure provides deep space contributions to mean motion dot due
*    to geopotential resonance with half day and one day orbits.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    cosim, sinim-
*    emsq        - eccentricity squared
*    argpo       - argument of perigee
*    s1, s2, s3, s4, s5      -
*    ss1, ss2, ss3, ss4, ss5 -
*    sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
*    t           - time
*    tc          -
*    gsto        - greenwich sidereal time                   rad
*    mo          - mean anomaly
*    mdot        - mean anomaly dot (rate)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*    nodedot     - right ascension of ascending node dot (rate)
*    xpidot      -
*    z1, z3, z11, z13, z21, z23, z31, z33 -
*    eccm        - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    xn          - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right ascension of ascending node
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    atime       -
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433    -
*    dedt        -
*    didt        -
*    dmdt        -
*    dndt        -
*    dnodt       -
*    domdt       -
*    del1, del2, del3        -
*    ses  , sghl , sghs , sgs  , shl  , shs  , sis  , sls
*    theta       -
*    xfact       -
*    xlamo       -
*    xli         -
*    xni
*
*  locals        :
*    ainv2       -
*    aonv        -
*    cosisq      -
*    eoc         -
*    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543  -
*    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533  -
*    sini2       -
*    temp        -
*    temp1       -
*    theta       -
*    xno2        -
*
*  coupling      :
*    getgravconst
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

void dsinit
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
     )
{
     /* --------------------- local variables ------------------------ */
     const double twopi = 2.0 * PI;

     double ainv2 , aonv=0.0, cosisq, eoc, f220 , f221  , f311  ,
          f321  , f322  , f330  , f441  , f442  , f522  , f523  ,
          f542  , f543  , g200  , g201  , g211  , g300  , g310  ,
          g322  , g410  , g422  , g520  , g521  , g532  , g533  ,
          ses   , sgs   , sghl  , sghs  , shs   , shll  , sis   ,
          sini2 , sls   , temp  , temp1 , theta , xno2  , q22   ,
          q31   , q33   , root22, root44, root54, rptim , root32,
          root52, x2o3  , xke   , znl   , emo   , zns   , emsqo,
          tumin, mu, radiusearthkm, j2, j3, j4, j3oj2;

     q22    = 1.7891679e-6;
     q31    = 2.1460748e-6;
     q33    = 2.2123015e-7;
     root22 = 1.7891679e-6;
     root44 = 7.3636953e-9;
     root54 = 2.1765803e-9;
     rptim  = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
     root32 = 3.7393792e-7;
     root52 = 1.1428639e-7;
     x2o3   = 2.0 / 3.0;
     znl    = 1.5835218e-4;
     zns    = 1.19459e-5;

     // sgp4fix identify constants and allow alternate values
     // TODO: replace by macros
     tumin = TUMIN;
     mu = GM;
     radiusearthkm = RE;
     xke = XKE;
     j2 = J2;
     j3 = J3;
     j4 = J4;
     j3oj2 = J3DIVJ2;

     /* -------------------- deep space initialization ------------ */
     *irez = 0;
     if ((*nm < 0.0052359877) && (*nm > 0.0034906585))
         *irez = 1;
     if ((*nm >= 8.26e-3) && (*nm <= 9.24e-3) && (*em >= 0.5))
         *irez = 2;

     /* ------------------------ do solar terms ------------------- */
     ses  =  ss1 * zns * ss5;
     sis  =  ss2 * zns * (sz11 + sz13);
     sls  = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
     sghs =  ss4 * zns * (sz31 + sz33 - 6.0);
     shs  = -zns * ss2 * (sz21 + sz23);
     // sgp4fix for 180 deg incl
     if ((*inclm < 5.2359877e-2) || (*inclm > PI - 5.2359877e-2))
       shs = 0.0;
     if (sinim != 0.0)
       shs = shs / sinim;
     sgs  = sghs - cosim * shs;

     /* ------------------------- do lunar terms ------------------ */
     *dedt = ses + s1 * znl * s5;
     *didt = sis + s2 * znl * (z11 + z13);
     *dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
     sghl = s4 * znl * (z31 + z33 - 6.0);
     shll = -znl * s2 * (z21 + z23);
     // sgp4fix for 180 deg incl
     if ((*inclm < 5.2359877e-2) || (*inclm > PI - 5.2359877e-2))
         shll = 0.0;
     *domdt = sgs + sghl;
     *dnodt = shs;
     if (sinim != 0.0)
       {
         *domdt = *domdt - cosim / sinim * shll;
         *dnodt = *dnodt + shll / sinim;
       }

     /* ----------- calculate deep space resonance effects -------- */
     *dndt   = 0.0;
     theta  = fmod(gsto + tc * rptim, twopi);
     *em     = *em + *dedt * t;
     *inclm  = *inclm + *didt * t;
     *argpm  = *argpm + *domdt * t;
     *nodem  = *nodem + *dnodt * t;
     *mm     = *mm + *dmdt * t;
     //   sgp4fix for negative inclinations
     //   the following if statement should be commented out
     //if (*inclm < 0.0)
     //  {
     //    *inclm  = -*inclm;
     //    *argpm  = *argpm - pi;
     //    *nodem = *nodem + pi;
     //  }

     /* -------------- initialize the resonance terms ------------- */
     if (*irez != 0)
       {
         aonv = pow(*nm / xke, x2o3);

         /* ---------- geopotential resonance for 12 hour orbits ------ */
         if (*irez == 2)
           {
             cosisq = cosim * cosim;
             emo    = *em;
             *em     = ecco;
             emsqo  = emsq;
             emsq   = eccsq;
             eoc    = *em * emsq;
             g201   = -0.306 - (*em - 0.64) * 0.440;

             if (*em <= 0.65)
               {
                 g211 =    3.616  -  13.2470 * *em +  16.2900 * emsq;
                 g310 =  -19.302  + 117.3900 * *em - 228.4190 * emsq +  156.5910 * eoc;
                 g322 =  -18.9068 + 109.7927 * *em - 214.6334 * emsq +  146.5816 * eoc;
                 g410 =  -41.122  + 242.6940 * *em - 471.0940 * emsq +  313.9530 * eoc;
                 g422 = -146.407  + 841.8800 * *em - 1629.014 * emsq + 1083.4350 * eoc;
                 g520 = -532.114  + 3017.977 * *em - 5740.032 * emsq + 3708.2760 * eoc;
               }
               else
               {
                 g211 =   -72.099 +   331.819 * *em -   508.738 * emsq +   266.724 * eoc;
                 g310 =  -346.844 +  1582.851 * *em -  2415.925 * emsq +  1246.113 * eoc;
                 g322 =  -342.585 +  1554.908 * *em -  2366.899 * emsq +  1215.972 * eoc;
                 g410 = -1052.797 +  4758.686 * *em -  7193.992 * emsq +  3651.957 * eoc;
                 g422 = -3581.690 + 16178.110 * *em - 24462.770 * emsq + 12422.520 * eoc;
                 if (*em > 0.715)
                     g520 =-5149.66 + 29936.92 * *em - 54087.36 * emsq + 31324.56 * eoc;
                   else
                     g520 = 1464.74 -  4664.75 * *em +  3763.64 * emsq;
               }
             if (*em < 0.7)
               {
                 g533 = -919.22770 + 4988.6100 * *em - 9064.7700 * emsq + 5542.21  * eoc;
                 g521 = -822.71072 + 4568.6173 * *em - 8491.4146 * emsq + 5337.524 * eoc;
                 g532 = -853.66600 + 4690.2500 * *em - 8624.7700 * emsq + 5341.4  * eoc;
               }
               else
               {
                 g533 =-37995.780 + 161616.52 * *em - 229838.20 * emsq + 109377.94 * eoc;
                 g521 =-51752.104 + 218913.95 * *em - 309468.16 * emsq + 146349.42 * eoc;
                 g532 =-40023.880 + 170470.89 * *em - 242699.48 * emsq + 115605.82 * eoc;
               }

             sini2=  sinim * sinim;
             f220 =  0.75 * (1.0 + 2.0 * cosim+cosisq);
             f221 =  1.5 * sini2;
             f321 =  1.875 * sinim  *  (1.0 - 2.0 * cosim - 3.0 * cosisq);
             f322 = -1.875 * sinim  *  (1.0 + 2.0 * cosim - 3.0 * cosisq);
             f441 = 35.0 * sini2 * f220;
             f442 = 39.3750 * sini2 * sini2;
             f522 =  9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim- 5.0 * cosisq) +
                     0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq) );
             f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim +
                    10.0 * cosisq) + 6.56250012 * (1.0+2.0 * cosim - 3.0 * cosisq));
             f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim+cosisq *
                    (-12.0 + 8.0 * cosim + 10.0 * cosisq));
             f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim+cosisq *
                    (12.0 + 8.0 * cosim - 10.0 * cosisq));
             xno2  =  *nm * *nm;
             ainv2 =  aonv * aonv;
             temp1 =  3.0 * xno2 * ainv2;
             temp  =  temp1 * root22;
             *d2201 =  temp * f220 * g201;
             *d2211 =  temp * f221 * g211;
             temp1 =  temp1 * aonv;
             temp  =  temp1 * root32;
             *d3210 =  temp * f321 * g310;
             *d3222 =  temp * f322 * g322;
             temp1 =  temp1 * aonv;
             temp  =  2.0 * temp1 * root44;
             *d4410 =  temp * f441 * g410;
             *d4422 =  temp * f442 * g422;
             temp1 =  temp1 * aonv;
             temp  =  temp1 * root52;
             *d5220 =  temp * f522 * g520;
             *d5232 =  temp * f523 * g532;
             temp  =  2.0 * temp1 * root54;
             *d5421 =  temp * f542 * g521;
             *d5433 =  temp * f543 * g533;
             *xlamo =  fmod(mo + nodeo + nodeo-theta - theta, twopi);
             *xfact =  mdot + *dmdt + 2.0 * (nodedot + *dnodt - rptim) - no;
             *em    = emo;
             emsq  = emsqo;
           }

         /* ---------------- synchronous resonance terms -------------- */
         if (*irez == 1)
           {
             g200  = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
             g310  = 1.0 + 2.0 * emsq;
             g300  = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
             f220  = 0.75 * (1.0 + cosim) * (1.0 + cosim);
             f311  = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
             f330  = 1.0 + cosim;
             f330  = 1.875 * f330 * f330 * f330;
             *del1  = 3.0 * *nm * *nm * aonv * aonv;
             *del2  = 2.0 * *del1 * f220 * g200 * q22;
             *del3  = 3.0 * *del1 * f330 * g300 * q33 * aonv;
             *del1  = *del1 * f311 * g310 * q31 * aonv;
             *xlamo = fmod(mo + nodeo + argpo - theta, twopi);
             *xfact = mdot + xpidot - rptim + *dmdt + *domdt + *dnodt - no;
           }

         /* ------------ for sgp4, initialize the integrator ---------- */
         *xli   = *xlamo;
         *xni   = no;
         *atime = 0.0;
         *nm    = no + *dndt;
       }

//#include "debug3.cpp"
}  // end dsinit

/*-----------------------------------------------------------------------------
*
*                           procedure dspace
*
*  this procedure provides deep space contributions to mean elements for
*    perturbing third body.  these effects have been averaged over one
*    revolution of the sun and moon.  for earth resonance effects, the
*    effects have been averaged over no revolutions of the satellite.
*    (mean motion)
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
*    dedt        -
*    del1, del2, del3  -
*    didt        -
*    dmdt        -
*    dnodt       -
*    domdt       -
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    argpo       - argument of perigee
*    argpdot     - argument of perigee dot (rate)
*    t           - time
*    tc          -
*    gsto        - gst
*    xfact       -
*    xlamo       -
*    no          - mean motion
*    atime       -
*    em          - eccentricity
*    ft          -
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    atime       -
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         -
*    nodem       - right ascension of ascending node
*    dndt        -
*    nm          - mean motion
*
*  locals        :
*    delt        -
*    ft          -
*    theta       -
*    x2li        -
*    x2omi       -
*    xl          -
*    xldot       -
*    xnddt       -
*    xndt        -
*    xomi        -
*
*  coupling      :
*    none        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

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
{
     const double twopi = 2.0 * PI;
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
     rptim = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
     stepp =    720.0;
     stepn =   -720.0;
     step2 = 259200.0;

     /* ----------- calculate deep space resonance effects ----------- */
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

     /* - update resonances : numerical (euler-maclaurin) integration - */
     /* ------------------------- epoch restart ----------------------  */
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
             /* ------------------- dot terms calculated ------------- */
             /* ----------- near - synchronous resonance terms ------- */
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
                 /* --------- near - half-day resonance terms -------- */
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

             /* ----------------------- integrator ------------------- */
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

//#include "debug4.cpp"
}  // end dsspace

/*-----------------------------------------------------------------------------
*
*                           procedure initl
*
*  this procedure initializes the spg4 propagator. all the initialization is
*    consolidated here instead of having multiple loops inside other routines.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    ecco        - eccentricity                           0.0 - 1.0
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    inclo       - inclination of satellite
*    no          - mean motion of satellite
*    satn        - satellite number
*
*  outputs       :
*    ainv        - 1.0 / a
*    ao          - semi major axis
*    con41       -
*    con42       - 1.0 - 5.0 cos(i)
*    cosio       - cosine of inclination
*    cosio2      - cosio squared
*    eccsq       - eccentricity squared
*    method      - flag for deep space                    'd', 'n'
*    omeosq      - 1.0 - ecco * ecco
*    posq        - semi-parameter squared
*    rp          - radius of perigee
*    rteosq      - square root of (1.0 - ecco*ecco)
*    sinio       - sine of inclination
*    gsto        - gst at time of observation               rad
*    no          - mean motion of satellite
*
*  locals        :
*    ak          -
*    d1          -
*    del         -
*    adel        -
*    po          -
*
*  coupling      :
*    getgravconst
*    gstime      - find greenwich sidereal time from the julian date
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

void initl
     (
       int satn,      int whichconst,
       double ecco,   double epoch,  double inclo,   double* no,
       char* method,
       double* ainv,  double* ao,    double* con41,  double* con42, double* cosio,
       double* cosio2,double* eccsq, double* omeosq, double* posq,
       double* rp,    double* rteosq,double* sinio , double* gsto,
       char opsmode
     )
{
     /* --------------------- local variables ------------------------ */
     double ak, d1, del, adel, po, x2o3, j2, xke,
            tumin, mu, radiusearthkm, j3, j4, j3oj2;

     // sgp4fix use old way of finding gst
     double ds70;
     double ts70, tfrac, c1, thgr70, fk5r, c1p2p;
     const double twopi = 2.0 * PI;

     /* ----------------------- earth constants ---------------------- */
     // sgp4fix identify constants and allow alternate values
     //getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
     // TODO: replace by macros
     tumin = TUMIN;
     mu = GM;
     radiusearthkm = RE;
     xke = XKE;
     j2 = J2;
     j3 = J3;
     j4 = J4;
     j3oj2 = J3DIVJ2;

     x2o3   = 2.0 / 3.0;

     /* ------------- calculate auxillary epoch quantities ---------- */
     *eccsq  = ecco * ecco;
     *omeosq = 1.0 - *eccsq;
     *rteosq = sqrt(*omeosq);
     *cosio  = cos(inclo);
     *cosio2 = *cosio * *cosio;

     /* ------------------ un-kozai the mean motion ----------------- */
     ak    = pow(xke / *no, x2o3);
     d1    = 0.75 * j2 * (3.0 * *cosio2 - 1.0) / (*rteosq * *omeosq);
     del   = d1 / (ak * ak);
     adel  = ak * (1.0 - del * del - del *
             (1.0 / 3.0 + 134.0 * del * del / 81.0));
     del   = d1/(adel * adel);
     *no    = *no / (1.0 + del);

     *ao    = pow(xke / *no, x2o3);
     *sinio = sin(inclo);
     po    = *ao * *omeosq;
     *con42 = 1.0 - 5.0 * *cosio2;
     *con41 = -*con42-*cosio2-*cosio2;
     *ainv  = 1.0 / *ao;
     *posq  = po * po;
     *rp    = *ao * (1.0 - ecco);
     *method = 'n';

     // sgp4fix modern approach to finding sidereal time
     if (opsmode == 'a')
        {
         // sgp4fix use old way of finding gst
         // count integer number of days from 0 jan 1970
         ts70  = epoch - 7305.0;
         ds70 = floor(ts70 + 1.0e-8);
         tfrac = ts70 - ds70;
         // find greenwich location at epoch
         c1    = 1.72027916940703639e-2;
         thgr70= 1.7321343856509374;
         fk5r  = 5.07551419432269442e-15;
         c1p2p = c1 + twopi;
         *gsto  = fmod( thgr70 + c1*ds70 + c1p2p*tfrac + ts70*ts70*fk5r, twopi);
         if ( *gsto < 0.0 )
             *gsto = *gsto + twopi;
       }
       else
        //*gsto = gstime(epoch + 2433281.5);
         *gsto = jul2gst(epoch + 2433281.5);


//#include "debug5.cpp"
}  // end initl

/*-----------------------------------------------------------------------------
*
*                             procedure sgp4init
*
*  this procedure initializes variables for sgp4.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    opsmode     - mode of operation afspc or improved 'a', 'i'
*    whichconst  - which set of constants to use  72, 84
*    satn        - satellite number
*    bstar       - sgp4 type drag coefficient              kg/m2er
*    ecco        - eccentricity
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    argpo       - argument of perigee (output if ds)
*    inclo       - inclination
*    mo          - mean anomaly (output if ds)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*
*  outputs       :
*    satrec      - common values for subsequent calls
*    return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    cnodm  , snodm  , cosim  , sinim  , cosomm , sinomm
*    cc1sq  , cc2    , cc3
*    coef   , coef1
*    cosio4      -
*    day         -
*    dndt        -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    eeta        -
*    etasq       -
*    gam         -
*    argpm       - argument of perigee
*    nodem       -
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    perige      - perigee
*    pinvsq      -
*    psisq       -
*    qzms24      -
*    rtemsq      -
*    s1, s2, s3, s4, s5, s6, s7          -
*    sfour       -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
*    sz1, sz2, sz3
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    tc          -
*    temp        -
*    temp1, temp2, temp3       -
*    tsi         -
*    xpidot      -
*    xhdot1      -
*    z1, z2, z3          -
*    z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*
*  coupling      :
*    getgravconst-
*    initl       -
*    dscom       -
*    dpper       -
*    dsinit      -
*    sgp4        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

bool sgp4init
     (
       int whichconst, char opsmode,   const int satn,     const double epoch,
       const double xbstar,  const double xecco, const double xargpo,
       const double xinclo,  const double xmo,   const double xno,
       const double xnodeo,  sat* sat
     )
{
     /* --------------------- local variables ------------------------ */
     double ao, ainv,   con42, cosio, sinio, cosio2, eccsq,
          omeosq, posq,   rp,     rteosq,
          cnodm , snodm , cosim , sinim , cosomm, sinomm, cc1sq ,
          cc2   , cc3   , coef  , coef1 , cosio4, day   , dndt  ,
          em    , emsq  , eeta  , etasq , gam   , argpm , nodem ,
          inclm , mm    , nm    , perige, pinvsq, psisq , qzms24,
          rtemsq, s1    , s2    , s3    , s4    , s5    , s6    ,
          s7    , sfour , ss1   , ss2   , ss3   , ss4   , ss5   ,
          ss6   , ss7   , sz1   , sz2   , sz3   , sz11  , sz12  ,
          sz13  , sz21  , sz22  , sz23  , sz31  , sz32  , sz33  ,
          tc    , temp  , temp1 , temp2 , temp3 , tsi   , xpidot,
          xhdot1, z1    , z2    , z3    , z11   , z12   , z13   ,
          z21   , z22   , z23   , z31   , z32   , z33,
          qzms2t, ss, j2, j3oj2, j4, x2o3, r[3], v[3],
          tumin, mu, radiusearthkm, xke, j3, delmotemp, qzms2ttemp, qzms24temp;

     /* ------------------------ initialization --------------------- */
     // sgp4fix divisor for divide by zero check on inclination
     // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
     // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
     const double temp4    =   1.5e-12;
     // TODO: Get rid of zeroing out and use a compount initializer
     /* ----------- set all near earth variables to zero ------------ */
     sat->isimp   = 0;   sat->method = 'n'; sat->aycof    = 0.0;
     sat->x3thm1   = 0.0; sat->cc1    = 0.0; sat->cc4      = 0.0;
     sat->cc5     = 0.0; sat->d2     = 0.0; sat->d3       = 0.0;
     sat->d4      = 0.0; sat->delmo  = 0.0; sat->eta      = 0.0;
     sat->argpdot = 0.0; sat->omgcof = 0.0; sat->sinmao   = 0.0;
     sat->t       = 0.0; sat->t2cof  = 0.0; sat->t3cof    = 0.0;
     sat->t4cof   = 0.0; sat->t5cof  = 0.0; sat->x1mth2   = 0.0;
     sat->x7thm1  = 0.0; sat->mdot   = 0.0; sat->nodedot  = 0.0;
     sat->xlcof   = 0.0; sat->xmcof  = 0.0; sat->nodecf   = 0.0;

     /* ----------- set all deep space variables to zero ------------ */
     sat->irez  = 0;   sat->d2201 = 0.0; sat->d2211 = 0.0;
     sat->d3210 = 0.0; sat->d3222 = 0.0; sat->d4410 = 0.0;
     sat->d4422 = 0.0; sat->d5220 = 0.0; sat->d5232 = 0.0;
     sat->d5421 = 0.0; sat->d5433 = 0.0; sat->dedt  = 0.0;
     sat->del1  = 0.0; sat->del2  = 0.0; sat->del3  = 0.0;
     sat->didt  = 0.0; sat->dmdt  = 0.0; sat->dnodt = 0.0;
     sat->domdt = 0.0; sat->e3    = 0.0; sat->ee2   = 0.0;
     sat->peo   = 0.0; sat->pgho  = 0.0; sat->pho   = 0.0;
     sat->pinco = 0.0; sat->plo   = 0.0; sat->se2   = 0.0;
     sat->se3   = 0.0; sat->sgh2  = 0.0; sat->sgh3  = 0.0;
     sat->sgh4  = 0.0; sat->sh2   = 0.0; sat->sh3   = 0.0;
     sat->si2   = 0.0; sat->si3   = 0.0; sat->sl2   = 0.0;
     sat->sl3   = 0.0; sat->sl4   = 0.0; sat->gsto  = 0.0;
     sat->xfact = 0.0; sat->xgh2  = 0.0; sat->xgh3  = 0.0;
     sat->xgh4  = 0.0; sat->xh2   = 0.0; sat->xh3   = 0.0;
     sat->xi2   = 0.0; sat->xi3   = 0.0; sat->xl2   = 0.0;
     sat->xl3   = 0.0; sat->xl4   = 0.0; sat->xlamo = 0.0;
     sat->zmol  = 0.0; sat->zmos  = 0.0; sat->atime = 0.0;
     sat->xli   = 0.0; sat->xni   = 0.0;

     // TODO: Only to test legacy math!
     sat->argpdot = sat->omegaprime; // omegaprime
     //sat->gsto = sat->GSTo;    // GSTo
     sat->bstar = sat->Bstar;   // Bstar
     sat->ecco = sat->eccentricity;    // e
     sat->argpo = sat->argument_perigee;   // omega
     sat->inclo = sat->inclination;   // i
     sat->mo = sat->mean_anomaly;      // Mo
     sat->nodeo = sat->right_asc_node;   // alpha

     // Convert to SGP4 units
     // ---- find no, ndot, nddot ----
     sat->mean_motion   = sat->mean_motion / RPD2RADPM; //* rad/min

     /* ------------------------ earth constants ----------------------- */
     // sgp4fix identify constants and allow alternate values
     // TODO: replace by macros
     tumin = TUMIN;
     mu = GM;
     radiusearthkm = RE;
     xke = XKE;
     j2 = J2;
     j3 = J3;
     j4 = J4;
     j3oj2 = J3DIVJ2;

     // ---- convert to sgp4 units ----
     sat->a    = pow( sat->mean_motion*tumin , (-2.0/3.0) );
     sat->mean_motion_dt2 = sat->mean_motion_dt2  / (RPD2RADPM*1440.0);  //* ? * minperday
     sat->mean_motion_ddt6= sat->mean_motion_ddt6 / (RPD2RADPM*1440.0*1440);

     // ---- find standard orbital elements ----
     sat->inclo = sat->inclo  * DEG2RAD;
     sat->nodeo = sat->nodeo  * DEG2RAD;
     sat->argpo = sat->argpo  * DEG2RAD;
     sat->mo    = sat->mo     * DEG2RAD;

     sat->alta = sat->a*(1.0 + sat->ecco) - 1.0;
     sat->altp = sat->a*(1.0 - sat->ecco) - 1.0;
     // sgp4fix - note the following variables are also passed directly via sat->
     // it is possible to streamline the sgp4init call by deleting the "x"
     // variables, but the user would need to set the sat->* values first. we
     // include the additional assignments in case twoline2rv is not used.
     //sat->bstar   = xbstar;
     //sat->ecco    = xecco;
     //sat->argpo   = xargpo;
     //sat->inclo   = xinclo;
     //sat->mo      = xmo;
     //sat->no      = xno;
     //sat->nodeo   = xnodeo;

     // sgp4fix add opsmode
     sat->operationmode = opsmode;

     ss     = 78.0 / radiusearthkm + 1.0;
     // sgp4fix use multiply for speed instead of pow
     qzms2ttemp = (120.0 - 78.0) / radiusearthkm;
     qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
     x2o3   =  2.0 / 3.0;

     sat->init = 'y';
     sat->t  = 0.0;

     //print_orbit(sat, "sgp4init pre-initl");
     initl
         (
           satn, whichconst, sat->ecco, epoch, sat->inclo, &sat->mean_motion, &sat->method,
           &ainv, &ao, &sat->x3thm1, &con42, &cosio, &cosio2, &eccsq, &omeosq,
           &posq, &rp, &rteosq, &sinio, &sat->gsto, sat->operationmode
         );
     sat->error = 0;

     //print_orbit(sat, "sgp4init post-initl");

     // sgp4fix remove this check as it is unnecessary
     // the mrt check in sgp4 handles decaying satellite cases even if the starting
     // condition is below the surface of te earth
//     if (rp < 1.0)
//       {
//         printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
//         sat->error = 5;
//       }

     if ((omeosq >= 0.0 ) || ( sat->mean_motion >= 0.0))
       {
         sat->isimp = 0;
         if (rp < (220.0 / radiusearthkm + 1.0))
             sat->isimp = 1;
         sfour  = ss;
         qzms24 = qzms2t;
         perige = (rp - 1.0) * radiusearthkm;

         /* - for perigees below 156 km, s and qoms2t are altered - */
         if (perige < 156.0)
           {
             sfour = perige - 78.0;
             if (perige < 98.0)
                 sfour = 20.0;
             // sgp4fix use multiply for speed instead of pow
             qzms24temp =  (120.0 - sfour) / radiusearthkm;
             qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
             sfour  = sfour / radiusearthkm + 1.0;
           }
         pinvsq = 1.0 / posq;

         tsi  = 1.0 / (ao - sfour);
         sat->eta  = ao * sat->ecco * tsi;
         etasq = sat->eta * sat->eta;
         eeta  = sat->ecco * sat->eta;
         psisq = fabs(1.0 - etasq);
         coef  = qzms24 * pow(tsi, 4.0);
         coef1 = coef / pow(psisq, 3.5);
         cc2   = coef1 * sat->mean_motion * (ao * (1.0 + 1.5 * etasq + eeta *
                        (4.0 + etasq)) + 0.375 * j2 * tsi / psisq * sat->x3thm1 *
                        (8.0 + 3.0 * etasq * (8.0 + etasq)));
         sat->cc1   = sat->bstar * cc2;
         cc3   = 0.0;
         if (sat->ecco > 1.0e-4)
             cc3 = -2.0 * coef * tsi * j3oj2 * sat->mean_motion * sinio / sat->ecco;
         sat->x1mth2 = 1.0 - cosio2;
         sat->cc4    = 2.0* sat->mean_motion * coef1 * ao * omeosq *
                           (sat->eta * (2.0 + 0.5 * etasq) + sat->ecco *
                           (0.5 + 2.0 * etasq) - j2 * tsi / (ao * psisq) *
                           (-3.0 * sat->x3thm1 * (1.0 - 2.0 * eeta + etasq *
                           (1.5 - 0.5 * eeta)) + 0.75 * sat->x1mth2 *
                           (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * sat->argpo)));
         sat->cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
                        (etasq + eeta) + eeta * etasq);
         cosio4 = cosio2 * cosio2;
         temp1  = 1.5 * j2 * pinvsq * sat->mean_motion;
         temp2  = 0.5 * temp1 * j2 * pinvsq;
         temp3  = -0.46875 * j4 * pinvsq * pinvsq * sat->mean_motion;
         sat->mdot     = sat->mean_motion + 0.5 * temp1 * rteosq * sat->x3thm1 + 0.0625 *
                            temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
         sat->argpdot  = -0.5 * temp1 * con42 + 0.0625 * temp2 *
                             (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                             temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
         xhdot1            = -temp1 * cosio;
         sat->nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
                              2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
         xpidot            =  sat->argpdot+ sat->nodedot;
         sat->omgcof   = sat->bstar * cc3 * cos(sat->argpo);
         sat->xmcof    = 0.0;
         if (sat->ecco > 1.0e-4)
             sat->xmcof = -x2o3 * coef * sat->bstar / eeta;
         sat->nodecf = 3.5 * omeosq * xhdot1 * sat->cc1;
         sat->t2cof   = 1.5 * sat->cc1;
         // sgp4fix for divide by zero with xinco = 180 deg
         if (fabs(cosio+1.0) > 1.5e-12)
             sat->xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
           else
             sat->xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
         sat->aycof   = -0.5 * j3oj2 * sinio;
         // sgp4fix use multiply for speed instead of pow
         delmotemp = 1.0 + sat->eta * cos(sat->mo);
         sat->delmo   = delmotemp * delmotemp * delmotemp;
         sat->sinmao  = sin(sat->mo);
         sat->x7thm1  = 7.0 * cosio2 - 1.0;
         //print_orbit(sat, "sgp4init pre-deepspace");
         /* --------------- deep space initialization ------------- */
         if ((2*PI / sat->mean_motion) >= 225.0)
           {
             sat->method = 'd';
             sat->isimp  = 1;
             tc    =  0.0;
             inclm = sat->inclo;

             dscom
                 (
                   epoch, sat->ecco, sat->argpo, tc, sat->inclo, sat->nodeo,
                   sat->mean_motion, &snodm, &cnodm,  &sinim, &cosim, &sinomm,     &cosomm,
                   &day, &sat->e3, &sat->ee2, &em,         &emsq, &gam,
                   &sat->peo,  &sat->pgho,   &sat->pho, &sat->pinco,
                   &sat->plo,  &rtemsq,        &sat->se2, &sat->se3,
                   &sat->sgh2, &sat->sgh3,   &sat->sgh4,
                   &sat->sh2,  &sat->sh3,    &sat->si2, &sat->si3,
                   &sat->sl2,  &sat->sl3,    &sat->sl4, &s1, &s2, &s3, &s4, &s5,
                   &s6,   &s7,   &ss1,  &ss2,  &ss3,  &ss4,  &ss5,  &ss6,  &ss7, &sz1, &sz2, &sz3,
                   &sz11, &sz12, &sz13, &sz21, &sz22, &sz23, &sz31, &sz32, &sz33,
                   &sat->xgh2, &sat->xgh3,   &sat->xgh4, &sat->xh2,
                   &sat->xh3,  &sat->xi2,    &sat->xi3,  &sat->xl2,
                   &sat->xl3,  &sat->xl4,    &nm, &z1, &z2, &z3, &z11,
                   &z12, &z13, &z21, &z22, &z23, &z31, &z32, &z33,
                   &sat->zmol, &sat->zmos
                 );
             //print_orbit(sat, "sgp4init post-dscom");
             dpper
                 (
                   sat->e3, sat->ee2, sat->peo, sat->pgho,
                   sat->pho, sat->pinco, sat->plo, sat->se2,
                   sat->se3, sat->sgh2, sat->sgh3, sat->sgh4,
                   sat->sh2, sat->sh3, sat->si2,  sat->si3,
                   sat->sl2, sat->sl3, sat->sl4,  sat->t,
                   sat->xgh2,sat->xgh3,sat->xgh4, sat->xh2,
                   sat->xh3, sat->xi2, sat->xi3,  sat->xl2,
                   sat->xl3, sat->xl4, sat->zmol, sat->zmos, inclm, sat->init,
                   &sat->ecco, &sat->inclo, &sat->nodeo, &sat->argpo, &sat->mo,
                   sat->operationmode
                 );
             //print_orbit(sat, "sgp4init post-dpper");
             argpm  = 0.0;
             nodem  = 0.0;
             mm     = 0.0;

             dsinit
                 (
                   whichconst,
                   cosim, emsq, sat->argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
                   ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, sat->t, tc,
                   sat->gsto, sat->mo, sat->mdot, sat->mean_motion, sat->nodeo,
                   sat->nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
                   sat->ecco, eccsq, &em, &argpm, &inclm, &mm, &nm, &nodem,
                   &sat->irez,  &sat->atime,
                   &sat->d2201, &sat->d2211, &sat->d3210, &sat->d3222 ,
                   &sat->d4410, &sat->d4422, &sat->d5220, &sat->d5232,
                   &sat->d5421, &sat->d5433, &sat->dedt,  &sat->didt,
                   &sat->dmdt,  &dndt,         &sat->dnodt, &sat->domdt ,
                   &sat->del1,  &sat->del2,  &sat->del3,  &sat->xfact,
                   &sat->xlamo, &sat->xli,   &sat->xni
                 );
             //print_orbit(sat, "sgp4init post-dsinit");
           }

       /* ----------- set variables if not deep space ----------- */
       if (sat->isimp != 1)
         {
           cc1sq          = sat->cc1 * sat->cc1;
           sat->d2    = 4.0 * ao * tsi * cc1sq;
           temp           = sat->d2 * tsi * sat->cc1 / 3.0;
           sat->d3    = (17.0 * ao + sfour) * temp;
           sat->d4    = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) *
                            sat->cc1;
           sat->t3cof = sat->d2 + 2.0 * cc1sq;
           sat->t4cof = 0.25 * (3.0 * sat->d3 + sat->cc1 *
                            (12.0 * sat->d2 + 10.0 * cc1sq));
           sat->t5cof = 0.2 * (3.0 * sat->d4 +
                            12.0 * sat->cc1 * sat->d3 +
                            6.0 * sat->d2 * sat->d2 +
                            15.0 * cc1sq * (2.0 * sat->d2 + cc1sq));
         }
       } // if omeosq = 0 ...
     //print_orbit(sat, "sgp4init post-simple model");
       /* finally propogate to zero epoch to initialize all others. */
       // sgp4fix take out check to let satellites process until they are actually below earth surface
//       if(sat->error == 0)
       sgp4(1, sat, 0.0, r, v);

       sat->init = 'n';
       //print_orbit(sat, "sgp4init post-sgp4 call");
//#include "debug6.cpp"
       //sgp4fix return boolean. sat->error contains any error codes
       return true;
}  // end sgp4init

/*-----------------------------------------------------------------------------
*
*                             procedure sgp4
*
*  this procedure is the sgp4 prediction model from space command. this is an
*    updated and combined version of sgp4 and sdp4, which were originally
*    published separately in spacetrack report #3. this version follows the
*    methodology from the aiaa paper (2006) describing the history and
*    development of the code.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satrec  - initialised structure from sgp4init() call.
*    tsince  - time eince epoch (minutes)
*
*  outputs       :
*    r           - position vector                     km
*    v           - velocity                            km/sec
*  return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    am          -
*    axnl, aynl        -
*    betal       -
*    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
*    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
*    cosisq  , cossu   , sinsu   , cosu    , sinu
*    delm        -
*    delomg      -
*    dndt        -
*    eccm        -
*    emsq        -
*    ecose       -
*    el2         -
*    eo1         -
*    eccp        -
*    esine       -
*    argpm       -
*    argpp       -
*    omgadf      -
*    pl          -
*    r           -
*    rtemsq      -
*    rdotl       -
*    rl          -
*    rvdot       -
*    rvdotl      -
*    su          -
*    t2  , t3   , t4    , tc
*    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
*    u   , ux   , uy    , uz     , vx     , vy     , vz
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right asc of ascending node
*    xinc        -
*    xincp       -
*    xl          -
*    xlm         -
*    mp          -
*    xmdf        -
*    xmx         -
*    xmy         -
*    nodedf      -
*    xnode       -
*    nodep       -
*    np          -
*
*  coupling      :
*    getgravconst-
*    dpper
*    dpspace
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

bool sgp4
     (
       int whichconst, sat* sat,  double tsince,
       double r[3],  double v[3]
     )
{
     double am   , axnl  , aynl , betal ,  cosim , cnod  ,
         cos2u, coseo1, cosi , cosip ,  cosisq, cossu , cosu,
         delm , delomg, em   , emsq  ,  ecose , el2   , eo1 ,
         ep   , esine , argpm, argpp ,  argpdf, pl,     mrt = 0.0,
         mvt  , rdotl , rl   , rvdot ,  rvdotl, sinim ,
         sin2u, sineo1, sini , sinip ,  sinsu , sinu  ,
         snod , su    , t2   , t3    ,  t4    , tem5  , temp,
         temp1, temp2 , tempa, tempe ,  templ , u     , ux  ,
         uy   , uz    , vx   , vy    ,  vz    , inclm , mm  ,
         nm   , nodem, xinc , xincp ,  xl    , xlm   , mp  ,
         xmdf , xmx   , xmy  , nodedf, xnode , nodep, tc  , dndt,
         twopi, x2o3  , j2   , j3    , tumin, j4 , xke   , j3oj2, radiusearthkm,
         mu, vkmpersec, delmtemp;
     int ktr;

     /* ------------------ set mathematical constants --------------- */
     // sgp4fix divisor for divide by zero check on inclination
     // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
     // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
     const double temp4 =   1.5e-12;
     twopi = 2.0 * PI;
     x2o3  = 2.0 / 3.0;
     // sgp4fix identify constants and allow alternate values
     // TODO: replace by macros
     tumin = TUMIN;
     mu = GM;
     radiusearthkm = RE;
     xke = XKE;
     j2 = J2;
     j3 = J3;
     j4 = J4;
     j3oj2 = J3DIVJ2;
     vkmpersec     = radiusearthkm * xke/60.0;

     /* --------------------- clear sgp4 error flag ----------------- */
     sat->t     = tsince;
     sat->error = 0;

     /* ------- update for secular gravity and atmospheric drag ----- */
     xmdf    = sat->mo + sat->mdot * sat->t;
     argpdf  = sat->argpo + sat->argpdot * sat->t;
     nodedf  = sat->nodeo + sat->nodedot * sat->t;
     argpm   = argpdf;
     mm      = xmdf;
     t2      = sat->t * sat->t;
     nodem   = nodedf + sat->nodecf * t2;
     tempa   = 1.0 - sat->cc1 * sat->t;
     tempe   = sat->bstar * sat->cc4 * sat->t;
     templ   = sat->t2cof * t2;

     if (sat->isimp != 1)
       {
         delomg = sat->omgcof * sat->t;
         // sgp4fix use mutliply for speed instead of pow
         delmtemp =  1.0 + sat->eta * cos(xmdf);
         delm   = sat->xmcof *
                  (delmtemp * delmtemp * delmtemp -
                  sat->delmo);
         temp   = delomg + delm;
         mm     = xmdf + temp;
         argpm  = argpdf - temp;
         t3     = t2 * sat->t;
         t4     = t3 * sat->t;
         tempa  = tempa - sat->d2 * t2 - sat->d3 * t3 -
                          sat->d4 * t4;
         tempe  = tempe + sat->bstar * sat->cc5 * (sin(mm) -
                          sat->sinmao);
         templ  = templ + sat->t3cof * t3 + t4 * (sat->t4cof +
                          sat->t * sat->t5cof);
       }

     nm    = sat->mean_motion;
     em    = sat->ecco;
     inclm = sat->inclo;
     if (sat->method == 'd')
       {
         tc = sat->t;
         dspace
             (
               sat->irez,
               sat->d2201, sat->d2211, sat->d3210,
               sat->d3222, sat->d4410, sat->d4422,
               sat->d5220, sat->d5232, sat->d5421,
               sat->d5433, sat->dedt,  sat->del1,
               sat->del2,  sat->del3,  sat->didt,
               sat->dmdt,  sat->dnodt, sat->domdt,
               sat->argpo, sat->argpdot, sat->t, tc,
               sat->gsto, sat->xfact, sat->xlamo,
               sat->mean_motion, &sat->atime,
               &em, &argpm, &inclm, &sat->xli, &mm, &sat->xni,
               &nodem, &dndt, &nm
             );
       } // if method = d

     if (nm <= 0.0)
       {
//         printf("# error nm %f\n", nm);
         sat->error = 2;
         // sgp4fix add return
         return false;
       }
     am = pow((xke / nm),x2o3) * tempa * tempa;
     nm = xke / pow(am, 1.5);
     em = em - tempe;

     // fix tolerance for error recognition
     // sgp4fix am is fixed from the previous nm check
     if ((em >= 1.0) || (em < -0.001)/* || (am < 0.95)*/ )
       {
//         printf("# error em %f\n", em);
         sat->error = 1;
         // sgp4fix to return if there is an error in eccentricity
         return false;
       }
     // sgp4fix fix tolerance to avoid a divide by zero
     if (em < 1.0e-6)
         em  = 1.0e-6;
     mm     = mm + sat->mean_motion * templ;
     xlm    = mm + argpm + nodem;
     emsq   = em * em;
     temp   = 1.0 - emsq;

     nodem  = fmod(nodem, twopi);
     argpm  = fmod(argpm, twopi);
     xlm    = fmod(xlm, twopi);
     mm     = fmod(xlm - argpm - nodem, twopi);

     /* ----------------- compute extra mean quantities ------------- */
     sinim = sin(inclm);
     cosim = cos(inclm);

     /* -------------------- add lunar-solar periodics -------------- */
     ep     = em;
     xincp  = inclm;
     argpp  = argpm;
     nodep  = nodem;
     mp     = mm;
     sinip  = sinim;
     cosip  = cosim;
     if (sat->method == 'd')
       {
         dpper
             (
               sat->e3,   sat->ee2,  sat->peo,
               sat->pgho, sat->pho,  sat->pinco,
               sat->plo,  sat->se2,  sat->se3,
               sat->sgh2, sat->sgh3, sat->sgh4,
               sat->sh2,  sat->sh3,  sat->si2,
               sat->si3,  sat->sl2,  sat->sl3,
               sat->sl4,  sat->t,    sat->xgh2,
               sat->xgh3, sat->xgh4, sat->xh2,
               sat->xh3,  sat->xi2,  sat->xi3,
               sat->xl2,  sat->xl3,  sat->xl4,
               sat->zmol, sat->zmos, sat->inclo,
               'n', &ep, &xincp, &nodep, &argpp, &mp, sat->operationmode
             );
         if (xincp < 0.0)
           {
             xincp  = -xincp;
             nodep = nodep + PI;
             argpp  = argpp - PI;
           }
         if ((ep < 0.0 ) || ( ep > 1.0))
           {
//            printf("# error ep %f\n", ep);
             sat->error = 3;

             // sgp4fix add return
             return false;
           }
       } // if method = d

     /* -------------------- long period periodics ------------------ */
     if (sat->method == 'd')
       {
         sinip =  sin(xincp);
         cosip =  cos(xincp);
         sat->aycof = -0.5*j3oj2*sinip;
         // sgp4fix for divide by zero for xincp = 180 deg
         if (fabs(cosip+1.0) > 1.5e-12)
             sat->xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
           else
             sat->xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
       }
     axnl = ep * cos(argpp);
     temp = 1.0 / (am * (1.0 - ep * ep));
     aynl = ep* sin(argpp) + temp * sat->aycof;
     xl   = mp + argpp + nodep + temp * sat->xlcof * axnl;

     /* --------------------- solve kepler's equation --------------- */
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

     /* ------------- short period preliminary quantities ----------- */
     ecose = axnl*coseo1 + aynl*sineo1;
     esine = axnl*sineo1 - aynl*coseo1;
     el2   = axnl*axnl + aynl*aynl;
     pl    = am*(1.0-el2);
     if (pl < 0.0)
       {
//         printf("# error pl %f\n", pl);
         sat->error = 4;
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

         /* -------------- update for short period periodics ------------ */
         if (sat->method == 'd')
           {
             cosisq                 = cosip * cosip;
             sat->x3thm1  = 3.0*cosisq - 1.0;
             sat->x1mth2 = 1.0 - cosisq;
             sat->x7thm1 = 7.0*cosisq - 1.0;
           }
         mrt   = rl * (1.0 - 1.5 * temp2 * betal * sat->x3thm1) +
                 0.5 * temp1 * sat->x1mth2 * cos2u;
         su    = su - 0.25 * temp2 * sat->x7thm1 * sin2u;
         xnode = nodep + 1.5 * temp2 * cosip * sin2u;
         xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
         mvt   = rdotl - nm * temp1 * sat->x1mth2 * sin2u / xke;
         rvdot = rvdotl + nm * temp1 * (sat->x1mth2 * cos2u +
                 1.5 * sat->x3thm1) / xke;

         /* --------------------- orientation vectors ------------------- */
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

         /* --------- position and velocity (in km and km/sec) ---------- */
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
//         printf("# decay condition %11.6f \n",mrt);
         sat->error = 6;
         return false;
       }

//#include "debug7.cpp"
     return true;
}  // end sgp4

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
       double r[3], double v[3], double mu,
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
         c1 = magv*magv - mu  / magr;
         rdotv = dot2( r,v );
         for (i= 0; i <= 2; i++)
             ebar[i]= (c1*r[i] - rdotv*v[i]) / mu;
         *ecc = mag2( ebar );

         // ------------  find *a e and semi-latus rectum   ----------
         sme= ( magv*magv*0.5  ) - ( mu  / magr );
         if ( fabs( sme ) > small )
             *a= -mu  / (2.0 *sme);
           else
             *a= infinite;
         *p = magh*magh / mu;

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
