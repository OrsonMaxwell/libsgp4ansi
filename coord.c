/*
 * transform.c - Coordinate transformation routines for libsgp4ansi.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#include <math.h>

#include "libsgp4ansi.h"
#include "const.h"
#include "coord.h"
#include "epoch.h"
#include "vector.h"

// ************************************************************************* //
//                                 MISC MATH                                 //
// ************************************************************************* //
/*
 * Return unity multiplier with the sign of the argument
 *
 * Inputs:  arg    - Source of the sign information
 * Outputs: None
 * Returns: unity with argument sign
 */
int
signof(double arg)
{
  return (arg < 0.0)? -1 : 1;
}

// ************************************************************************* //
//                                   TIME                                    //
// ************************************************************************* //


// ************************************************************************* //
//                               COORDINATES                                 //
// ************************************************************************* //

/*
 * Transform position and velocity vectors from TEME to ECEF frame of reference
 *
 * Inputs:  posteme - Position vector in TEME frame, km
 *          velteme - Velocity vector in TEME frame, km/s
 *          julian  - Julian time of interest
 * Outputs: posteme - Position vector in ECEF frame, km
 *          velteme - Velocity vector in ECEF frame, km/s
 * Returns: None
 */
void
teme2ecef
(
    vec3* posteme,
    vec3* velteme,
    double julian,
    vec3* posecef,
    vec3* velecef
)
{
  // Greenwich Siderial Time, rad
  double GST = jul2gst(julian);

  // Pef - tod matrix
  double pef_tod[3][3] =
  {
    {cos(GST), -sin(GST),  0.0},
    {sin(GST),  cos(GST),  0.0},
    {0.0,       0.0,       1.0}
  };

  // Pseudo Earth fixed position vector
  double ppef[3];
  ppef[0] = pef_tod[0][0] * posteme->x + pef_tod[1][0] * posteme->y +
              pef_tod[2][0] * posteme->z;
  ppef[1] = pef_tod[0][1] * posteme->x + pef_tod[1][1] * posteme->y +
              pef_tod[2][1] * posteme->z;
  ppef[2] = pef_tod[0][2] * posteme->x + pef_tod[1][2] * posteme->y +
              pef_tod[2][2] * posteme->z;

  // Refer to IERS Bulletin - A (Vol. XXVIII No. 030)
  double MJD = julian - 2400000.5;
  double A   = 2 * PI * (MJD - 57226) / 365.25;
  double C   = 2 * PI * (MJD - 57226) / 435;

  // Polar motion coefficients, rad
  double xp;
  double yp;
  xp = (0.1033 + 0.0494*cos(A) + 0.0482*sin(A) + 0.0297*cos(C) + 0.0307*sin(C))
       * 4.84813681e-6;
  yp = (0.3498 + 0.0441*cos(A) - 0.0393*sin(A) + 0.0307*cos(C) - 0.0297*sin(C))
       * 4.84813681e-6;

  double pm[3][3] =
  {
    {cos(xp),             0.0,      -sin(xp)},
    {sin(xp) * sin(yp),   cos(yp),   cos(xp) * sin(yp)},
    {sin(xp) * cos(yp),  -sin(yp),   cos(xp) * cos(yp)}
  };

  // ECEF postion vector
  double omegaearth[3];
  posecef->i = pm[0][0] * ppef[0] + pm[1][0] * ppef[1] + pm[2][0] * ppef[2];
  posecef->j = pm[0][1] * ppef[0] + pm[1][1] * ppef[1] + pm[2][1] * ppef[2];
  posecef->k = pm[0][2] * ppef[0] + pm[1][2] * ppef[1] + pm[2][2] * ppef[2];

  // Earth angular rotation vector
  omegaearth[0] = 0.0;
  omegaearth[1] = 0.0;
  omegaearth[2] = 7.29211514670698e-05 * (1.0  - 0.002 / 86400.0);

  // Pseudo Earth Fixed velocity vector
  double vpef[3];
  vpef[0] = pef_tod[0][0] * velteme->x + pef_tod[1][0] * velteme->y
          + pef_tod[2][0] * velteme->z - (omegaearth[1] * ppef[2]
          - omegaearth[2] * ppef[1]);
  vpef[1] = pef_tod[0][1] * velteme->x + pef_tod[1][1] * velteme->y
          + pef_tod[2][1] * velteme->z - (omegaearth[2] * ppef[0]
          - omegaearth[0] * ppef[2]);
  vpef[2] = pef_tod[0][2] * velteme->x + pef_tod[1][2] * velteme->y
          + pef_tod[2][2] * velteme->z - (omegaearth[0] * ppef[1]
          - omegaearth[1] * ppef[0]);

  // ECEF velocty vector
  velecef->i = pm[0][0] * vpef[0] + pm[1][0] * vpef[1] + pm[2][0] * vpef[2];
  velecef->j = pm[0][1] * vpef[0] + pm[1][1] * vpef[1] + pm[2][1] * vpef[2];
  velecef->k = pm[0][2] * vpef[0] + pm[1][2] * vpef[1] + pm[2][2] * vpef[2];
}

/*
 * Transform ECEF position to geodetic latitude, longitude, and altitude
 *
 * Inputs:  posecef   - Position vector in TEME frame, km
 *          velteme   - Velocity vector in TEME frame, km/s
 *          julian    - Julian time of interest
 *          maxiter   - Maximum iteration count for geodetic latitude
 *          tolerance - Desired precision threshold for geodetic latitude
 * Outputs: latlonalt - Position vector in ECEF frame, rad
 * Returns: None
 */
void
ecef2latlonalt
(
    vec3* posecef,
    double julian,
    unsigned int maxiter,
    double tolerance,
    vec3* latlonalt
)
{
  if ((tolerance >= 10.0) || (tolerance < 0.0))
  {
    tolerance = 1.0e-6; // Default tolerance
  }

  // Longitude
  double ijsq = sqrt(posecef->i * posecef->i + posecef->j * posecef->j);

  if (fabs(ijsq) < tolerance)
  {
    latlonalt->lon= signof(posecef->k) * PIDIV2;
  }
  else
  {
    latlonalt->lon= atan2( posecef->j, posecef->i );
  }

  // Wrap around
  if (fabs(latlonalt->lon) >= PI)
  {
    if (latlonalt->lon < 0.0)
    {
      latlonalt->lon += TWOPI;
    }
    else
    {
      latlonalt->lon -= TWOPI;
    }
  }

  // Latitude
  double posmag = magvec3(posecef);
  latlonalt->lat = asin(posecef->k / posmag);

  // Converge latitude to the goid over 10 iterations or less
  const double eesqrd = 0.006694385000; // Earth eccentricity squared
  double c, latsine;
  int i = 1;
  double delta = latlonalt->lat + 10.0;

  while ((fabs(delta - latlonalt->lat) >= tolerance) && (i < maxiter))
  {
    delta   = latlonalt->lat;
    latsine = sin(latlonalt->lat);
    c       = RE / (sqrt(1.0 - eesqrd * latsine * latsine));
    latlonalt->lat = atan((posecef->k + c * eesqrd * latsine) / ijsq);
    i++;
  }

  // Altitude
  if ((PIDIV2 - fabs(latlonalt->lat)) > DEG2RAD)
  {
    latlonalt->alt = (ijsq / cos(latlonalt->lat)) - c;
  }
  else
  {
    latlonalt->alt = posecef->k / sin(latlonalt->lat) - c * (1.0 - eesqrd);
  }
}

/*
 * Transform geodetic latitude, longitude, and altitude to ECEF position vector
 *
 * Inputs:  latlonalt - Geodetic latitude, longitude and altitude vector
 * Outputs: posecef   - Position vector in ECEF frame, rad
 * Returns: None
 */
void
latlonalt2ecef
(
    vec3* latlonalt,
    vec3* posecef
)
{
  // Geocentric latitude
  double gclat = atan(pow(1 - FLATT, 2) * tan(latlonalt->lat));

  // Radius of Earth ad surface point
  double rsurf = sqrt(pow(RE, 2) / ((1 / pow(1.0 - FLATT, 2) - 1) *
                 pow(sin(gclat), 2) + 1));

  // ECEF vector
  posecef->i = rsurf * cos(gclat) * cos(latlonalt->lon)
             + latlonalt->alt * cos(latlonalt->lat) * cos(latlonalt->lon);
  posecef->j = rsurf * cos(gclat) * sin(latlonalt->lon)
             + latlonalt->alt * cos(latlonalt->lat) * sin(latlonalt->lon);
  posecef->k = rsurf * sin(gclat) + latlonalt->alt * sin (latlonalt->lat);
}

/*------------------------------------------------------------------------------
*
*                           procedure rv_tradec
*
*  this procedure converts topocentric right-ascension declination with
*    position and velocity vectors. uses velocity vector to find the
*    solution of singular cases.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    rijk        - ijk position vector            er
*    vijk        - ijk velocity vector            er/tu
*    rsijk       - ijk site position vector       er
*    direct      -  direction to convert          eFrom  eTo
*
*  outputs       :
*    rho         - top radius of the sat          er
*    trtasc      - top right ascension            rad
*    tdecl       - top declination                rad
*    drho        - top radius of the sat rate     er/tu
*    tdrtasc     - top right ascension rate       rad/tu
*    tddecl      - top declination rate           rad/tu
*
*  locals        :
*    rhov        - ijk range vector from site     er
*    drhov       - ijk velocity vector from site  er / tu
*    temp        - temporary extended value
*    temp1       - temporary extended value
*    i           - index
*
*  coupling      :
*    astMath::mag         - astMath::magnitude of a vector
*    atan2       - arc tangent function that resolves the quadrant ambiguities
*    arcsin      - arc sine function
*    lncom2      - linear combination of 2 vectors
*    addvec      - add two vectors
*    dot         - dot product of two vectors
*
*  references    :
*    vallado       2013, 260, alg 26
-----------------------------------------------------------------------------*/

void rv_tradec
(
vec3* rijk, vec3* vijk, vec3* rsijk,
int direct,
double* rho, double* trtasc, double* tdecl,
double* drho, double* dtrtasc, double* dtdecl
)
{
  const double small = 0.00000001;
  const double omegaearth = 0.05883359221938136;  // earth rot rad/tu

  vec3 earthrate, rhov, drhov, vsijk;
  double   latgc, temp, temp1;

  latgc = asin(rsijk->k / magvec3(rsijk));
  earthrate.x = 0.0;
  earthrate.y = 0.0;
  earthrate.z = omegaearth;
  crossvec3(&earthrate, rsijk, &vsijk);

/*  if (direct == 1) //from
  {
    // --------  calculate topocentric vectors ------------------
    rhov.x = (*rho * cos(*tdecl) * cos(*trtasc));
    rhov.y = (*rho * cos(*tdecl) * sin(*trtasc));
    rhov.z = (*rho * sin(*tdecl));

    drhov.x = (*drho * cos(*tdecl) * cos(*trtasc) -
      *rho * sin(*tdecl) * cos(*trtasc) * *dtdecl -
      *rho * cos(*tdecl) * sin(*trtasc) * *dtrtasc);
    drhov.y = (*drho * cos(*tdecl) * sin(*trtasc) -
      *rho * sin(*tdecl) * sin(*trtasc) * *dtdecl +
      *rho * cos(*tdecl) * cos(*trtasc) * *dtrtasc);
    drhov.z = (*drho * sin(*tdecl) + *rho * cos(*tdecl) * *dtdecl);

    // ------ find ijk range vector from site to satellite ------
    addvect(1.0, &rhov, 1.0, rsijk, rijk);
    addvect(1.0, &drhov, cos(latgc), &vsijk, vijk);
  }
  else //to
  {*/
    /* ------ find ijk range vector from site to satellite ------ */
    addvec3(1.0, rijk, -1.0, rsijk, &rhov);
    addvec3(1.0, vijk, -cos(latgc), &vsijk, &drhov);

    /* -------- calculate topocentric angle and rate values ----- */
    *rho = magvec3(&rhov);
    temp = sqrt(rhov.x * rhov.x + rhov.y * rhov.y);
    if (temp < small)
    {
      temp1 = sqrt(drhov.x * drhov.x + drhov.y * drhov.y);
      *trtasc = atan2(drhov.y / temp1, drhov.x / temp1);
    }
    else
      *trtasc = atan2(rhov.y / temp, rhov.x / temp);

    *tdecl = asin(rhov.z / magvec3(&rhov));

    temp1 = -rhov.y * rhov.y - rhov.x * rhov.x;
    *drho = dotvec3(&rhov, &drhov) / *rho;
    if (fabs(temp1) > small)
      *dtrtasc = (drhov.x * rhov.y - drhov.y * rhov.x) / temp1;
    else
      *dtrtasc = 0.0;
    if (fabs(temp) > small)
      *dtdecl = (drhov.z - *drho * sin(*tdecl)) / temp;
    else
      *dtdecl = 0.0;
//  }
}

double
ecef2range(vec3* obsposecef, vec3* satposecef)
{
  // Observer to satellite vector
  vec3 obs2sat;
  addvec3(1.0, satposecef, -1.0, obsposecef, &obs2sat);

  return magvec3(&obs2sat);
}

void ecef2azel
(
    vec3* posecef,
    vec3* vecef,
    vec3* obsecef,
    double lat,
    double lon,
    double* range,
    double* az,
    double* el,
    double* rrate,
    double* azrate,
    double* elrate
)
{

}

