/*
 * transform.c - Coordinate and time transformation routines for libsgp4ansi.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#include <math.h>
#include <time.h>

#include "libsgp4ansi.h"
#include "const.h"
#include "transform.h"

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

/*
 * Find magnitude of a 3D vector
 *
 * Inputs:  arg    - Source vector
 * Outputs: None
 * Returns: vector arg magnitude
 */
double
magnitude(vect* arg)
{
  return sqrt(arg->x * arg->x + arg->y * arg->y + arg->z * arg->z);
}

/*
 * Copy a 3D vector
 *
 * Inputs:  from - Source vector
 * Outputs: to   - Destination vector
 * Returns: None
 */
void
copyvect(vect* from, vect* to)
{
  to->x = from->x;
  to->y = from->y;
  to->z = from->z;
}

/*
 * Add two 3D vectors with coefficients
 *
 * Inputs:  c1     - 1st vector coefficient
 *          vect1  - 1st vector
 *          c2     - 2nd vector coefficient
 *          vect2  - 2nd vector
 * Outputs: result - Resulting vector
 * Returns: None
 */
void
addvect(double c1, vect* vect1, double c2, vect* vect2, vect* result)
{
  result->x = c1 * vect1->x + c2 * vect2->x;
  result->y = c1 * vect1->y + c2 * vect2->y;
  result->z = c1 * vect1->z + c2 * vect2->z;
}

/*
 * Dot product of two 3D vectors
 *
 * Inputs:  vect1  - 1st vector
 *          vect2  - 2nd vector
 * Outputs: None
 * Returns: dot product
 */
double
dot(vect* vect1, vect* vect2)
{
  return vect1->x * vect2->x + vect1->y * vect2->y + vect1->z * vect2->z;
}

void rotate2(vect* src, double angle, vect* result)
{
  double tempz = src->z;
  double s     = sin(angle);
  double c     = cos(angle);

  result->x = src->x * c - tempz * s;
  result->y = src->y;
  result->z = src->z * c + src->x * s;
}

void rotate3(vect* src, double angle, vect* result)
{
  double tempy  = src->y;
  double sine   = sin(angle);
  double cosine = cos(angle);

  result->y = src->y * cosine - src->x * sine;
  result->x = src->x * cosine + tempy * sine;
  result->z = src->z;
}

// ************************************************************************* //
//                                   TIME                                    //
// ************************************************************************* //

/*
 * Convert unix time to Julian date
 *
 * Inputs:  time - Timestamp in unix time
 *          usec - Fracitonal second part, us
 * Returns: Julian date
 */
double
unix2jul(time_t* time, unsigned int usec)
{
  struct tm* t;
  t = gmtime(time);

  return 367.0 * (t->tm_year + 1900)
  - floor((7 * ((t->tm_year + 1900) + floor((t->tm_mon + 10) / 12.0))) * 0.25)
  + floor(275 * (t->tm_mon + 1) / 9.0 )
  + t->tm_mday + 1721013.5
  + ((((double)t->tm_sec + usec / 1000) / 60.0L + t->tm_min) / 60.0
  + t->tm_hour) / 24.0;
}

/*
 * Convert Julian date to Greenwich Sidereal Time
 *
 * Inputs:  julian - Julian date
 * Returns: GST time, rad
 */
double
jul2gst(double julian)
{
  double result, tempUT1;

  tempUT1 = (julian - 2451545.0) / 36525.0;
  result = -6.2e-6* tempUT1 * tempUT1 * tempUT1 + 0.093104 * tempUT1 * tempUT1 +
      (876600.0*3600 + 8640184.812866) * tempUT1 + 67310.54841;

  result = fmod(result * deg2rad / 240.0, twopi);

  // Check quadrants
  if (result < 0.0)
    result += twopi;

  return result;
}

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
    vect* posteme,
    vect* velteme,
    double julian,
    vect* posecef,
    vect* velecef
)
{
  // Greenwich Sidereal Time, rad
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
  double A   = 2 * pi * (MJD - 57226) / 365.25;
  double C   = 2 * pi * (MJD - 57226) / 435;

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
    vect* posecef,
    double julian,
    unsigned int maxiter,
    double tolerance,
    vect* latlonalt
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
    latlonalt->lon= signof(posecef->k) * pidiv2;
  }
  else
  {
    latlonalt->lon= atan2( posecef->j, posecef->i );
  }

  if (fabs(latlonalt->lon) >= pi)
  {
    if (latlonalt->lon < 0.0)
    {
      latlonalt->lon += twopi;
    }
    else
    {
      latlonalt->lon -= twopi;
    }
  }

  // Latitude
  double posmag = magnitude(posecef);
  latlonalt->lat = asin(posecef->k / posmag);

  // Converge latitude to a geodetic ellipsoid over 10 iterations or less
  const double eesqrd = 0.006694385000; // Earth eccentricity squared
  double c, latsine;
  int i = 1;
  double delta = latlonalt->lat + 10.0;

  while ((fabs(delta - latlonalt->lat) >= tolerance) && (i < maxiter))
  {
    delta   = latlonalt->lat;
    latsine = sin(latlonalt->lat);
    c       = Re / (sqrt(1.0 - eesqrd * latsine * latsine));
    latlonalt->lat = atan((posecef->k + c * eesqrd * latsine) / ijsq);
    i++;
  }

  // Altitude
  if ((pidiv2 - fabs(latlonalt->lat)) > deg2rad)
  {
    latlonalt->alt = (ijsq / cos(latlonalt->lat)) - c;
  }
  else
  {
    latlonalt->alt = posecef->k / sin(latlonalt->lat) - c * (1.0 - eesqrd);
  }
}

/*------------------------------------------------------------------------------
*
*                           procedure rv_razel
*
*  this procedure converts range, azimuth, and elevation and their rates with
*    the geocentric equatorial (ecef) position and velocity vectors.  notice the
*    value of small as it can affect rate term calculations. uses velocity
*    vector to find the solution of singular cases.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    recef       - ecef position vector           km
*    vecef       - ecef velocity vector           km/s
*    rsecef      - ecef site position vector      km
*    latgd       - geodetic latitude              -pi/2 to pi/2 rad
*    lon         - geodetic longitude             -2pi to pi rad
*    direct      -  direction to convert          eFrom  eTo
*
*  outputs       :
*    rho         - satellite range from site      km
*    az          - azimuth                        0.0 to 2pi rad
*    el          - elevation                      -pi/2 to pi/2 rad
*    drho        - range rate                     km/s
*    daz         - azimuth rate                   rad/s
*    del         - elevation rate                 rad/s
*
*  locals        :
*    rhovecef    - ecef range vector from site    km
*    drhovecef   - ecef velocity vector from site km/s
*    rhosez      - sez range vector from site     km
*    drhosez     - sez velocity vector from site  km
*    tempvec     - temporary vector
*    temp        - temporary extended value
*    temp1       - temporary extended value
*    i           - index
*
*  coupling      :
*    astMath::mag         - astMath::magnitude of a vector
*    addvec      - add two vectors
*    rot3        - rotation about the 3rd axis
*    rot2        - rotation about the 2nd axis
*    atan2       - arc tangent function which also resloves quadrants
*    dot         - dot product of two vectors
*    rvsez_razel - find r and v from site in topocentric horizon (sez) system
*    lncom2      - combine two vectors and constants
*    arcsin      - arc sine function
*    sgn         - returns the sign of a variable
*
*  references    :
*    vallado       2013, 265, alg 27
-----------------------------------------------------------------------------*/

void ecef2azel
(
    vect* posecef,
    vect* vecef,
    vect* obsecef,
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
  const double small = 0.0000001;

  double temp, temp1;
  double drhoecef[3];

  // ECEF range vector from observer to satellite
  vect obs2sat;
  vect rrateecef;
  addvect(1.0, posecef, -1.0, obsecef, &obs2sat);
  copyvect(vecef, &rrateecef);
  *range = magnitude(&obs2sat);

  // ------------ convert to sez for calculations -------------
  vect tempvec, rangesez, rratesez;
  rotate3(&obs2sat, lon, &tempvec);
  rotate2(&tempvec, pidiv2 - lat, &rangesez);

  rotate3(&rrateecef, lon, &tempvec);
  rotate2(&tempvec, pidiv2 - lat, &rratesez);

  // ------------ calculate azimuth and elevation -------------
  temp = sqrt(rangesez.x * rangesez.x + rangesez.y * rangesez.y);
  if (fabs(rangesez.y) < small)
    if (temp < small)
    {
      temp1 = sqrt(rratesez.x * rratesez.x +
                   rratesez.y * rratesez.y);
      *az = atan2(rratesez.y / temp1, -rratesez.x / temp1);
    }
    else
      if (rangesez.x > 0.0)
        *az = pi;
      else
        *az = 0.0;
  else
    *az = atan2(rangesez.y / temp, -rangesez.x / temp);

  if (temp < small)  // directly over the north pole
    *el = signof(rangesez.z) * pidiv2; // +- 90
  else
    *el = asin(rangesez.z / magnitude(&rangesez));

  // ----- calculate range, azimuth and elevation rates -------
  *rrate = dot(&rangesez, &rratesez) / *range;
  if (fabs(temp * temp) > small)
    *azrate = (rratesez.x * rangesez.y - rratesez.y * rangesez.x) /
    (temp * temp);
  else
    *azrate = 0.0;

  if (fabs(temp) > 0.00000001)
    *elrate = (rratesez.z - *rrate * sin(*el)) / temp;
  else
    *elrate = 0.0;
}

