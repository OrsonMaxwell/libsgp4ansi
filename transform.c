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
  omegaearth[2] = 7.29211514670698e-05 * (1.0  - 0.002/86400.0);

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
 * Inputs:  posecef - Position vector in TEME frame, km
 *          velteme - Velocity vector in TEME frame, km/s
 *          julian  - Julian time of interest
 * Outputs: posteme - Position vector in ECEF frame, km
 *          velteme - Velocity vector in ECEF frame, km/s
 * Returns: None
 */
void
ecef2latlonalt(vect* posecef, double julian, vect* latlonalt)
{
  const double tolerance = 1.0e-12;
  const double eesqrd    = 0.006694385000; // Earth eccentricity squared

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
  double c, latsine;
  int i = 1;
  double delta = latlonalt->lat + 10.0;

  while ((fabs(delta - latlonalt->lat) >= tolerance) && (i < 10))
  {
    delta   = latlonalt->lat;
    latsine = sin(latlonalt->lat);
    c       = Re / (sqrt(1.0 - eesqrd * latsine * latsine));
    latlonalt->lat = atan((posecef->k + c * eesqrd * latsine) / ijsq);
    i++;
  }

  // Altitude
  if ((pidiv2 - fabs(latlonalt->lat)) > pi / 180.0)
  {
    latlonalt->alt = (ijsq / cos(latlonalt->lat)) - c;
  }
  else
  {
    latlonalt->alt = posecef->k / sin(latlonalt->lat) - c * (1.0 - eesqrd);
  }
}
