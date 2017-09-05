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
mag(vect* arg)
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

/*
 * Crossing of two 3D vectors
 *
 * Inputs:  vect1  - 1st vector
 *          vect2  - 2nd vector
 * Outputs: result - Resulting vector
 * Returns: None
 */
void cross
(
    vect* vect1, vect* vect2, vect* result
)
{
  result->x= vect1->y * vect2->z - vect1->z * vect2->y;
  result->y= vect1->z * vect2->x - vect1->x * vect2->z;
  result->z= vect1->x * vect2->y - vect1->y * vect2->x;
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
 * Convert year and fractional day to unix time
 *
 * Inputs:  year   - Year
 *          days   - Decimal day since year start
 * Outputs: result - Unix time
 * Returns: None
 */
void
fractday2unix(unsigned int year, double days, time_t* result)
{
  struct tm res_tm;

  int mon_len[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  int day_of_year = (int)floor(days);

  res_tm.tm_year = 100 + year;

  // Month and day of month
  if ((year % 4) == 0) // Leap year?
    mon_len[1] = 29;

  int i = 1, j = 0;
  while ((day_of_year > j + mon_len[i - 1]) && (i < 12))
  {
    j = j + mon_len[i-1];
    i++;
  }
  res_tm.tm_mon = i;
  res_tm.tm_mday = day_of_year - j;

  // Hours, minutes, and seconds
  double temp;
  temp           = (days - day_of_year) * 24.0;
  res_tm.tm_hour = (int)floor(temp);
  temp           = (temp - res_tm.tm_hour) * 60.0;
  res_tm.tm_min  = (int)floor(temp);
  res_tm.tm_sec  = (temp - res_tm.tm_min) * 60.0;

  // TODO: Make cross-platform
  result = mktime(&res_tm) - timezone;
}

/*
 * Convert unix time to Julian date
 *
 * Inputs:  time - Timestamp in unix time
 *          msec - Fracitonal second part, us
 * Returns: Julian date
 */
double
unix2jul(time_t* time, unsigned int msec)
{
  struct tm* t;
  t = gmtime(time);

  return 367.0 * (t->tm_year + 1900)
  - floor((7 * ((t->tm_year + 1900) + floor((t->tm_mon + 10) / 12.0))) * 0.25)
  + floor(275 * (t->tm_mon + 1) / 9.0 )
  + t->tm_mday + 1721013.5
  + ((((double)t->tm_sec + msec / 1000) / 60.0L + t->tm_min) / 60.0
  + t->tm_hour) / 24.0;
}

/*
 * Convert Julian date to Greenwich Siderial Time
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

  // Wrap around
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
  double posmag = mag(posecef);
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
    vect* latlonalt,
    vect* posecef
)
{
  // Geocentric latitude
  double gclat = atan(pow(1 - flatt, 2) * tan(latlonalt->lat));

  // Radius of Earth ad surface point
  double rsurf = sqrt(pow(Re, 2) / ((1 / pow(1.0 - flatt, 2) - 1) *
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
vect* rijk, vect* vijk, vect* rsijk,
int direct,
double* rho, double* trtasc, double* tdecl,
double* drho, double* dtrtasc, double* dtdecl
)
{
  const double small = 0.00000001;
  const double omegaearth = 0.05883359221938136;  // earth rot rad/tu

  vect earthrate, rhov, drhov, vsijk;
  double   latgc, temp, temp1;

  latgc = asin(rsijk->k / mag(rsijk));
  earthrate.x = 0.0;
  earthrate.y = 0.0;
  earthrate.z = omegaearth;
  cross(&earthrate, rsijk, &vsijk);

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
    addvect(1.0, rijk, -1.0, rsijk, &rhov);
    addvect(1.0, vijk, -cos(latgc), &vsijk, &drhov);

    /* -------- calculate topocentric angle and rate values ----- */
    *rho = mag(&rhov);
    temp = sqrt(rhov.x * rhov.x + rhov.y * rhov.y);
    if (temp < small)
    {
      temp1 = sqrt(drhov.x * drhov.x + drhov.y * drhov.y);
      *trtasc = atan2(drhov.y / temp1, drhov.x / temp1);
    }
    else
      *trtasc = atan2(rhov.y / temp, rhov.x / temp);

    *tdecl = asin(rhov.z / mag(&rhov));

    temp1 = -rhov.y * rhov.y - rhov.x * rhov.x;
    *drho = dot(&rhov, &drhov) / *rho;
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
ecef2range(vect* obsposecef, vect* satposecef)
{
  // Observer to satellite vector
  vect obs2sat;
  addvect(1.0, satposecef, -1.0, obsposecef, &obs2sat);

  return mag(&obs2sat);
}

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

}

