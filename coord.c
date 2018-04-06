/*
 * coord.c - Coordinate transformation routines for libsgp4ansi.
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
#include <stdlib.h>

#include "libsgp4ansi.h"
#include "const.h"
#include "coord.h"
#include "epoch.h"
#include "vector.h"
#include "solar.h"

/*
 * Solve Kepler's equation for known true anomaly
 *
 * Inputs:  ecc - Eccentricity
 *          nu  - True anomaly, [-2pi; 2pi) rad
 * Returns: m   - mean anomaly
 */
double
kepler_newton
(
  const double ecc,
  const double nu
)
{
  double m = INFINITY;
  double small = 1.0e-8;

  double sine, cose, e0;

  // Circular
  if (fabs(ecc) < small)
  {
    m = nu;
  }
  else
  {
    // Elliptical
    if (ecc < 1 - small)
    {
      sine = (sqrt(1 - pow(ecc, 2)) * sin(nu)) / (1 + ecc * cos(nu));
      cose = (ecc + cos(nu)) / (1 + ecc * cos(nu));
      e0   = atan2(sine, cose);
      m   = e0 - ecc * sin(e0);
    }
    else
    {
      // Hyperbolic
      if (ecc > 1 + small)
      {
        if ((ecc > 1) && (fabs(nu) + 0.00001 < PI - acos(1 / ecc)))
        {
          sine = (sqrt(pow(ecc, 2) - 1) * sin(nu)) / (1 + ecc * cos(nu) );
          e0   = log(sine + sqrt(pow(sine, 2) + 1));
          m   = ecc * sinh(e0) - e0;
        }
      }
      else
      {
        // Parabolic
        if (fabs(nu) < 168 * PI / 180) // Arbitrary
        {
          e0 = tan(nu * 0.5);
          m = e0 + pow(e0, 3) / 3;
        }
      }
    }
  }

  if (ecc < 1)
  {
    m = fmod(m, TAU);

    if (m < 0)
    {
      m += TAU;
    }
  }

  return m;
}

/*
 * Transform position and velocity vectors from TEME to ECEF frame of reference
 *
 * Inputs:  posteme - Position vector in TEME frame, km
 *          velteme - Velocity vector in TEME frame, km/s
 *          julian  - Julian time of interest
 * Outputs: posecef - Position vector in ECEF frame, km
 *          velecef - Velocity vector in ECEF frame, km/s
 * Returns: None
 */
void
teme2ecef
(
  const vec3*  posteme,
  const vec3*  velteme,
        double julian,
        vec3*  posecef,
        vec3*  velecef
)
{
  // Greenwich Siderial Time, rad
  double GST = jul2gmst(julian);

  // Pef - tod matrix
  double pef_tod[3][3] =
  {
    {cos(GST), -sin(GST),  0},
    {sin(GST),  cos(GST),  0},
    {0,         0,         1}
  };

  // Pseudo Earth fixed position vector
  double ppef[3];
  ppef[0] = pef_tod[0][0] * posteme->x + pef_tod[1][0] * posteme->y
          + pef_tod[2][0] * posteme->z;
  ppef[1] = pef_tod[0][1] * posteme->x + pef_tod[1][1] * posteme->y
          + pef_tod[2][1] * posteme->z;
  ppef[2] = pef_tod[0][2] * posteme->x + pef_tod[1][2] * posteme->y
          + pef_tod[2][2] * posteme->z;

  // Refer to IERS Bulletin - A (Vol. XXVIII No. 030)
  double MJD = julian - 2400000.5;
  double A   = TAU * (MJD - 57226) / 365.25;
  double C   = TAU * (MJD - 57226) / 435;

  // Polar motion coefficients, rad
  double xp;
  double yp;
  xp = (0.1033 + 0.0494*cos(A) + 0.0482*sin(A) + 0.0297*cos(C)
     + 0.0307*sin(C)) * 4.84813681e-6;
  yp = (0.3498 + 0.0441*cos(A) - 0.0393*sin(A) + 0.0307*cos(C)
     - 0.0297*sin(C)) * 4.84813681e-6;

  double pm[3][3] =
  {
    {cos(xp),             0,        -sin(xp)},
    {sin(xp) * sin(yp),   cos(yp),   cos(xp) * sin(yp)},
    {sin(xp) * cos(yp),  -sin(yp),   cos(xp) * cos(yp)}
  };

  // ECEF postion vector
  double omegaearth[3];
  posecef->i = pm[0][0] * ppef[0] + pm[1][0] * ppef[1] + pm[2][0] * ppef[2];
  posecef->j = pm[0][1] * ppef[0] + pm[1][1] * ppef[1] + pm[2][1] * ppef[2];
  posecef->k = pm[0][2] * ppef[0] + pm[1][2] * ppef[1] + pm[2][2] * ppef[2];

  // Earth angular rotation vector
  omegaearth[0] = 0;
  omegaearth[1] = 0;
  omegaearth[2] = OMEGAE * (1  - 0.002 / 86400.0);

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
 * Outputs: None
 * Returns: geo     - Geodetic coordinates vector
 */
vec3
ecef2geo
(
  const vec3* posecef
)
{
  vec3 geo;

  unsigned int maxiter   = 5;
  double       tolerance = 1.0e-12;

  // Longitude
  double ijsq = sqrt(posecef->i * posecef->i + posecef->j * posecef->j);

  if (fabs(ijsq) < tolerance)
  {
    geo.lon = ((posecef->k < 0)?-1:1) * PIDIV2;
  }
  else
  {
    geo.lon = atan2(posecef->j, posecef->i);
  }

  // Wrap around
  if (fabs(geo.lon) >= PI)
  {
    if (geo.lon < 0)
    {
      geo.lon += TAU;
    }
    else
    {
      geo.lon -= TAU;
    }
  }

  // Latitude
  double posmag  = vec3_mag(posecef);
         geo.lat = asin(posecef->k / posmag);

  // Converge latitude to the goid over 10 iterations or less
  double c, latsine;
  int    i     = 1;
  double delta = geo.lat + 10;

  while ((fabs(delta - geo.lat) >= tolerance) && (i < maxiter))
  {
    delta   = geo.lat;
    latsine = sin(geo.lat);
    c       = RE / (sqrt(1 - ECC * ECC * latsine * latsine));
    geo.lat = atan((posecef->k + c * ECC * ECC * latsine) / ijsq);
    i++;
  }

  // Altitude
  if ((PIDIV2 - fabs(geo.lat)) > DEG2RAD)
  {
    geo.alt = (ijsq / cos(geo.lat)) - c;
  }
  else
  {
    geo.alt = posecef->k / sin(geo.lat) - c * (1 - ECC * ECC);
  }

  return geo;
}

/*
 * Transform geodetic latitude, longitude, and altitude to ECEF position vector
 *
 * Inputs:  geo  - Geodetic latitude-longitude-altitude vector, rad, rad, km
 * Outputs: None
 * Returns: ecef - Position vector in ECEF frame, km
 */
vec3
geo2ecef
(
    const vec3* geo
)
{
  vec3 ecef;

  // Geocentric latitude
  double gclat = atan(pow(1 - FLATT, 2) * tan(geo->lat));

  // Radius of Earth at surface point
  double rsurf = sqrt(pow(RE, 2) / ((1 / pow(1.0 - FLATT, 2) - 1) *
                 pow(sin(gclat), 2) + 1));

  // ECEF position vector
  ecef.x = rsurf * cos(gclat) * cos(geo->lon)
             + geo->alt * cos(geo->lat) * cos(geo->lon);
  ecef.y = rsurf * cos(gclat) * sin(geo->lon)
             + geo->alt * cos(geo->lat) * sin(geo->lon);
  ecef.z = rsurf * sin(gclat) + geo->alt * sin (geo->lat);

  return ecef;
}

/*
 * Find azimuth, elevation and range from ECEF vectors
 *
 * Inputs:  op      - Observer position vector in ECEF frame
 *          dp      - Observer to satellite position vector in ECEF frame
 * Outputs: None
 * Returns: azelrng - Azimuth, elevation, range vector, rad, rad, km
 */
vec3
ecef2azelrng
(
  const vec3* op,
  const vec3* dp
)
{
  vec3 azelrng;

  // Find azimuth
  double cosaz = (-op->z * op->x * dp->x - op->z * op->y * dp->y
               + (pow(op->x, 2) + pow(op->y, 2)) * dp->z)
          / sqrt((pow(op->x, 2) + pow(op->y, 2))
               * (pow(op->x, 2) + pow(op->y, 2) + pow(op->z, 2))
               * (pow(dp->x, 2) + pow(dp->y, 2) + pow(dp->z, 2)));

  double sinaz = (-op->y * dp->x + op->x * dp->y)
         / sqrt((pow(op->x, 2) + pow(op->y, 2))
              * (pow(dp->x, 2) + pow(dp->y, 2) + pow(dp->z, 2)));

  azelrng.az = PIDIV2 - atan2(cosaz, sinaz);
  if (azelrng.az < 0)
  {
    azelrng.az += TAU;
  }

  // Find elevation
  double cosel = (op->x * dp->x + op->y * dp->y + op->z * dp->z)
               / sqrt((pow(op->x, 2) + pow(op->y, 2) + pow(op->z, 2))
                    * (pow(dp->x, 2) + pow(dp->y, 2) + pow(dp->z, 2)));

  azelrng.el = PIDIV2 - acos(cosel);

  // Find range
  azelrng.rng = vec3_mag(dp);

  return azelrng;
}

/*
 * Convert equatorial vector to az-el-rng vector from ground st-n at given time
 *
 * Inputs:  radecrv   - Equatorial coordinates vector
 *          obs_geo   - Geodetic coordinates of the ground station, rad, rad, km
 *          timestamp - Unix time
 *          time_ms   - Millisecond portion if the time
 * Outputs: None
 * Returns: azelrng   - Azimuth-elevation-range vector, rad, rad, km
 */
vec3
eq2azelrng
(
  const vec3*  radecrv,
  const vec3*  obs_geo,
        time_t timestamp,
        float  time_ms
)
{
  vec3 azelrng, tc_radecrv;

//  vec3 *obs_geo, *radecrv;
//  obs_geo      = malloc(sizeof(vec3));
//  obs_geo->lat = 33.35611 * DEG2RAD;
//  obs_geo->lon = -116.8625 * DEG2RAD;
//  obs_geo->alt = 1.706;
//  radecrv      = malloc(sizeof(vec3));
//  radecrv->ra  = 339.530208 * DEG2RAD;
//  radecrv->dec = -15.771083 * DEG2RAD;
//  radecrv->rv  = 0.37276 * AU;
//  timestamp    = 1062040620;
//  time_ms      = 0;

  // Julian date
  double julian_date = unix2jul(timestamp, time_ms);

  // find Earth nutation correction
  double D, M, Mdot, F, Omega, Ldot, dpsi, depsilon;
  nutation((julian_date - J2000) / 36525, &D, &M, &Mdot, &F, &Omega, &Ldot,
           &dpsi, &depsilon);

  // Sidereal time
  double GMST = jul2gmst(julian_date);
  double AGST = GMST - dpsi * cos(depsilon);
  // Hour angle
  double H    = fmod(AGST - radecrv->ra + obs_geo->lon, TAU);

//  printf("GMST:%lf\n", GMST * RAD2DEG);
//  printf("AGST:%lf\n", AGST * RAD2DEG);
//  printf("a:   %lf\n", radecrv->ra * RAD2DEG);
//  printf("L:   %lf\n", obs_geo->lon * RAD2DEG);
//  printf("H:   %lf\n", H * RAD2DEG + 360);
  // ECEF observer vector
  vec3 obsposecef = geo2ecef(obs_geo);

  // Switching to topocentric equatorial frame
  // Geocentric latitude
  double gclat  = atan(pow(1 - FLATT, 2) * tan(obs_geo->lat));
  // Mean equatorial parallax
  double pi       = asin(RE / radecrv->rv);
  // Parallax compensation
  double ro       = vec3_mag(&obsposecef) / RE;
  double dalpha   = atan2(-ro * cos(gclat) * sin(pi) * sin(H),
                          cos(radecrv->dec) - ro * cos(gclat) * sin(pi)
                          * sin(H));

//  printf("pi:  %lf\n", pi * RAD2DEG);
//  printf("Ro:  %lf\n", vec3_mag(&obsposecef));
//  printf("Gcl: %lf\n", gclat * RAD2DEG);
//  printf("da:  %lf\n", dalpha * RAD2DEG);

  tc_radecrv.ra   = radecrv->ra + dalpha;
  tc_radecrv.dec  = atan2((sin(radecrv->dec) - ro * sin(gclat) * sin(pi))
                          * cos(dalpha),
                           cos(radecrv->dec) - ro * cos(gclat) * sin(pi)
                          * cos(H));

  printf("Ra:      %.6f\n", radecrv->ra * RAD2DEG);
  printf("Dec:     %.6f\n", radecrv->dec * RAD2DEG);
  printf("Ra':     %.6f\n", tc_radecrv.ra * RAD2DEG);
  printf("Dec':    %.6f\n", tc_radecrv.dec * RAD2DEG);

  // Equation 4-12 (Elevation Deg)
  azelrng.el = asin(sin(obs_geo->lat) * sin(tc_radecrv.dec) + cos(obs_geo->lat)
             * cos(tc_radecrv.dec) * cos(H));

  // Equation 4-13 / 4-14 (Adaptation) (Azimuth Deg)
  azelrng.az = fmod(atan2(-sin(H) * cos(tc_radecrv.dec) / cos(azelrng.el),
               (sin(tc_radecrv.dec) - sin(azelrng.el) * sin(obs_geo->lat))
             / (cos(azelrng.el) * cos(obs_geo->lat))), TAU);

  if (azelrng.az < 0)
  {
    azelrng.az += TAU;
  }

  azelrng.rng = radecrv->rv - RE * ro;

//  free(radecrv);
//  free(obs_geo);
  return azelrng;
}

/*
 * Transform position vector from equatorial to TEME frame
 *
 * Inputs:  radecrv - Position vector in equatorial frame
 * Outputs: None
 * Returns: posteme - Position vector in TEME frame, km
 */
vec3
eq2teme
(
  const vec3* radecrv
)
{
  vec3 posteme;

  posteme.x = (radecrv->rv * cos(radecrv->dec) * cos(radecrv->ra));
  posteme.y = (radecrv->rv * cos(radecrv->dec) * sin(radecrv->ra));
  posteme.z = (radecrv->rv * sin(radecrv->dec));

  return posteme;
}

/*
 * Cast a ray vector onto the WGS ellipsoid centered around origin
 *
 * Inputs:  origin - Ray origin vector
 *          dir    - direction unit vector
 * Outputs: None
 * Returns: shadow - point vector of closest intersection (zero when no
 *                   solutions found)
 */
vec3
cast2ellipsoid
(
  vec3* origin,
  vec3* dir
)
{
  vec3 result = {0};

  // For clarity
  double x  = origin->x;
  double y  = origin->y;
  double z  = origin->z;
  double u  = dir->x - origin->x;
  double v  = dir->y - origin->y;
  double w  = dir->z - origin->z;

  double x2 = pow(x, 2);
  double y2 = pow(y, 2);
  double z2 = pow(z, 2);
  double u2 = pow(u, 2);
  double v2 = pow(v, 2);
  double w2 = pow(w, 2);

  // Semiaxes of the ellipsoid, squared
  double a2 = pow(RE, 2);
  double b2 = pow(RE / (1 + FLATT), 2);

  double d  = (b2 * (u2 + v2) + a2 * w2);

  // Closest intersection solution
  double t;
  if (d != 0)
  {
    t = 1 / d * (b2 * (u * x + v * y)
               + a2 * w * z + 0.5 * sqrt(4 * pow((b2 * (u * x + v * y)
               + a2 * w * z), 2) - 4 * d * (b2 * (-a2 + x2 + y2) + a2 * z2)));
    t *= -1;


    result.x = x + u * t;
    result.y = y + v * t;
    result.z = z + w * t;
  }

  // Back substitute for intersection point vector
  return result;
}

