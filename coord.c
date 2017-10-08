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
#include <string.h> // TODO: get rid of this
#include <stdio.h> // TODO: get rid of this

#include "libsgp4ansi.h"
#include "const.h"
#include "coord.h"
#include "epoch.h"
#include "vector.h"

/*
 * Solve Kepler's equation for known true anomaly
 *
 * Inputs:  ecc - Eccentricity
 *          nu  - True anomaly, (-2pi; 2pi) rad
 * Returns: mean anomaly
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
 * Get classical orbital elements from TEME vectors
 *
 * Inputs:  posteme - Position vector in TEME frame, km
 *          velteme - Velocity vector in TEME frame, km/s
 * Outputs: None
 * Returns: coe     - Success - struct containint the elements
 *          NULL    - Decayed satellite
 */
coe
teme2coe
(
  const vec3* posteme,
  const vec3* velteme
)
{
  coe e = {0};
  char orbit_type;

  double tolerance = 1.0e-8;

  double magr = vec3_mag(posteme);
  double magv = vec3_mag(velteme);

  // Find h n and e vectors
  vec3 hbar, nbar, ebar;

  hbar = vec3_cross(posteme, velteme);
  double magh = vec3_mag(&hbar);

  if (magh < tolerance)
  {
    return e;
  }

  nbar.x = -hbar.y;
  nbar.y =  hbar.x;
  nbar.z =  0;

  double magn  = vec3_mag(&nbar);
  double c1    = pow(magv, 2) - GM  / magr;
  double rdotv = vec3_dot(posteme, velteme);

  ebar.x = (c1 * posteme->x - rdotv * velteme->x) / GM;
  ebar.y = (c1 * posteme->y - rdotv * velteme->y) / GM;
  ebar.z = (c1 * posteme->z - rdotv * velteme->z) / GM;

  e.ecc = vec3_mag(&ebar);

  // Find a e and semi-latus rectum
  double sme = (pow(magv, 2) * 0.5) - (GM  / magr);

  if (fabs(sme) > tolerance)
  {
    e.a = -GM / (2 * sme);
  }
  else
  {
    e.a = INFINITY;
  }
  e.p = pow(magh, 2) / GM;

  // Find inclination
  double hk = hbar.z / magh;
  e.incl= acos( hk );

  // Determine type of orbit
  // Elliptical, parabolic, hyperbolic inclined
  //strcpy(orbit_type, "ei");
  orbit_type = 1;
  if (e.ecc < tolerance)
  {
    // Circular equatorial
    if ((e.incl < tolerance) || (fabs(e.incl - PI) < tolerance))
    {
      //strcpy(orbit_type,"ce");
      orbit_type = -2;
    }
    else
    {
      // Circular inclined
      //strcpy(orbit_type,"ci");
      orbit_type = -1;
    }
  }
  else
  {
    // Elliptical, parabolic, hyperbolic equatorial
    if ((e.incl < tolerance) || (fabs(e.incl - PI) < tolerance))
    {
      //strcpy(orbit_type,"ee");
      orbit_type = 2;
    }
  }

  // Find longitude of ascending node
  if (magn > tolerance)
  {
    double temp = nbar.x / magn;

    if (fabs(temp) > 1) // TODO: WTF?
    {
      temp = (temp >= 0)?(1.0):(-1.0);
    }

    e.omega = acos(temp);
    if (nbar.y < 0)
    {
      e.omega = TAU - e.omega;
    }
  }
  else
    e.omega = NAN;

  // Find argument of perigee
  if (orbit_type == 1)
  {
    e.argp = vec3_angle(&nbar, &ebar);
    if (ebar.z < 0)
    {
      e.argp= TAU - e.argp;
    }
  }
  else
    e.argp= NAN;

  // Find true anomaly at epoch
  if (orbit_type > 0)
  {
    e.nu = vec3_angle(&ebar, posteme);
    if ( rdotv < 0)
    {
      e.nu= TAU - e.nu;
    }
  }
  else
    e.nu= NAN;

  // Find argument of latitude - circular inclined
  if (orbit_type == -1)
  {
    e.arglat = vec3_angle(&nbar, posteme);
    if (posteme->z < 0)
    {
      e.arglat = TAU - e.arglat;
    }
    e.m = e.arglat;
  }
  else
    e.arglat = NAN;

  // Find longitude of perigee - elliptical equatorial
  if ((e.ecc > tolerance) && (orbit_type == 2))
  {
    double temp = ebar.x / e.ecc;

    if (fabs(temp) > 1)
    {
      temp = (temp >= 0)?(1):(-1);
    }

    e.lonper = acos( temp );

    if (ebar.y < 0)
    {
      e.lonper = TAU - e.lonper;
    }
    if (e.incl > PIDIV2)
    {
      e.lonper = TAU - e.lonper;
    }
  }
  else
    e.lonper = NAN;

  // Find true longitude - circular equatorial
  if  ((magr > tolerance) && (orbit_type == -2))
  {
    double temp = posteme->x / magr;

    if ( fabs(temp) > 1)
    {
      temp = (temp >= 0)?(1):(-1);
    }

    e.truelon = acos(temp);

    if (posteme->y < 0)
    {
      e.truelon = TAU - e.truelon;
    }
    if (e.incl > PIDIV2)
    {
      e.truelon = TAU - e.truelon;
    }
    e.m = e.truelon;
  }
  else
    e.truelon = NAN;

  // Find mean anomaly for all orbits
  if (orbit_type > 0)
  {
    e.m = kepler_newton(e.ecc, e.nu);
  }

  return e;
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
  double GST = jul2gst(julian);

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
 * Inputs:  posecef   - Position vector in TEME frame, km
 * Outputs: None
 * Returns: geodetic coordinates vector
 */
vec3
ecef2geo
(
  const vec3* posecef
)
{
  vec3 geo;

  unsigned int maxiter = 5;
  double tolerance     = 1.0e-12;

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
  int i = 1;
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
 * Inputs:  latlonalt - Geodetic latitude, longitude and altitude vector
 * Outputs: posecef   - Position vector in ECEF frame, km
 *          velecef   - Velocity vector in ECEF frame, km/s
 * Returns: None
 */
void
geo2ecef
(
    const vec3* latlonalt,
          vec3* posecef,
          vec3* velecef
)
{
  // Geocentric latitude
  double gclat = atan(pow(1 - FLATT, 2) * tan(latlonalt->lat));

  // Radius of Earth at surface point
  double rsurf = sqrt(pow(RE, 2) / ((1 / pow(1.0 - FLATT, 2) - 1) *
                 pow(sin(gclat), 2) + 1));

  // ECEF position vector
  posecef->x = rsurf * cos(gclat) * cos(latlonalt->lon)
             + latlonalt->alt * cos(latlonalt->lat) * cos(latlonalt->lon);
  posecef->y = rsurf * cos(gclat) * sin(latlonalt->lon)
             + latlonalt->alt * cos(latlonalt->lat) * sin(latlonalt->lon);
  posecef->z = rsurf * sin(gclat) + latlonalt->alt * sin (latlonalt->lat);

  // ECEF velocity vector
  double factor = TAU * RPSID / 86400;
  velecef->x = -factor * posecef->y;
  velecef->y = factor * posecef->x;
  velecef->z = 0;
}

/*
 * Find elevation from ECEF vectors
 *
 * Inputs:  op - Observer position vector in ECEF frame
 *          dp - Observer to satellite position vector in ECEF frame
 * Outputs: None
 * Returns: elevation angle, [-pi/2; pi/2] rad
 */
double
ecef2el
(
  const vec3*  op,
  const vec3*  dp
)
{
  // Cosine of elevation
  double cosel = (op->x * dp->x + op->y * dp->y + op->z * dp->z)
               / sqrt((pow(op->x, 2) + pow(op->y, 2) + pow(op->z, 2))
                    * (pow(dp->x, 2) + pow(dp->y, 2) + pow(dp->z, 2)));

  return PIDIV2 - acos(cosel);
}

/*
 * Find azimuth from ECEF vectors
 *
 * Inputs:  op - Observer position vector in ECEF frame
 *          dp - Observer to satellite position vector in ECEF frame
 * Outputs: None
 * Returns: azimuth angle, [0; 2pi] rad
 */
double
ecef2az
(
  const vec3*  op,
  const vec3*  dp
)
{
  double cosaz = (-op->z * op->x * dp->x - op->z * op->y * dp->y
               + (pow(op->x, 2) + pow(op->y, 2)) * dp->z)
          / sqrt((pow(op->x, 2) + pow(op->y, 2))
               * (pow(op->x, 2) + pow(op->y, 2) + pow(op->z, 2))
               * (pow(dp->x, 2) + pow(dp->y, 2) + pow(dp->z, 2)));

  double sinaz = (-op->y * dp->x + op->x * dp->y)
         / sqrt((pow(op->x, 2) + pow(op->y, 2))
              * (pow(dp->x, 2) + pow(dp->y, 2) + pow(dp->z, 2)));

  return PIDIV2 - atan2(cosaz, sinaz);
}
