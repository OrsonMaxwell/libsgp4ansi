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
kepler_newton(double ecc, double nu)
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
    m = fmod(m, TWOPI);

    if (m < 0)
    {
      m += TWOPI;
    }
  }

  return m;
}

/*
 * Get classical orbital elements from TEME vectors
 *
 * Inputs:  posteme - Position vector in TEME frame, km
 *          velteme - Velocity vector in TEME frame, km/s
 * Outputs: p       - semilatus rectum               km
 *          a       - semimajor axis                 km
 *          ecc     - eccentricity
 *          incl    - inclination                    [0; pi)  rad
 *          omega   - longitude of ascending node    [0; 2pi) rad
 *          argp    - argument of perigee            [0; 2pi) rad
 *          nu      - true anomaly                   [0; 2pi) rad
 *          m       - mean anomaly                   [0; 2pi) rad
 *          arglat  - argument of latitude      (ci) [0; 2pi) rad
 *          truelon - true longitude            (ce) [0; 2pi) rad
 *          lonper  - longitude of periapsis    (ee) [0; 2pi) rad
 * Returns: 0       - Success
 *         -1       - Decayed satellite
 */
int
teme2coe
(
    vec3* posteme, vec3 *velteme,
    double* p,    double* a,  double* ecc, double* incl,   double* omega,
    double* argp, double* nu, double* m,   double* arglat, double* truelon,
    double* lonper
)
{
  int i;
  char typeorbit[3]; // TODO: Change to numeric to get rid of strings

  double tolerance = 1.0e-8; // TODO: Const?

  double magr = vec3_mag(posteme);
  double magv = vec3_mag(velteme);

  // Find h n and e vectors
  vec3 hbar, nbar, ebar;

  vec3_cross(posteme, velteme, &hbar);
  double magh = vec3_mag(&hbar);

  if (magh < tolerance)
  {
    return -1;
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

  *ecc = vec3_mag(&ebar);

  // Find *a e and semi-latus rectum
  double sme = (pow(magv, 2) * 0.5) - (GM  / magr);

  if (fabs(sme) > tolerance)
  {
    *a = -GM / (2 * sme);
  }
  else
  {
    *a = INFINITY;
  }
  *p = pow(magh, 2) / GM;

  // Find inclination
  double hk = hbar.z / magh;
  *incl= acos( hk );

  // Determine type of orbit
  // Elliptical, parabolic, hyperbolic inclined
  strcpy(typeorbit, "ei");
  if (*ecc < tolerance)
  {
    // Circular equatorial
    if ((*incl < tolerance) || (fabs(*incl - PI) < tolerance))
    {
      strcpy(typeorbit,"ce");
    }
    else
    {
      // Circular inclined
      strcpy(typeorbit,"ci");
    }
  }
  else
  {
    // Elliptical, parabolic, hyperbolic equatorial
    if ((*incl < tolerance) || (fabs(*incl - PI) < tolerance))
    {
      strcpy(typeorbit,"ee");
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

    *omega = acos(temp);
    if (nbar.y < 0)
    {
      *omega = TWOPI - *omega;
    }
  }
  else
    *omega = NAN;

  // Find argument of perigee
  if (strcmp(typeorbit, "ei") == 0)
  {
    *argp = vec3_angle(&nbar, &ebar);
    if (ebar.z < 0)
    {
      *argp= TWOPI - *argp;
    }
  }
  else
    *argp= NAN;

  // Find true anomaly at epoch
  if (typeorbit[0] == 'e')
  {
    *nu = vec3_angle(&ebar, posteme);
    if ( rdotv < 0)
    {
      *nu= TWOPI - *nu;
    }
  }
  else
    *nu= NAN;

  // Find argument of latitude - circular inclined
  if (strcmp(typeorbit,"ci") == 0)
  {
    *arglat = vec3_angle(&nbar, posteme);
    if (posteme->z < 0)
    {
      *arglat = TWOPI - *arglat;
    }
    *m = *arglat;
  }
  else
    *arglat = NAN;

  // Find longitude of perigee - elliptical equatorial
  if ((*ecc > tolerance) && (strcmp(typeorbit,"ee") == 0))
  {
    double temp = ebar.x / *ecc;

    if (fabs(temp) > 1)
    {
      temp = (temp >= 0)?(1.0):(-1.0); // TODO: Bleh
    }

    *lonper = acos( temp );

    if (ebar.y < 0)
    {
      *lonper = TWOPI - *lonper;
    }
    if (*incl > PIDIV2)
    {
      *lonper = TWOPI - *lonper;
    }
  }
  else
    *lonper = NAN;

  // Find true longitude - circular equatorial
  if  ((magr > tolerance) && (strcmp(typeorbit,"ce") == 0))
  {
    double temp = posteme->x / magr;

    if ( fabs(temp) > 1)
    {
      temp = (temp >= 0)?(1.0):(-1.0); // TODO: Bleh
    }

    *truelon = acos(temp);

    if (posteme->y < 0)
    {
      *truelon = TWOPI - *truelon;
    }
    if (*incl > PIDIV2)
    {
      *truelon = TWOPI - *truelon;
    }
    *m = *truelon;
  }
  else
    *truelon = NAN;

  // Find mean anomaly for all orbits
  if (typeorbit[0] == 'e')
  {
    *m = kepler_newton(*ecc, *nu);
  }

  return 0;
}

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
    latlonalt->lon = ((posecef->k < 0)?-1:1) * PIDIV2;
  }
  else
  {
    latlonalt->lon = atan2( posecef->j, posecef->i );
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
  double posmag  = vec3_mag(posecef);
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
    const vec3* latlonalt,
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

  latgc = asin(rsijk->k / vec3_mag(rsijk));
  earthrate.x = 0.0;
  earthrate.y = 0.0;
  earthrate.z = omegaearth;
  vec3_cross(&earthrate, rsijk, &vsijk);

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
    vec3_add(1.0, rijk, -1.0, rsijk, &rhov);
    vec3_add(1.0, vijk, -cos(latgc), &vsijk, &drhov);

    /* -------- calculate topocentric angle and rate values ----- */
    *rho = vec3_mag(&rhov);
    temp = sqrt(rhov.x * rhov.x + rhov.y * rhov.y);
    if (temp < small)
    {
      temp1 = sqrt(drhov.x * drhov.x + drhov.y * drhov.y);
      *trtasc = atan2(drhov.y / temp1, drhov.x / temp1);
    }
    else
      *trtasc = atan2(rhov.y / temp, rhov.x / temp);

    *tdecl = asin(rhov.z / vec3_mag(&rhov));

    temp1 = -rhov.y * rhov.y - rhov.x * rhov.x;
    *drho = vec3_dot(&rhov, &drhov) / *rho;
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
  vec3_add(1.0, satposecef, -1.0, obsposecef, &obs2sat);

  return vec3_mag(&obs2sat);
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

