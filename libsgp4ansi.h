/*
 * libsgp4ansi.h - an ANSI C-11 SGP4/SDP4 implementation library for sgp4ansid.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#ifndef LIBSGP4ANSI_H_
#define LIBSGP4ANSI_H_

#include <stdbool.h>
#include <stdint.h>
#include <time.h>

#include "const.h"

// ************************************************************************* //
//                                VERSION                                    //
// ************************************************************************* //

#define VERSION "0.2"
extern const char libsgp4ansi_version[];

// ************************************************************************* //
//                               CUSTOM TYPES                                //
// ************************************************************************* //

typedef struct _vec3
{
  union {
    double a, i, l, u, x, lat;
  };
  union {
    double b, j, m, v, y, lon;
  };
  union {
    double c, k, n, w, z, alt;
  };
} vec3;

/*
 * Satellite orbital element set
 */
typedef struct _sat
{
  // NORAD TLE
  char         name[25];          // Satellite name, 24 chars + \0
  char         sec_class;         // Security classification, 1 char
  char         int_designator[9]; // International designator, 8 chars + \0
  time_t       epoch;             // Epoch of the TLE
  float        epoch_ms;          // Fractional seconds portion of epoch, ms
  double       julian_epoch;      // Julian time at epoch
  double       mean_motion_dt2;   // 1st deriv. of mean motion div2, rev/day2
  double       mean_motion_ddt6;  // 2nd deriv. of mean motion div6, rev/day3
  double       Bstar;             // Pseudo-ballistic drag coefficient, 1/AE
  double       inclination;       // Orbital inclination, [0;180) deg
  double       right_asc_node;    // Right ascension of asc. node, [0;360) deg
  double       eccentricity;      // Orbital eccentricity, [0;1)
  double       argument_perigee;  // Argument of perigee, [0;360) deg
  double       mean_anomaly;      // Mean anomaly at epoch, [0;360) deg
  double       mean_motion;       // Mean motion at epoch, rev/day
  unsigned int norad_number;      // Catalogue number
  unsigned int orbit_number;      // Number of revolutions at epoch
  // Flags
  bool is_deep_space, use_simple_model, is_24h_resonant, is_12h_resonant;
  // Standard orbital elements
  double GSTo;                     // Greenwich Sidereal Time at epoch
  double xnodp;                   // Original mean motion recovered from TLE
  double aodp;                    // Semimajor axis, AE
  double perigee;                 // Perigee, AE
  double perigee_alt;             // Altitude of perigee from surface, km
  double period;                  // Orbital period
  // Common constants
  double aycof, C1, C4, con41,  eta, omgdot, t2cof, x1mth2, x1m5th2, x7thm1,
         xlcof, xnodcf, xnodot, xmdot;
  // Near space constants
  double C5, D2, D3, D4, delmo, omgcof, sinmo, t3cof, t4cof, t5cof, xmcof;
  // Deep space solar terms
  double se2, se3, si2, si3, sl2, sl3, sl4, sgh2, sgh3, sgh4, sh2, sh3;
  double zmos;
  // Deep space lunar terms
  double ee2, e3,  xi2, xi3, xl2, xl3, xl4, xgh2, xgh3, xgh4, xh2, xh3;
  double zmol;
  // Deep space lunar solar terms
  double peo, pinco, plo, pgho, pho;
  // Deep space long period last perturbed elements
  double inclination_lp;
  double eccentricity_lp;
  double right_asc_node_lp;
  double argument_perigee_lp;
  double mean_anomaly_lp;
  // Deep space resonance terms
  double dedt, didt, dmdt, dndt, dnodt, domdt;
  double xlamo, xfact;
  double d2201, d2211, d3210, d3222 , d4410, d4422, d5220, d5232, d5421, d5433;
  double del1, del2, del3; // TODO: Rename
  // Deep space integrator terms
  double xli, xni, atime;
} sat;

/*
 * Satellite observational data at given time from a given location
 */
typedef struct _obs
{
  char   name[25];  // Satellite name, 24 chars + \0
  vec3   latlonalt; // Satellite projected geodetic coordinates
  double velocity;  // Ground velocity
  double azimuth;   // Azimuth from north
  double elevation; // Elevation over horizon
  double az_rate;   // Instantaneous azimuth rate
  double el_rate;   // Instantaneous elevation rate
  double range;     // Direct range, km
  double rng_rate;  // Range rate, km/s
  bool   is_illum;  // Is the satellite illuminated by the Sun?
} obs;

// ************************************************************************* //
//                                    API                                    //
// ************************************************************************* //

// Initialize SGP4/SDP4 orbit model from raw NORAD TLE lines
extern int
sat_load_tle
(
  const char* tlestr0,
  const char* tlestr1,
  const char* tlestr2,
  sat* s
);

// Expand SGP4/SDP4 orbit elements from an orbit containing NORAD TLE portion
extern int
sat_init(sat*);

// Get position and velocity vectors in the TEME frame at given time since epoch
extern int
sat_propagate(sat*, double, unsigned int, double, vec3*, vec3*);

extern int
sat_observe
(
        sat*    s,
  const time_t* time,
        double  time_ms,
  const vec3*   obs_lla,
        obs*    response
);

// Get classical orbital elements from TEME vectors
extern int
teme2coe
(
  vec3*, vec3*, double*, double*, double*, double*, double*, double*, double*,
  double*, double*, double*, double*
);

#endif /* LIBSGP4ANSI_H_ */
