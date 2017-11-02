/*
 * libsgp4ansi.h - an ANSI C-11 SGP4/SDP4 implementation library for sgp4ansid.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 * Fundamentals of Astrodynamics and Applications, D. Vallado, Second Edition
 * Astronomical Algorithms, Jean Meeus
 *
 * Copyright � 2017 Orson J. Maxwell. Please see LICENSE for details.
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

#define LIBSGP4ANSI_VERSION "0.4"
extern const char libsgp4ansi_version[];

// ************************************************************************* //
//                               CUSTOM TYPES                                //
// ************************************************************************* //

/*
 * Enumerated skylight types for convenience
 */
typedef enum _skylight
{
  Nighttime,
  Astronomical,
  Nautical,
  Civil,
  Daytime
} skylight;

typedef struct _vec3
{
  union {
    double a, i, l, u, x, lat, az, ra;
  };
  union {
    double b, j, m, v, y, lon, el, dec;
  };
  union {
    double c, k, n, w, z, alt, rng, rv;
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
  double GSTo;                    // Greenwich Sidereal Time at epoch
  double xnodp;                   // Original mean motion recovered from TLE
  double aodp;                    // Semimajor axis, AE
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
  // Deep space resonance terms
  double dedt, didt, dmdt, dndt, dnodt, domdt;
  double xlamo, xfact;
  double d2201, d2211, d3210, d3222 , d4410, d4422, d5220, d5232, d5421, d5433;
  double del1, del2, del3;
  // Deep space integrator terms
  double xli, xni, atime;
} sat;

/*
 * Satellite classical orbital elements
 */
typedef struct _coe
{
  double p;       // semilatus rectum, km
  double a;       // semimajor axis, km
  double ecc;     // eccentricity
  double incl;    // inclination                    [0; pi)  rad
  double omega;   // longitude of ascending node    [0; 2pi) rad
  double argp;    // argument of perigee            [0; 2pi) rad
  double nu;      // true anomaly                   [0; 2pi) rad
  double m;       // mean anomaly                   [0; 2pi) rad
  double arglat;  // argument of latitude      (ci) [0; 2pi) rad
  double truelon; // true longitude            (ce) [0; 2pi) rad
  double lonper;  // longitude of periapsis    (ee) [0; 2pi) rad
} coe;

/*
 * Satellite observational data at given time from a given location
 */
typedef struct _obs
{
  vec3   latlonalt;    // Satellite projected geodetic coordinates
  vec3   azelrng;      // Azimuth-Elevation-Range vector
  double velocity;     // Satellite velocity, km/s
  double rng_rate;     // Distance change rate, km/s
  bool   is_illum;     // Is the satellite illuminated by the Sun?
  vec3   sun_azelrng;  // Azimuth-Elevation-Range of The Sun
  vec3   moon_azelrng; // Azimuth-Elevation-Range of The Moon
} obs;

/*
 * Satellite pass data at given time from a given location
 */
typedef struct _pass
{
  time_t   aos_t;     // Acquisition of signal unix time
  time_t   tca_t;     // Time if closest approach unix time
  time_t   los_t;     // Loss of signal unix time
  time_t   flare_t;   // Unix time of illumination
  time_t   eclipse_t; // Unix time of (penumbral) eclipse
  vec3     aos;       // Acquisition of signal azimuth-elevation-range vector
  vec3     tca;       // Time if closest approach azimuth-elevation-range vector
  vec3     los;       // Loss of signal azimuth-elevation-range vector
  vec3     flare;     // Illumination azimuth-elevation-range vector
  vec3     eclipse;   // Umbral eclipse azimuth-elevation-range vector
  skylight sky;       // Skylight type
} pass;

// ************************************************************************* //
//                                    API                                    //
// ************************************************************************* //

// Initialize SGP4/SDP4 orbit model from raw NORAD TLE lines
int
sat_load_tle
(
  const char* tlestr0,
  const char* tlestr1,
  const char* tlestr2,
  sat* s
);

// Initialize SGP4/SDP4 orbit model from NORAD parametres
int
sat_load_params
(
  const char         name[25],
        char         sec_class,
  const char         int_designator[9],
        time_t       epoch,
        float        epoch_ms,
        double       mean_motion_dt2,
        double       mean_motion_ddt6,
        double       Bstar,
        double       inclination,
        double       right_asc_node,
        double       eccentricity,
        double       argument_perigee,
        double       mean_anomaly,
        double       mean_motion,
        unsigned int norad_number,
        unsigned int orbit_number,
        sat*         s
);

// Expand SGP4/SDP4 orbit elements from an orbit containing NORAD TLE portion
int
sat_init
(
  sat* s
);

// Get classical orbital elements from TEME vectors TODO: make a wrapper
extern coe
teme2coe
(
  const vec3* posteme,
  const vec3* velteme
);

// Get position and velocity vectors in the TEME frame at given time since epoch
int
sat_propagate
(
  sat*         s,
  double       tdelta,
  unsigned int maxiter,
  double       tolerance,
  vec3*        p,
  vec3*        v
);

// Get observational data about the satellite from ground station
int
sat_observe
(
        sat*    s,
        time_t  timestamp,
        float   time_ms,
  const vec3*   obs_geo,
        obs*    result
);

// Find satellite passes over ground station inside given timeframe
int
sat_find_passes
(
        sat*         s,
  const time_t*      start_time,
  const time_t*      stop_time,
  const vec3*        obs_geo,
        unsigned int delta_t,
        double       horizon,
        pass*        passes
);

#endif /* LIBSGP4ANSI_H_ */
