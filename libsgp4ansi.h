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
#include <time.h>

// ************************************************************************* //
//                                VERSION                                    //
// ************************************************************************* //

#define VERSION "0.1"
extern const char library_version[];

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
 * Satellite deep space orbital integrator values
 */
struct intv
{
  double xndot, xnddt, xldot;
};

/*
 * Satellite orbital element set
 */
typedef struct _sat {
  struct tle {
    char   name[25];
    char   int_designator[9];
    char   sec_class;
    time_t epoch;
    float  epoch_ms;
    double epoch_jul;
    double mean_motion_dt2;
    double mean_motion_ddt6;
    double Bstar;
    double inclination;
    double right_asc_node;
    double eccentricity;
    double argument_perigee;
    double mean_anomaly;
    double mean_motion;
    unsigned int norad_number;
    unsigned int orbit_number;
  } tle;
  // Common constants
  struct comm
  {
    double cosio, sinio, eta, t2cof, a3ovk2, x1mth2, x3thm1, x7thm1, aycof,
           xlcof, xnodcf, c1, c4, omgdot, xnodot, xmdot,  xnodp,  aodp,
           perigee, period;
    bool use_simple_model;
    bool is_deep_space;
  } comm;
  // Near space constants
  struct near
  {
    double c5, omgcof, xmcof, delmo, sinmo,
    d2, d3,     d4,    t3cof, t4cof, t5cof;
  } near;
  // Deep space constants
  struct deep
  {
    double gsto, zmol, zmos;
    // Lunar and Solar terms for epoch
    double sse, ssi, ssl, ssg, ssh;
    // Lunar and Solar secular terms
    double se2,  si2,  sl2,  sgh2, sh2,  se3,  si3,  sl3,
    sgh3, sh3,  sl4,  sgh4, ee2,  e3,   xi2,  xi3,
    xl2,  xl3,  xl4,  xgh2, xgh3, xgh4, xh2,  xh3;
    // Lunar and Solar dot terms
    double d2201, d2211, d3210, d3222, d4410, d4422, d5220,
    d5232, d5421, d5433, del1,  del2,  del3;
    // Geopotential resonance (12h)
    bool is_resonant;
    // Geosynchronous resonance (24h)
    bool is_synchronous;
  } deep;
  // Integrator values for epoch
  struct intv intv0;
  // Integrator values for current a_time
  struct intv intvt;
  // Integrator constants
  struct intc
  {
    double xfact, xlamo;
  } intc;
  // Integrator parametres
  struct intp
  {
    double xli, xni, atime;
  } intp;
} sat;

// ************************************************************************* //
//                                    API                                    //
// ************************************************************************* //

// Printout sat struct
void
sat_print(sat*, const char*);

// Initialize SGP4/SDP4 orbit model from a raw NORAD TLE lines
extern int
sat_load_tle(char*, char*, char*, sat*);

// Expand SGP4/SDP4 orbit elements from an orbit containing NORAD TLE portion
extern int
sat_init(sat*);

#endif /* LIBSGP4ANSI_H_ */
