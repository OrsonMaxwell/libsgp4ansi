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

// ************************************************************************* //
//                               CUSTOM TYPES                                //
// ************************************************************************* //
/*
 * 3D vector
 */
typedef struct _vect
{
  union {
    double i;
    double l;
    double u;
    double x;
    double lat;
  };
  union {
    double j;
    double m;
    double v;
    double y;
    double lon;
  };
  union {
    double k;
    double n;
    double w;
    double z;
    double alt;
  };
} vect;

/*
 * Satellite orbital element set
 */
typedef struct _orbit
{
  // NORAD TLE portion
  char         name[24];      // Satellite name
  unsigned int number;        // Catalogue number
  char         sec_class;     // Security classification
  char         designator[10]; // International designator
  time_t       epoch;         // Epoch of the TLE
  unsigned int epoch_ms;      // Fractional seconds portion of epoch, ms
  double       nprimediv2;    // First derivative of mean motion div2, rev/day2
  double       ndprimediv6;   // Second derivative of mean motion div6, rev/day3
  double       Bstar;         // Pseudo-ballistic drag coefficient, 1/Earth r
  uint8_t      ephem_type;    // Ephemeris type
  unsigned int elset_number;  // Current element set number
  double       i;             // Orbital inclination, 0..180deg
  double       alpha;         // Right ascension of ascension node, 0..360deg
  double       e;             // Orbital eccentricity, 0.0..1.0
  double       omega;         // Argument of perigee, 0..360deg
  double       Mo;            // Mean anomaly at epoch, 0..360deg
  double       no;            // Mean motion at epoch, rev/day
  unsigned int rev_number;    // Number of revolutions at epoch
  // Time
  double julepoch, GSTo;
  // Flags
  bool isdeepspace, islowperigee;
  // Standard orbital terms
  double a, altapoR, altperR, aycof, C1, C4, C5, con41, cosi, d2, d3, d4,
         delMo, eta,   mdot,  nodecf, nodedot, omegaprime, omgcof, sinMo, sini,
         t2cof, t3cof, t4cof, t5cof,  x1mth2,  x7thm1,     xlcof,  xmcof;
  // Deep space terms
  double e3,  ee2,  peo, pgho, pho,  pinco, plo, se2,  se3,  sgh2, sgh3, sgh4,
         sh2, sh3,  si2, si3,  sl2,  sl3,   sl4, xgh2, xgh3, xgh4, xh2,  xh3,
         xi2, xi3,  xl2, xl3,  xl4,  zmol,  zmos;
  // Resonant terms
  double d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433,
         dedt,  didt,  dmdt,  dnodt, domdt, del1,  del2,  del3,  xfact, xlamo,
         xli,   xni;
} orbit;

// ************************************************************************* //
//                                 INTERFACE                                 //
// ************************************************************************* //

// Initialize SGP4/SDP4 orbit model from a raw NORAD TLE lines
extern int
tle2orbit(char*, char*, orbit*);

// Expand SGP4/SDP4 orbit elements from an orbit containing NORAD TLE portion
extern int
orbit_init(orbit*);

// Get position and velocity vectors in the TEME frame at given time
extern int
orbit_prop(orbit*, double, unsigned int, double, vect*, vect*);
//orbit_prop(orbit*, time_t*, unsigned int, unsigned int, double, vect*, vect*);

// TODO: Make private!
extern void
print_orbit(orbit*, char*);

#endif /* LIBSGP4ANSI_H_ */
