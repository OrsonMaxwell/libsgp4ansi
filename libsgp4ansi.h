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
 * Satellite orbital element set
 */
typedef struct _sat
{
  // NORAD TLE portion
  char         name[25];          // Satellite name, 24 chars + \0
  char         sec_class;         // Security classification
  char         int_designator[9]; // International designator, 8 chars + \0
  time_t       epoch;             // Epoch of the TLE
  float        epoch_ms;          // Fractional seconds portion of epoch, ms
  double       mean_motion_dt2;   // 1st derivative of mean motion div2, rev/day2
  double       mean_motion_ddt6;  // 2nd derivative of mean motion div6, rev/day3
  double       Bstar;             // Pseudo-ballistic drag coefficient, 1/Earth r
  unsigned int elset_number;      // Current element set number
  double       inclination;       // Orbital inclination, 0..180deg
  double       right_asc_node;    // Right ascension of ascension node, 0..360deg
  double       eccentricity;      // Orbital eccentricity, 0.0..1.0
  double       argument_perigee;  // Argument of perigee, 0..360deg
  double       mean_anomaly;      // Mean anomaly at epoch, 0..360deg
  double       mean_motion;       // Mean motion at epoch, rev/day
  unsigned int norad_number;      // Catalogue number
  unsigned int orbit_number;      // Number of revolutions at epoch
  // Flags
  bool is_deep_space, use_simple_model, is_resonant;
  // Common terms
  double a3ovk2, aycof, C1, C4, cosio,  eta,   omgdot, perigee, period, sinio,
         t2cof, x1mth2, x3thm1, x7thm1, xlcof, xmdot,  xnodcf,  xnodot;
  // Near space terms
  double a, altapoR, altperR, C5 , cosi, d2, d3, d4,
         delMo, mdot,  nodecf, nodedot, omegaprime, omgcof, sinMo, sini,
         t3cof, t4cof, t5cof, xmcof;
  // Deep space terms
  double julepoch, GSTo;
  double e3,  ee2,  peo, pgho, pho,  pinco, plo, se2,  se3,  sgh2, sgh3, sgh4,
         sh2, sh3,  si2, si3,  sl2,  sl3,   sl4, xgh2, xgh3, xgh4, xh2,  xh3,
         xi2, xi3,  xl2, xl3,  xl4,  zmol,  zmos;
  // Resonant terms
  double d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433,
         dedt,  didt,  dmdt,  dnodt, domdt, del1,  del2,  del3,  xfact, xlamo,
         xli,   xni;
  // TODO: REMOVE THIS! DONE FOR BACKWARDS COMPATIBILITY DURING MATH-OVERHAUL
  int isimp;      // use simple model
  int error;      //
  char method;    // unused
  char operationmode; // unused
  char init;      // unused
  double cc1;     // C1
  double cc4;     // C4
  double cc5;     // C5
  double delmo;   // ?
  double argpdot; // omegaprime
  double sinmao;  // ?
  double t;       // ?
  int    irez;    // ?
  double gsto;    // GSTo
  double atime;   // ?
  double bstar;   // Bstar
  double ecco;    // e
  double argpo;   // omega
  double inclo;   // i
  double mo;      // Mo
  double nodeo;   // alpha
  double alta;    // ?
  double altp;    // ?
} sat;


// ************************************************************************* //
//                                    API                                    //
// ************************************************************************* //

// Initialize SGP4/SDP4 orbit model from a raw NORAD TLE lines
extern int
tle2orbit(char*, char*, char*, sat*);

// Expand SGP4/SDP4 orbit elements from an orbit containing NORAD TLE portion
extern int
orbit_init(sat*);

// Get position and velocity vectors in the TEME frame at given time since epoch
extern int
orbit_prop(sat*, double, unsigned int, double, vec3*, vec3*);

// Get position and velocity vectors in the TEME frame at given unix time
extern int
orbit_at(sat*, time_t*, unsigned int, unsigned int, double, vec3*, vec3*);

/* ----------- local functions - only ever used internally by sgp4 ---------- */
void dpper
     (
       double e3,     double ee2,    double peo,     double pgho,   double pho,
       double pinco,  double plo,    double se2,     double se3,    double sgh2,
       double sgh3,   double sgh4,   double sh2,     double sh3,    double si2,
       double si3,    double sl2,    double sl3,     double sl4,    double t,
       double xgh2,   double xgh3,   double xgh4,    double xh2,    double xh3,
       double xi2,    double xi3,    double xl2,     double xl3,    double xl4,
       double zmol,   double zmos,   double inclo,
       char init,
       double* ep,    double* inclp, double* nodep,  double* argpp, double* mp,
       char opsmode
     );

void dscom
     (
       double epoch,  double ep,     double argpp,   double tc,     double inclp,
       double nodep,  double np,
       double* snodm, double* cnodm, double* sinim,  double* cosim, double* sinomm,
       double* cosomm,double* day,   double* e3,     double* ee2,   double* em,
       double* emsq,  double* gam,   double* peo,    double* pgho,  double* pho,
       double* pinco, double* plo,   double* rtemsq, double* se2,   double* se3,
       double* sgh2,  double* sgh3,  double* sgh4,   double* sh2,   double* sh3,
       double* si2,   double* si3,   double* sl2,    double* sl3,   double* sl4,
       double* s1,    double* s2,    double* s3,     double* s4,    double* s5,
       double* s6,    double* s7,    double* ss1,    double* ss2,   double* ss3,
       double* ss4,   double* ss5,   double* ss6,    double* ss7,   double* sz1,
       double* sz2,   double* sz3,   double* sz11,   double* sz12,  double* sz13,
       double* sz21,  double* sz22,  double* sz23,   double* sz31,  double* sz32,
       double* sz33,  double* xgh2,  double* xgh3,   double* xgh4,  double* xh2,
       double* xh3,   double* xi2,   double* xi3,    double* xl2,   double* xl3,
       double* xl4,   double* nm,    double* z1,     double* z2,    double* z3,
       double* z11,   double* z12,   double* z13,    double* z21,   double* z22,
       double* z23,   double* z31,   double* z32,    double* z33,   double* zmol,
       double* zmos
     );

void dsinit
     (
       int whichconst,
       double cosim,  double emsq,   double argpo,   double s1,     double s2,
       double s3,     double s4,     double s5,      double sinim,  double ss1,
       double ss2,    double ss3,    double ss4,     double ss5,    double sz1,
       double sz3,    double sz11,   double sz13,    double sz21,   double sz23,
       double sz31,   double sz33,   double t,       double tc,     double gsto,
       double mo,     double mdot,   double no,      double nodeo,  double nodedot,
       double xpidot, double z1,     double z3,      double z11,    double z13,
       double z21,    double z23,    double z31,     double z33,    double ecco,
       double eccsq,  double* em,    double* argpm,  double* inclm, double* mm,
       double* nm,    double* nodem,
       int* irez,
       double* atime, double* d2201, double* d2211,  double* d3210, double* d3222,
       double* d4410, double* d4422, double* d5220,  double* d5232, double* d5421,
       double* d5433, double* dedt,  double* didt,   double* dmdt,  double* dndt,
       double* dnodt, double* domdt, double* del1,   double* del2,  double* del3,
       double* xfact, double* xlamo, double* xli,    double* xni
     );

void dspace
     (
       int irez,
       double d2201,  double d2211,  double d3210,   double d3222,  double d4410,
       double d4422,  double d5220,  double d5232,   double d5421,  double d5433,
       double dedt,   double del1,   double del2,    double del3,   double didt,
       double dmdt,   double dnodt,  double domdt,   double argpo,  double argpdot,
       double t,      double tc,     double gsto,    double xfact,  double xlamo,
       double no,
       double* atime, double* em,    double* argpm,  double* inclm, double* xli,
       double* mm,    double* xni,   double* nodem,  double* dndt,  double* nm
     );

void initl
     (
       int satn,      int whichconst,
       double ecco,   double epoch,  double inclo,   double* no,
       char* method,
       double* ainv,  double* ao,    double* con41,  double* con42, double* cosio,
       double* cosio2,double* eccsq, double* omeosq, double* posq,
       double* rp,    double* rteosq,double* sinio , double* gsto, char opsmode
     );

bool sgp4init
     (
       int whichconst, char opsmode,   const int satn,     const double epoch,
       const double xbstar,  const double xecco, const double xargpo,
       const double xinclo,  const double xmo,   const double xno,
       const double xnodeo,  sat* sat
     );

bool sgp4
     (
       int whichconst, sat* sat,  double tsince,
       double r[3],  double v[3]
     );


#endif /* LIBSGP4ANSI_H_ */
