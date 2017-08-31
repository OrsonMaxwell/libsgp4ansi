// ************************************************************************* //
//                                CUSTOM TYPES                               //
// ************************************************************************* //

#ifndef TYPES_H_
#define TYPES_H_

#include <time.h>
#include <stdbool.h>

/*
 * Satellite orbital element set
 */
typedef struct _orbit
{
  // NORAD TLE portion
  char         name[24];      // Satellite name
  unsigned int number;        // Catalogue number
  char         class;         // Security classification
  char         designator[8]; // International designator
  time_t       epoch;         // Epoch of the TLE
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
  // Standard orbital elements
  double a,     alta,   altp,    aycof,      C1,
         C4,    C5,     con41,   delMo,      eta,
         mdot,  nodecf, nodedot, omegaprime, omgcof,
         sinMo, t2cof,  x1mth2,  x7thm1,     xlcof,
         xmcof;
  // Deep space perturbations elements
  double e3,  ee2,  peo,  pgho, pho,  pinco, plo, se2,
         se3, sgh2, sgh3, sgh4, sh2,  sh3,   si2, si3,
         sl2, sl3,  sl4,  xgh2, xgh3, xgh4,  xh2, xh3,
         xi2, xi3,  xl2,  xl3,  xl4,  zmol,  zmos;

  double gsto;   // GST at epoch, rad
  // Misc flags
  bool isdeepspace;
  bool isimp;    // ???
} orbit;

#endif /* TYPES_H_ */
