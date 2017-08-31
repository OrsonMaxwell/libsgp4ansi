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
typedef struct sat
{
  // NORAD TLE portion
  char         name[24];      // Satellite name
  unsigned int number;        // Catalogue number
  char         class;         // Classification
  char         designator[8]; // International designator
  time_t       epoch;         // Epoch of the TLE
  double       nprimediv2;    // First derivative of mean motion div2, rev/day2
  double       ndprimediv6;   // Second derivative of mean motion div6, rev/day3
  double       bstar;         // Pseude-ballistic drag coefficient, 1/Earth r
  uint8_t      ephem_type;    // Ephemeris type
  unsigned int elset_number;  // Current element set number
  double       i;             // Orbital inclination, 0..180deg
  double       alpha;         // Right ascension of ascension node, 0..360deg
  double       e;             // Orbital eccentricity, 0.0..1.0
  double       omega;         // Argument of perigee, 0..360deg
  double       Mo;            // Mean anomaly at epoch, 0..360deg
  double       no;            // Mean motion at epoch, rev/day
  unsigned int rev_number;    // Number of revolutions at epoch
  // Standard orbital elements and aux epoch quantities
  double a;
  double alta;
  double altp;
  double eta;
  double C1;
  double C4;
  double C5;
  double x1mth2;
  double x7thm1;
  double mdot;
  double con41;
  double omegaprime;
  double nodedot;
  double omgcof;
  double xmcof;
  double xlcof;
  double nodecf;
  double t2cof;
  double aycof;
  double delMo;
  double sinMo;

  double gsto;   // GST at epoch, rad
  // Misc flags
  bool isimp;    // ???
} sat;

#endif /* TYPES_H_ */
