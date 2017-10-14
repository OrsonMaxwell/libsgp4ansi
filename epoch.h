/*
 * epoch.h - Time routines for libsgp4ansi.
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

#ifndef EPOCH_H_
#define EPOCH_H_

#include <time.h>

// Convert year and fractional day to unix time
int
fractday2unix
(
  unsigned int year,
  double       days,
  time_t*      unix,
  float*       ms
);

// Convert unix time to Julian date
double
unix2jul
(
  time_t time,
  float  ms
);

// Convert Julian date to Greenwich Sidereal Time
double
jul2gst
(
  double julian
);

#endif /* EPOCH_H_ */
