/*
 * libsgp4ansi.h - an ANSI C-11 SGP4/SDP4 implementation library for sgp4ansid.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */


#ifndef LIBSGP4ANSI_H_
#define LIBSGP4ANSI_H_

#include "types.h"

// SGP4/SDP4 math engine initialization
extern int
orbit_init(orbit*);

// Get position and velocity vectors of the satellite at given time in TEME frame
extern int
posvel_teme(orbit*, time_t*, unsigned int, unsigned int, double, vect*, vect*);

// Get position and velocity vectors of the satellite at given time in ECEF frame
extern int
posvel_ecef(orbit*, time_t*, unsigned int, unsigned int, double, vect*, vect*);

// TODO: Make private!
void teme2ecef(vect*, vect*, double, vect*, vect*);
double unix2jul(time_t*, unsigned int);

#endif /* LIBSGP4ANSI_H_ */
