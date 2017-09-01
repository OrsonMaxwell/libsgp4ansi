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

// SGP4/SDP4 propagation function wrapper
extern int
orbit_prop(orbit*, time_t*, unsigned int, double, vect*, vect*);

#endif /* LIBSGP4ANSI_H_ */
