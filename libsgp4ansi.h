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

// ************************************************************************* //
//                             PRIVATE FUNCTIONS                             //
// ************************************************************************* //

// SGP4/SDP4 propagation function wrapper
int
orbit_prop(orbit*, time_t*, unsigned int, double, vect*, vect*);

// SGP4 propagation function implementation
int
orbit_sgp4(orbit*, double, unsigned int, double, vect*, vect*);

// SDP4 propagation function implementation
int
orbit_sdp4(orbit*, double, unsigned int, double, vect*, vect*);

// ************************************************************************* //
//                                 INTERFACE                                 //
// ************************************************************************* //

// SGP4/SDP4 math engine initialization
extern int
orbit_init(orbit*);

#endif /* LIBSGP4ANSI_H_ */
