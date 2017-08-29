/*
 * libtle2tcp.c
 *
 *  Created on: 28 рту. 2017 у.
 *      Author: Orson
 */

#include <stdio.h>
#include <math.h>

#include "libtle2tcp.h"
#include "sgp4unit.h"
#include "sgp4ext.h"
#include "sgp4io.h"

// Test TLE 29.08.2017
//ISS (ZARYA)
//1 25544U 98067A   17241.20968750  .00016118  00000-0  25119-3 0  9990
//2 25544  51.6408  34.2995 0004424 198.2687 159.0461 15.54014642 73102

// Orbitron simulation for 0-0-0 latlonalt ovserver at 2017-08-29 00:00:00UTC
// AZ:208.5 EL:-28.3 RRt: -3.998km/s (Downlink Doppler shift +1944Hz)
int
azel(void)
{
  printf("azel\n");
  gravconsttype grav = wgs72;
  elsetrec satrec;

  double startmfe, stopmfe, deltamin;


  twoline2rv
  (
      "1 25544U 98067A   17241.20968750  .00016118  00000-0  25119-3 0  9990",
      "2 25544  51.6408  34.2995 0004424 198.2687 159.0461 15.54014642 73102",
      'm',  'm', 'i', grav,
      &startmfe,
      &stopmfe,
      &deltamin,
      &satrec
  );
  return 0;
}
