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

/* Test TLE 29.08.2017
ISS (ZARYA)
1 25544U 98067A   17241.20968750  .00016118  00000-0  25119-3 0  9990
2 25544  51.6408  34.2995 0004424 198.2687 159.0461 15.54014642 73102

Epoch is August 29 2017 05:01:57 UTC?

Orbitron simulation for 0-0-0 latlonalt observer at 2017-08-29 06:01:57UTC

1ISS
Lon 176.8886А W
Lat 37.4990А S
Alt (km)  414.447
Azm 184.1А
Elv -70.7А
RA  249.5544А
Decl  -19.2665А
Range (km)  12 467.355
RRt (km/s)  -1.636
Vel (km/s)  7.664
Direction Descending
Eclipse No
MA (phase)  32.1А (23)
TA  32.2А
Orbit # 7 311
Mag (illum) Not visible
Constellation Oph
Downlink Doppler shift 795Hz
Uplink Doppler shift -664Hz
*/


int
azel(void)
{
  gravconsttype grav = wgs72;
  elsetrec satrec;

  double startmfe, stopmfe, deltamin;

  double ro[3];
  double vo[3];

  double tsince = 60.0;

  char longstr1[130] = "1 25544U 98067A   17241.20968750  .00016118  00000-0  25119-3 0  9990";
  char longstr2[130] = "2 25544  51.6408  34.2995 0004424 198.2687 159.0461 15.54014642 73102";

  twoline2rv
  (
      longstr1,
      longstr2,
      'm',  'c', 'i', grav,
      &startmfe,
      &stopmfe,
      &deltamin,
      &satrec
  );

  sgp4(grav, &satrec, tsince, ro,  vo);

  printf("TEME vectors:\nX:\t%fkm\nY:\t%fkm\nZ:\t%fkm\nVX:\t%fkm/s\nVY:\t%fkm/s\nVZ:\t%fkm/s\n",
         ro[0],ro[1],ro[2],vo[0],vo[1],vo[2]);

  return 0;
}
