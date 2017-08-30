/*
 * main.c
 *
 *  Created on: 28 рту. 2017 у.
 *      Author: Orson
 */

#include <stdio.h>
#include "libtle2tcp.h"

int
main ()
{
  tle zarya;
  zarya.name   = "ISS";
  zarya.number = 25544;
  zarya.class  = 'U';
  zarya.designator = "98067A";
  zarya.epoch = 0; // !!!!!!!!!!!!!!!!!!!
  zarya.ndotdiv2 = .00016118;
  zarya.nddotdiv6 = 0;
  zarya.bstar = 25119-3;
  zarya.ephem_type
  zarya.element_set
  zarya.inclination
  zarya.right_asc
  zarya.eccentricity
  zarya.arg_of_perigee
  zarya.mean_anomaly
  zarya.mean_motion
  zarya.rev_number = ;




  azel();
  return 0;
}

/* Test TLE 29.08.2017
ISS (ZARYA)
1 25544U 98067A   17241.20968750  .00016118  00000-0  25119-3 0  9990
2 25544  51.6408  34.2995 0004424 198.2687 159.0461 15.54014642 73102

Epoch is August 29 2017 05:01:57 UTC

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
Downlink Doppler shift 795Hz @ 145.8MHz
Uplink Doppler shift -792Hz @ 145.2MHz
*/

/*
Orbitron simulation for 0-0-0 latlonalt observer at 2017-08-29 15:01:57UTC

1ISS
Lon 2.4533А E
Lat 7.9010А N
Alt (km)  402.443
Azm 17.2А
Elv 19.0А
RA  244.2845А
Decl  64.5851А
Range (km)  1 025.207
RRt (km/s)  -3.928
Vel (km/s)  7.673
Direction Descending
Eclipse No
MA (phase)  330.1А (234)
TA  330.0А
Orbit # 7 316
Mag (illum) 1.8 (16%)
Constellation Dra
Downlink Doppler shift 1910Hz @ 145.8MHz
Uplink Doppler shift -1902Hz @ 145.2MHz
*/
