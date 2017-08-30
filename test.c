/*
 * main.c
 *
 *  Created on: 28 ���. 2017 �.
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
Lon 176.8886� W
Lat 37.4990� S
Alt (km)  414.447
Azm 184.1�
Elv -70.7�
RA  249.5544�
Decl  -19.2665�
Range (km)  12 467.355
RRt (km/s)  -1.636
Vel (km/s)  7.664
Direction Descending
Eclipse No
MA (phase)  32.1� (23)
TA  32.2�
Orbit # 7 311
Mag (illum) Not visible
Constellation Oph
Downlink Doppler shift 795Hz @ 145.8MHz
Uplink Doppler shift -792Hz @ 145.2MHz
*/

/*
Orbitron simulation for 0-0-0 latlonalt observer at 2017-08-29 15:01:57UTC

1ISS
Lon 2.4533� E
Lat 7.9010� N
Alt (km)  402.443
Azm 17.2�
Elv 19.0�
RA  244.2845�
Decl  64.5851�
Range (km)  1 025.207
RRt (km/s)  -3.928
Vel (km/s)  7.673
Direction Descending
Eclipse No
MA (phase)  330.1� (234)
TA  330.0�
Orbit # 7 316
Mag (illum) 1.8 (16%)
Constellation Dra
Downlink Doppler shift 1910Hz @ 145.8MHz
Uplink Doppler shift -1902Hz @ 145.2MHz
*/
