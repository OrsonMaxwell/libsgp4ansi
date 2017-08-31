/*
 * main.c
 *
 *  Created on: 28 рту. 2017 у.
 *      Author: Orson
 */

#include <stdio.h>
#include <string.h>
#include "libtle2tcp.h"

#include "sgp4unit.h"
#include "sgp4ext.h"
#include "sgp4io.h"

void old_test(void);

int
main ()
{
  orbit iss = {0};

  struct tm epoch_tm;

  epoch_tm.tm_year       = 117;
  epoch_tm.tm_mon        = 8;
  epoch_tm.tm_mday       = 29;
  epoch_tm.tm_hour       = 5;
  epoch_tm.tm_min        = 1;
  epoch_tm.tm_sec        = 57;

  strcpy(iss.name, "ISS");
  iss.number         = 25544;
  iss.class          = 'U';
  strcpy(iss.designator, "98067A  ");
  iss.epoch          = mktime(&epoch_tm);
  iss.nprimediv2     = 0.00016118;
  iss.ndprimediv6    = 0;
  iss.Bstar          = 25119e-3;
  iss.ephem_type     = 0;
  iss.elset_number   = 999;
  iss.i              = 51.6408;
  iss.alpha          = 34.2995;
  iss.e              = 0.0004424;
  iss.omega          = 198.2687;
  iss.Mo             = 159.0461;
  iss.no             = 15.54014642;
  iss.rev_number     = 7310;

  orbit_init(&iss);

  old_test();

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


void old_test(void){
  // Parsing TLE and initializing the math constants
  gravconsttype grav = wgs72;
  elsetrec satrec;

  double startmfe, stopmfe, deltamin;

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

  // Calculating propagation for +1 hour since TLE epoch
  double ro[3], recef[3];
  double vo[3], vecef[3];
  double tsince = 60.0; // 2017-08-29 06:01:57UTC
  //double tsince = 600.0; // 2017-08-29 15:01:57UTC
  sgp4(grav, &satrec, tsince, ro,  vo);

  //printf("TEME vectors:\nX:\t\t%fkm\nY:\t\t%fkm\nZ:\t\t%fkm\nVX:\t\t%fkm/s\nVY:\t\t%fkm/s\nVZ:\t\t%fkm/s\n\n",
  //       ro[0],ro[1],ro[2],vo[0],vo[1],vo[2]);

  // Rotating TEME to ECEF coordinates
  double jday1 = 0.0;
  jday(2017,8,29,6,1,57, &jday1);// 2017-08-29 06:01:57UTC
  //jday(2017,8,29,15,1,57, &jday1); // 2017-08-29 15:01:57UTC

  teme2ecef(ro, vo, jday1, recef, vecef);

  //printf("ECEF vectors:\nX:\t\t%fkm\nY:\t\t%fkm\nZ:\t\t%fkm\nVX:\t\t%fkm/s\nVY:\t\t%fkm/s\nVZ:\t\t%fkm/s\n\n",
  //       recef[0],recef[1],recef[2],vecef[0],vecef[1],vecef[2]);

  // Converting ECEF coodrinates to latlonalt
  double latgc, latgd, lon, hellp;

  ijk2ll(recef, jday1, &latgc, &latgd, &lon, &hellp);

  //printf("LATLONALT:\nLat (gd):\t%f\nLat (gc):\t%f\nLon:\t\t%f\nAlt:\t\t%f\n\n",
  //       latgc*180/pi,latgd*180/pi,lon*180/pi,hellp);

  // Calculate range, azimuth, elevation and their respective rates relative to the observer
  double rsecef[3] = {6378.137, 0.0, 0.0}; // observer to sat vector at latlonalt 0,0,0
  double rho, az, el, drho, daz, del;

  rv_razel(recef, vecef, rsecef, latgd, lon, eTo,
           &rho, &az, &el, &drho, &daz, &del);

  //printf("RAZEL:\nRange:\t\t%f\nAzimuth:\t%f\nElevation:\t%f\nRRate:\t\t%f\nAZRate:\t\t%f\nELRate:\t\t%f\n\n",
  //       rho,(az<0)?((180-az*180)/pi):(az*180/pi),el*180/pi,drho, daz*180/pi, del*180/pi);

  // Calculate doppler shift
  double c = 299792.458;
  double f0down = 145800000;
  double f0up = 145200000;
  double downshift = (-drho / c) * f0down;
  double upshift = (drho / c) * f0up;

  //printf("Doppler shifts:\nDownlink:\t%f\nUplink:\t\t%f\n\n", downshift, upshift);
}

