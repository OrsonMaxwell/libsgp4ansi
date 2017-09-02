/*
 * test.c
 *
 *  Created on: 28 рту. 2017 у.
 *      Author: Orson
 */

#include <stdio.h>
#include <string.h>
#include "libsgp4ansi.h"

#include "sgp4unit.h"
#include "sgp4ext.h"
#include "sgp4io.h"

void old_test_iss(double*, double*);
void new_test_iss(double*, double*);

void old_test_navstar53(void);
void new_test_navstar53(void);

int
main ()
{
  double oldr[3];
  double oldv[3];

  double newr[3];
  double newv[3];

  printf("-------------- ISS --------------\n");
  old_test_iss(oldr, oldv);
  new_test_iss(newr, newv);

  //12h non-resonant GPS (ecc < 0.5 ecc)
  //printf("----------- navstar53 -----------\n");
  //old_test_navstar53();
  //new_test_navstar53();
/*
  teme2ecef(ro, vo, jday1, recef, vecef);

    printf("ECEF vectors:\nX:\t\t%fkm\nY:\t\t%fkm\nZ:\t\t%fkm\nVX:\t\t%fkm/s\nVY:\t\t%fkm/s\nVZ:\t\t%fkm/s\n\n",
           recef[0],recef[1],recef[2],vecef[0],vecef[1],vecef[2]);

    // Converting ECEF coodrinates to latlonalt
    double latgc, latgd, lon, hellp;

    ijk2ll(recef, jday1, &latgc, &latgd, &lon, &hellp);

    printf("LATLONALT:\nLat (gd):\t%f\nLat (gc):\t%f\nLon:\t\t%f\nAlt:\t\t%f\n\n",
           latgc*180/pi,latgd*180/pi,lon*180/pi,hellp);

    // Calculate range, azimuth, elevation and their respective rates relative to the observer
    double rsecef[3] = {6378.137, 0.0, 0.0}; // observer to sat vector at latlonalt 0,0,0
    double rho, az, el, drho, daz, del;

    rv_razel(recef, vecef, rsecef, latgd, lon, eTo,
             &rho, &az, &el, &drho, &daz, &del);

    printf("RAZEL:\nRange:\t\t%f\nAzimuth:\t%f\nElevation:\t%f\nRRate:\t\t%f\nAZRate:\t\t%f\nELRate:\t\t%f\n\n",
           rho,(az<0)?((180-az*180)/pi):(az*180/pi),el*180/pi,drho, daz*180/pi, del*180/pi);

    // Calculate doppler shift
    double c = 299792.458;
    double f0down = 145800000;
    double f0up = 145200000;
    double downshift = (-drho / c) * f0down;
    double upshift = (drho / c) * f0up;

    printf("Doppler shifts:\nDownlink:\t%f\nUplink:\t\t%f\n\n", downshift, upshift);

*/

  printf("OLD pos: %f\t\t%f\t\t%f\nOLD vel: %f\t\t%f\t\t%f\n",
         oldr[0], oldr[1], oldr[2], oldv[0], oldv[1], oldv[2]);
  printf("NEW pos: %f\t\t%f\t\t%f\nNEW vel: %f\t\t%f\t\t%f\n",
           newr[0], newr[1], newr[2], newv[0], newv[1], newv[2]);

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
void old_test_iss(double* r, double* v){
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

  sgp4(wgs72, &satrec, 600.0, r, v);
}

void old_test_navstar53(void){
  // Parsing TLE and initializing the math constants
  gravconsttype grav = wgs72;
  elsetrec satrec;

  double startmfe, stopmfe, deltamin;

  char longstr1[130] = "1 28129U 03058A   06175.57071136 -.00000104  00000-0  10000-3 0   459";
  char longstr2[130] = "2 28129  54.7298 324.8098 0048506 266.2640  93.1663  2.00562768 18443";

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
}

void new_test_iss(double* r, double* v)
{
  orbit iss = {0};
  vect pos = {0}, vel = {0};
  struct tm epoch_tm, prop_tm;
  time_t prop_time;

  epoch_tm.tm_year       = 117;
  epoch_tm.tm_mon        = 7;
  epoch_tm.tm_mday       = 29;
  epoch_tm.tm_hour       = 5;
  epoch_tm.tm_min        = 1;
  epoch_tm.tm_sec        = 57;

  prop_tm.tm_year        = 117;
  prop_tm.tm_mon         = 7;
  prop_tm.tm_mday        = 29;
  prop_tm.tm_hour        = 15;
  prop_tm.tm_min         = 1;
  prop_tm.tm_sec         = 57;

  prop_time = mktime(&prop_tm) - timezone;

  strcpy(iss.name, "ISS");
  iss.number         = 25544;
  iss.sec_class      = 'U';
  strcpy(iss.designator, "98067A  ");
  iss.epoch          = mktime(&epoch_tm) - timezone;
  iss.epoch_ms       = 0;
  iss.nprimediv2     = 0.00016118;
  iss.ndprimediv6    = 0;
  iss.Bstar          = 0.25119e-3;
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
  orbit_prop(&iss, &prop_time, 0, 10, 1.0e-12, &pos, &vel);

  r[0] = pos.x;
  r[1] = pos.y;
  r[2] = pos.z;
  v[0] = vel.x;
  v[1] = vel.y;
  v[2] = vel.z;
}

void new_test_navstar53(void)
{
  orbit navstar53 = {0};

  struct tm epoch_tm;

  // 175.57071136
  epoch_tm.tm_year       = 106;
  epoch_tm.tm_mon        = 6;
  epoch_tm.tm_mday       = 24;
  epoch_tm.tm_hour       = 13;
  epoch_tm.tm_min        = 41;
  epoch_tm.tm_sec        = 49;

  strcpy(navstar53.name, "NAVSTAR 53");
  navstar53.number         = 28129;
  navstar53.sec_class          = 'U';
  strcpy(navstar53.designator, "03058A  ");
  navstar53.epoch          = mktime(&epoch_tm);
  navstar53.epoch_ms       = 461504;
  navstar53.nprimediv2     = -0.00000104;
  navstar53.ndprimediv6    = 0;
  navstar53.Bstar          = 1.0e-4;
  navstar53.ephem_type     = 0;
  navstar53.elset_number   = 461504;
  navstar53.i              = 54.7298;
  navstar53.alpha          = 324.8098;
  navstar53.e              = 0.0048506;
  navstar53.omega          = 266.2640;
  navstar53.Mo             = 93.1663;
  navstar53.no             = 2.00562768;
  navstar53.rev_number     = 1844;

  orbit_init(&navstar53);
}
