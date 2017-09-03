/*
 * test.c
 *
 *  Created on: 28 рту. 2017 у.
 *      Author: Orson
 */

//ISS (ZARYA)
//1 25544U 98067A   17246.76764087  .00004029  00000-0  68274-4 0  9991
//2 25544  51.6442   6.5972 0004102 223.1972 288.5525 15.54044462 73972

// Orbitron simulation from site 54.9246, 38.0475, 180m
/* -- 2017-09-04 04:45:55
Lon 47.3990А E
Lat 45.1120А N
Alt (km)  412.076
Azm 144.8А
Elv 11.6А
RA  84.1638А
Decl  -17.1536А
Range (km)  1 379.692
RRt (km/s)  -0.236
Vel (km/s)  7.668
Direction Ascending
Eclipse No
MA (phase)  200.0А (142)
TA  200.0А
Orbit # 7 402
Mag (illum) 1.1 (55%)
Constellation Lep
*/

/* -- 2017-09-04 05:45:55
1ISS
Lon 69.9492А W
Lat 44.1932А S
Alt (km)  418.842
Azm 252.6А
Elv -66.3А
RA  284.8039А
Decl  -54.9371А
Range (km)  12 109.854
RRt (km/s)  -2.745
Vel (km/s)  7.660
Direction Ascending
Eclipse Umbral
MA (phase)  73.1А (52)
TA  73.1А
Orbit # 7 403
Mag (illum) Not visible
Constellation Tel
*/

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "sgp4unit.h"
#include "sgp4ext.h"
#include "sgp4io.h"

#include "libsgp4ansi.h"
#include "const.h"
#include "transform.h"

void printdiff(char* caption, double old, double new)
{
  printf("| %-17s | %13.7lf | %13.7lf | %.13lf\% |\n", caption, old, new, (old == 0)?0:fabs((new - old) / old) * 100);
}

int
main ()
{
  vect oldposteme;
  vect oldvelteme;
  vect oldposecef;
  vect oldvelecef;
  vect oldlatlonalt;

  vect newposteme;
  vect newvelteme;
  vect newposecef;
  vect newvelecef;
  vect newlatlonalt;

  //******************************* O-L-D **************************************
  // Parsing TLE and initializing the math constants
  gravconsttype grav = wgs72;
  elsetrec satrec;

  double startmfe, stopmfe, deltamin;

  char longstr1[130] = "1 25544U 98067A   17246.76764087  .00004029  00000-0  68274-4 0  9991";
  char longstr2[130] = "2 25544  51.6442   6.5972 0004102 223.1972 288.5525 15.54044462 73972";

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

  double r[3], v[3];

  sgp4(wgs72, &satrec, 560.516666666, r, v); //9h 20m 31s

  oldposteme.x = r[0];
  oldposteme.y = r[1];
  oldposteme.z = r[2];
  oldvelteme.x = v[0];
  oldvelteme.y = v[1];
  oldvelteme.z = v[2];

  //******************************* N-E-W **************************************
  orbit iss = {0};
  struct tm epoch_tm, prop_tm;
  time_t prop_time;

  epoch_tm.tm_year   = 117;
  epoch_tm.tm_mon    = 8;
  epoch_tm.tm_mday   = 3;
  epoch_tm.tm_hour   = 18;
  epoch_tm.tm_min    = 25;
  epoch_tm.tm_sec    = 24;// 117ms

  prop_tm.tm_year    = 117;
  prop_tm.tm_mon     = 8;
  prop_tm.tm_mday    = 4;
  prop_tm.tm_hour    = 4;
  prop_tm.tm_min     = 45;
  prop_tm.tm_sec     = 55;

  prop_time = mktime(&prop_tm) - timezone;

  strcpy(iss.name, "ISS");
  strcpy(iss.designator, "98067A  ");
  iss.number         = 25544;
  iss.sec_class      = 'U';
  iss.epoch          = mktime(&epoch_tm) - timezone;
  iss.epoch_ms       = 117;
  iss.nprimediv2     = 0.00004029;
  iss.ndprimediv6    = 0;
  iss.Bstar          = 0.68274e-4;
  iss.ephem_type     = 0;
  iss.elset_number   = 999;
  iss.i              = 51.6442;
  iss.alpha          = 6.5972;
  iss.e              = 0.0004102;
  iss.omega          = 223.1972;
  iss.Mo             = 288.5525;
  iss.no             = 15.54044462;
  iss.rev_number     = 7397;

  orbit_init(&iss);
  orbit_prop(&iss, &prop_time, 0, 10, 1.0e-12, &newposteme, &newvelteme);
  // ***************************************************************************

  teme2ecef(&oldposteme, &oldposteme, unix2jul(&prop_time, 0), &oldposecef, &oldvelecef);
  teme2ecef(&newposteme, &newposteme, unix2jul(&prop_time, 0), &newposecef, &newvelecef);

  ecef2latlonalt(&oldposecef, unix2jul(&prop_time, 0), 10, 1.0e-12, &oldlatlonalt);
  ecef2latlonalt(&newposecef, unix2jul(&prop_time, 0), 10, 1.0e-12, &newlatlonalt);

  vect obsposlla = {0};  // observer lla vector at latlonalt 0,0,0

  // Dacha
  obsposlla.lat = 54.9246 * deg2rad;
  obsposlla.lon = 38.0475 * deg2rad;
  obsposlla.alt = 0.180;

  vect obsposecef; // observer ECEF vector at latlonalt 0,0,0
  //obsposecef.x = 6378.137;
  //obsposecef.y = 0;
  //obsposecef.z = 0;

  latlonalt2ecef(&obsposlla, &obsposecef);

  double orho, ortasc, odecl, odrho, odrtasc, oddecl;
  double nrho, nrtasc, ndecl, ndrho, ndrtasc, nddecl;

  rv_tradec(&oldposecef, &oldvelecef, &obsposecef, 0, &orho, &ortasc, &odecl, &odrho, &odrtasc, &oddecl);
  rv_tradec(&newposecef, &newvelecef, &obsposecef, 0, &nrho, &nrtasc, &ndecl, &ndrho, &ndrtasc, &nddecl);

  printf("+-------------------+----- OLD -----+----- NEW -----+------ DIFF ------+\n");
  printdiff("TEME pos x", oldposteme.x, newposteme.x);
  printdiff("TEME pos y", oldposteme.y, newposteme.y);
  printdiff("TEME pos z", oldposteme.z, newposteme.z);
  printdiff("TEME vel x", oldvelteme.x, newvelteme.x);
  printdiff("TEME vel y", oldvelteme.y, newvelteme.y);
  printdiff("TEME vel z", oldvelteme.z, newvelteme.z);
  printdiff("ECEF pos x", oldposecef.x, newposecef.x);
  printdiff("ECEF pos y", oldposecef.y, newposecef.y);
  printdiff("ECEF pos z", oldposecef.z, newposecef.z);
  printdiff("ECEF vel x", oldvelecef.x, newvelecef.x);
  printdiff("ECEF vel y", oldvelecef.y, newvelecef.y);
  printdiff("ECEF vel z", oldvelecef.z, newvelecef.z);
  printdiff("Latitude", oldlatlonalt.lat * rad2deg, newlatlonalt.lat * rad2deg);
  printdiff("Longitude", oldlatlonalt.lon * rad2deg, newlatlonalt.lon * rad2deg);
  printdiff("Altitude", oldlatlonalt.alt, newlatlonalt.alt);
  printdiff("ECEF observer i", obsposecef.i, obsposecef.i);
  printdiff("ECEF observer j", obsposecef.j, obsposecef.j);
  printdiff("ECEF observer k", obsposecef.k, obsposecef.k);
  printdiff("Range (my)", ecef2range(&obsposecef, &oldposecef), ecef2range(&obsposecef, &newposecef));
  printdiff("Range (Vallado)", orho, nrho);
  printdiff("Right ascension", ortasc, nrtasc);
  printdiff("Declination", odecl, ndecl);
  printdiff("R. rate (Vallado)", odrho, ndrho);
  printdiff("R.A. rate", odrtasc, ndrtasc);
  printdiff("Decl. rate", oddecl, nddecl);
  printf("+-------------------+---------------+---------------+------------------+\n");

/*
  // Calculate doppler shift
  double c = 299792.458;
  double f0down = 145800000;
  double f0up = 145200000;
  double downshift = (-drho / c) * f0down;
  double upshift = (drho / c) * f0up;

  printf("Doppler shifts:\nDownlink:\t%f\nUplink:\t\t%f\n\n", downshift, upshift);
*/

  return 0;
}

