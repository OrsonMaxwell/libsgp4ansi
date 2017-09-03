/*
 * test.c
 *
 *  Created on: 28 рту. 2017 у.
 *      Author: Orson
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

  double r[3], v[3];

  sgp4(wgs72, &satrec, 60.0, r, v);

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
  epoch_tm.tm_mon    = 7;
  epoch_tm.tm_mday   = 29;
  epoch_tm.tm_hour   = 5;
  epoch_tm.tm_min    = 1;
  epoch_tm.tm_sec    = 57;

  prop_tm.tm_year    = 117;
  prop_tm.tm_mon     = 7;
  prop_tm.tm_mday    = 29;
  prop_tm.tm_hour    = 6;
  prop_tm.tm_min     = 1;
  prop_tm.tm_sec     = 57;

  prop_time = mktime(&prop_tm) - timezone;

  strcpy(iss.name, "ISS");
  strcpy(iss.designator, "98067A  ");
  iss.number         = 25544;
  iss.sec_class      = 'U';
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
  orbit_prop(&iss, &prop_time, 0, 10, 1.0e-12, &newposteme, &newvelteme);
  // ***************************************************************************

  printf("TEME -------------------------------------------------\n");
  printf("OLD pos: %12.8f\t%12.8f\t%12.8f\nNEW pos: %12.8f\t%12.8f\t%12.8f\n",
         oldposteme.x, oldposteme.y, oldposteme.z, newposteme.x, newposteme.y, newposteme.z);
  printf("OLD vel: %12.8f\t%12.8f\t%12.8f\nNEW vel: %12.8f\t%12.8f\t%12.8f\n",
         oldvelteme.x, oldvelteme.y, oldvelteme.z, newvelteme.x, newvelteme.y, newvelteme.z);

  teme2ecef(&oldposteme, &oldposteme, unix2jul(&prop_time, 0), &oldposecef, &oldvelecef);
  teme2ecef(&newposteme, &newposteme, unix2jul(&prop_time, 0), &newposecef, &newvelecef);

  printf("ECEF -------------------------------------------------\n");
  printf("OLD pos: %12.8f\t%12.8f\t%12.8f\nNEW pos: %12.8f\t%12.8f\t%12.8f\n",
         oldposecef.i, oldposecef.j, oldposecef.k, newposecef.i, newposecef.j, newposecef.k);
  printf("OLD vel: %12.8f\t%12.8f\t%12.8f\nNEW vel: %12.8f\t%12.8f\t%12.8f\n",
         oldvelecef.i, oldvelecef.j, oldvelecef.k, newvelecef.i, newvelecef.j, newvelecef.k);

  ecef2latlonalt(&oldposecef, unix2jul(&prop_time, 0), 10, 1.0e-12, &oldlatlonalt);
  ecef2latlonalt(&newposecef, unix2jul(&prop_time, 0), 10, 1.0e-12, &newlatlonalt);

  printf("LATLONALT --------------------------------------------\n");
  printf("OLD lla: %12.8f\t%12.8f\t%12.8f\nNEW lla: %12.8f\t%12.8f\t%12.8f\n",
         oldlatlonalt.lat * rad2deg, oldlatlonalt.lon * rad2deg, oldlatlonalt.alt,
         newlatlonalt.lat * rad2deg, newlatlonalt.lon * rad2deg, newlatlonalt.alt);

  vect obsposlla = {0};  // observer lla vector at latlonalt 0,0,0
  vect dacha;
  dacha.lat = 54.924551 * deg2rad;
  dacha.lon = 38.047505 * deg2rad;
  dacha.alt = 0.180;
  vect obsposecef; // observer ECEF vector at latlonalt 0,0,0
  //obsposecef.x = 6378.137;
  //obsposecef.y = 0;
  //obsposecef.z = 0;

  latlonalt2ecef(&dacha, &obsposecef);

  printf("OBSECEF ----------------------------------------------\n");
  printf("NEW obp: %12.8f\t%12.8f\t%12.8f\n",
         obsposecef.i, obsposecef.j, obsposecef.k);
  printf("NEW rng: %12.8f\n", ecef2range(&obsposecef, &newposecef));

/*
  double orange, oaz, oel, orrate, oazrate, oelrate;
  double nrange, naz, nel, nrrate, nazrate, nelrate;

  ecef2azel(&oldposecef, &oldvelecef, &obsposecef, oldlatlonalt.lat, oldlatlonalt.lon,
            &orange, &oaz, &oel, &orrate, &oazrate, &oelrate);
  ecef2azel(&newposecef, &newvelecef, &obsposecef, newlatlonalt.lat, newlatlonalt.lon,
            &nrange, &naz, &nel, &nrrate, &nazrate, &nelrate);

  printf("AZELRANGE --------------------------------------------\n");
  printf("OLD aer: %12.8f\t%12.8f\t%12.8f\nNEW aer: %12.8f\t%12.8f\t%12.8f\n",
         oaz * rad2deg, oel * rad2deg, orange, naz * rad2deg, nel * rad2deg, nrange);
  printf("AZELRANGERATE ----------------------------------------\n");
  printf("OLD aer: %12.8f\t%12.8f\t%12.8f\nNEW aer: %12.8f\t%12.8f\t%12.8f\n",
         oazrate * rad2deg, oelrate * rad2deg, orrate, nazrate * rad2deg, nelrate * rad2deg, nrrate);
*/


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
