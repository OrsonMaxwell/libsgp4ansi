#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef __unix__
#include <unistd.h>
#elif defined(_WIN32) || defined(WIN32)
#include <windows.h>
#endif

#include "libsgp4ansi.h"

// Get fancy name for skylight type
void
skylight_name(char* string, skylight sky)
{
  switch (sky)
  {
    case Nighttime:
      strcpy(string, "Night time");
      break;
    case Astronomical:
      strcpy(string, "Astronomical twilight");
      break;
    case Nautical:
      strcpy(string, "Nautical twilight");
      break;
    case Civil:
      strcpy(string, "Civil twilight");
      break;
    default:
      strcpy(string, "Daytime");
      break;
  }
}

// Get fancy name for moon phases
void
moon_phase(char* string, double phase)
{
  if (phase > 0)
  {
    if (fabs(phase) >= 0.96)
    {
      strcpy(string, "Full moon");
    }
    else if (fabs(phase) >= 0.7)
    {
      strcpy(string, "Waxing gibbous");
    }
    else if (fabs(phase) >= 0.3)
    {
      strcpy(string, "First quarter");
    }
    else if (fabs(phase) >= 0.04)
    {
      strcpy(string, "Waxing crescent");
    }
    else
    {
      strcpy(string, "New moon");
    }
  }
  else
  {
    if (fabs(phase) >= 0.96)
    {
      strcpy(string, "Full moon");
    }
    else if (fabs(phase) >= 0.7)
    {
      strcpy(string, "Waning gibbous");
    }
    else if (fabs(phase) >= 0.3)
    {
      strcpy(string, "Last quarter");
    }
    else if (fabs(phase) >= 0.04)
    {
      strcpy(string, "Waning crescent");
    }
    else
    {
      strcpy(string, "New moon");
    }
  }
}

int
main (int argc, char** argv)
{
  if (argc < 2)
    return 0;
  if ((argv[1][0] != 'c')
      && (argv[1][0] != 'v')
      && (argv[1][0] != 't')
      && (argv[1][0] != 'o')
      && (argv[1][0] != 'p'))
  {
    return 1;
  }

  char tlestr0[130];
  char tlestr1[130];
  char tlestr2[130];

  char timedef[36];

  FILE* tle_file, * outfile;

  vec3 posteme = {0};
  vec3 velteme = {0};

  double t_start = -1440, t_stop = 1440, deltamin = 20;

  sat s = {0};

  if (argv[1][0] == 'o')
  {
    sat_load_tle("ISS (ZARYA)",
                 "1 25544U 98067A   17304.43179856  .00002963  00000-0  51675-4 0  9996",
                 "2 25544  51.6415  79.1213 0005738  79.0270 358.3133 15.54307707 82943",
                 &s);

//    sat_load_tle("FENGYUN 2E",
//                 "1 33463U 08066A   17281.80233449 -.00000213  00000-0  00000+0 0  9997",
//                 "2 33463   2.3909  67.6760 0005346 217.9202 107.6646  1.00271179 32268",
//                 &s);

    struct tm t = {
      .tm_year  = 117,
      .tm_mon   = 10,
      .tm_mday  = 8,
      .tm_hour  = 0,
      .tm_min   = 0,
      .tm_sec   = 0,
      .tm_isdst = 0
    };

    double time_ms = 0;

    time_t timestamp = mktime(&t) - TIMEZONE;
    vec3   obs_geo   = {54.9246 * DEG2RAD, 38.0475 * DEG2RAD, 0.180};
    obs    o = {0};

    while (true)
    {
    timestamp = time(0);

    sat_observe(&s, timestamp, time_ms, &obs_geo, &o);
    printf("Lat:    %13.3lf deg\n",  o.latlonalt.lat * RAD2DEG);
    printf("Lon:    %13.3lf deg\n",  o.latlonalt.lon * RAD2DEG);
    printf("Alt:    %13.3lf km\n",   o.latlonalt.alt);
    printf("Vel:    %13.3lf km/s\n", o.velocity);
    printf("Az:     %13.3lf deg\n",  o.azelrng.az * RAD2DEG);
    printf("El      %13.3lf deg\n",  o.azelrng.el * RAD2DEG);
    printf("Range:  %13.3lf km\n",   o.azelrng.rng);
    printf("RRate:  %13.3lf km/s\n", o.rng_rate);
    printf("Illum:  %9d\n",          o.is_illum);
    printf("-------------- The Sun --------------\n");
    printf("Az:     %13.3lf deg\n", o.sun_azelrng.az * RAD2DEG);
    printf("El:     %13.3lf deg\n", o.sun_azelrng.el * RAD2DEG);
    printf("Range:  %13.3lf km\n", o.sun_azelrng.rng);
    printf("-------------- The Moon -------------\n");
    printf("Az:     %13.3lf deg\n", o.moon_azelrng.az * RAD2DEG);
    printf("El:     %13.3lf deg\n", o.moon_azelrng.el * RAD2DEG);
    printf("Range:  %13.3lf km\n", o.moon_azelrng.rng);
    printf("Illum:  %13.0lf %%\n", fabs(o.moon_phase * 100));
    printf("Phase:  %13s\n", (o.moon_phase > 0)?("Waxing"):("Waning"));
    printf("-------------------------------------\n");

#ifdef __unix__
    usleep(1000000);
#elif defined(_WIN32) || defined(WIN32)
    Sleep(1000);
#endif

    }
    return 0;
  }

  if (argv[1][0] == 'p')
    {
      sat_load_tle("ISS (ZARYA)",
                 "1 25544U 98067A   17304.43179856  .00002963  00000-0  51675-4 0  9996",
                 "2 25544  51.6415  79.1213 0005738  79.0270 358.3133 15.54307707 82943",
                 &s);

//         sat_load_tle("JUGNU",
//                   "1 37839U 11058B   17281.88493397  .00000316  00000-0  27005-4 0  9993",
//                   "2 37839  19.9607 121.0705 0018917 356.6865 128.9420 14.12595841309799",
//                   &s);


//      sat_load_tle("???",
//                   "1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813",
//                   "2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656",
//                   &s);

//      sat_load_tle("FENGYUN 2E",
//                 "1 33463U 08066A   17281.80233449 -.00000213  00000-0  00000+0 0  9997",
//                 "2 33463   2.3909  67.6760 0005346 217.9202 107.6646  1.00271179 32268",
//                 &s);

      struct tm t = {
        .tm_year  = 117,
        .tm_mon   = 9,
        .tm_mday  = 10,
        .tm_hour  = 17,
        .tm_min   = 0,
        .tm_sec   = 0,
        .tm_isdst = 0
      };

//      struct tm t = {
//        .tm_year  = 117,
//        .tm_mon   = 10,
//        .tm_mday  = 1,
//        .tm_hour  = 0,
//        .tm_min   = 0,
//        .tm_sec   = 0,
//        .tm_isdst = 0
//      };

      vec3   observer_geo  = {54.9246 * DEG2RAD, 38.0475 * DEG2RAD, 0.180};
      time_t start_time    = mktime(&t) - TIMEZONE;
      time_t stop_time     = mktime(&t) - TIMEZONE + 7 * 1440 * 60;
      //time_t stop_time     = mktime(&t) - TIMEZONE + 3600;
      pass*  passes;
      char   buff[70];

      unsigned int delta_t = 60;
      double       horizon = 2 * DEG2RAD;

      // Important heuristic!
      unsigned int  max_passes  = (unsigned int)ceil((stop_time - start_time)
                                   / delta_t) / s.period * 2 + 1;

      passes = calloc(max_passes, sizeof(pass));

      sat_find_passes(&s, &start_time, &stop_time, &observer_geo, delta_t,
                      horizon, passes);

      printf("+--------------------------------------------------------------------+\n");
      printf("|                      Running pass prediction                       |\n");
      printf("+-----------+--------------------------------------------------------+\n");
      printf("| Satellite | %-55s|\n", s.name);
      strftime(buff, sizeof buff, "%Y-%m-%d %H:%M:%S", gmtime(&start_time));
      printf("| Start     | %-55s|\n", buff);
      strftime(buff, sizeof buff, "%Y-%m-%d %H:%M:%S", gmtime(&stop_time));
      printf("| Stop      | %-55s|\n", buff);
      printf("+-----------+-----------------------+---------+---------+------------+\n");
      printf("|   Event   |         Time          | Az, deg | El, deg |  Range, km |\n");
      printf("+-----------+-----------------------+---------+---------+------------+\n");

      unsigned int i = 0;
      while ((i < max_passes) && (passes[i].tca.el > horizon))
      {
        strftime(buff, sizeof buff, "%Y-%m-%d %H:%M:%S", gmtime(&passes[i].aos_t));
        printf("| AOS       | %-21s | %7.2lf | %7.2lf | %10.3lf |\n", buff, passes[i].aos.az * RAD2DEG, passes[i].aos.el * RAD2DEG, passes[i].aos.rng);
        strftime(buff, sizeof buff, "%Y-%m-%d %H:%M:%S", gmtime(&passes[i].tca_t));
        printf("| TCA       | %-21s | %7.2lf | %7.2lf | %10.3lf |\n", buff, passes[i].tca.az * RAD2DEG, passes[i].tca.el * RAD2DEG, passes[i].tca.rng);
        strftime(buff, sizeof buff, "%Y-%m-%d %H:%M:%S", gmtime(&passes[i].los_t));
        printf("| LOS       | %-21s | %7.2lf | %7.2lf | %10.3lf |\n", buff, passes[i].los.az * RAD2DEG, passes[i].los.el * RAD2DEG, passes[i].los.rng);
        strftime(buff, sizeof buff, "%Y-%m-%d %H:%M:%S", gmtime(&passes[i].flare_t));
        printf("| Flare     | %-21s | %7.2lf | %7.2lf | %10.3lf |\n", buff, passes[i].flare.az * RAD2DEG, passes[i].flare.el * RAD2DEG, passes[i].flare.rng);
        strftime(buff, sizeof buff, "%Y-%m-%d %H:%M:%S", gmtime(&passes[i].eclipse_t));
        printf("| Eclipse   | %-21s | %7.2lf | %7.2lf | %10.3lf |\n", buff, passes[i].eclipse.az * RAD2DEG, passes[i].eclipse.el * RAD2DEG, passes[i].eclipse.rng);
        skylight_name(buff, passes[i].sky);
        printf("| Skylight  | %-21s |         |         |            |\n", buff);
        moon_phase(buff, passes[i].moon_phase);
        printf("| Moon disc | %3.0lf%% %16s |         |         |            |\n", fabs(passes[i].moon_phase) * 100, buff);
        printf("+-----------+-----------------------+---------+---------+------------+\n");

        i++;
      }
      free(passes);

      return 0;
    }

  if (argv[1][0] == 't')
  {
    // Timing run
    printf("libsgp4ansi v%s: timing run\n", libsgp4ansi_version);

    strcpy(tlestr0, "DELTA 1");
    strcpy(tlestr1, "1 06251U 62025E   06176.82412014  .00008885  00000-0  12808-3 0  3985");
    strcpy(tlestr2, "2 06251  58.0579  54.0425 0030035 139.1568 221.1854 15.56387291  6774");

    for (int t = 0; t < 100000; t++) {
      sat_load_tle(tlestr0, tlestr1, tlestr2, &s);
    }

    for (int t = 0; t < 100000; t++) {
      t_start = t % 2880 - 1440;
      sat_propagate(&s, t_start, 10, 1.0e-12, &posteme, &velteme);
    }

    strcpy(tlestr0, "SL - 6 R / B(2)");
    strcpy(tlestr1, "1 16925U 86065D   06151.67415771  .02550794 -30915-6  18784-3 0  4486");
    strcpy(tlestr2, "2 16925  62.0906 295.0239 5596327 245.1593  47.9690  4.88511875148616");

    for (int t = 0; t < 100000; t++) {
      sat_load_tle(tlestr0, tlestr1, tlestr2, &s);
    }

    for (int t = 0; t < 100000; t++) {
      t_start = t % 2880 - 1440;
      sat_propagate(&s, t_start, 10, 1.0e-12, &posteme, &velteme);
    }

    strcpy(tlestr0, "MOLNIYA 2-14");
    strcpy(tlestr1, "1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813");
    strcpy(tlestr2, "2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656");

    for (int t = 0; t < 100000; t++) {
      sat_load_tle(tlestr0, tlestr1, tlestr2, &s);
    }

    for (int t = 0; t < 100000; t++) {
      t_start = t % 2880 - 1440;
      sat_propagate(&s, t_start, 10, 1.0e-12, &posteme, &velteme);
    }

    strcpy(tlestr0, "ITALSAT 2");
    strcpy(tlestr1, "1 24208U 96044A   06177.04061740 -.00000094  00000-0  10000-3 0  1600");
    strcpy(tlestr2, "2 24208   3.8536  80.0121 0026640 311.0977  48.3000  1.00778054 36119");

    for (int t = 0; t < 100000; t++) {
      sat_load_tle(tlestr0, tlestr1, tlestr2, &s);
    }

    for (int t = 0; t < 100000; t++) {
      t_start = t % 2880 - 1440;
      sat_propagate(&s, t_start, 10, 1.0e-12, &posteme, &velteme);
    }

    return 0;
  }

  tle_file = fopen(argv[2], "r");

  if (!tle_file) {
    printf("File opening failed!\n");
    return 0;
  }

  if (argv[1][0] == 'v')
  {
    printf("libsgp4ansi v%s: verification run\n", libsgp4ansi_version);
    outfile  = fopen("ansi_ver.out", "w");
  }
  else
  {
    printf("libsgp4ansi v%s: full catalogue run\n", libsgp4ansi_version);
    outfile  = fopen("ansi.out", "w");
  }

  double p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper;

  coe e;

  while (feof(tle_file) == 0)
  {
    fgets(tlestr0, 130, tle_file);
    fgets(tlestr1, 130, tle_file);
    fgets(tlestr2, 130, tle_file);

    if (argv[1][0] == 'v'){
      strncpy(timedef, tlestr2 + 69, 35);
      timedef[35] = '\0';
      sscanf(timedef, "%lf %lf %lf", &t_start, &t_stop, &deltamin);
      strncpy(tlestr2, tlestr2, 69);
    }

    sat_load_tle(tlestr0, tlestr1, tlestr2, &s);

    fprintf(outfile, "%u (%.8e)\n", s.norad_number, s.period);

    // Iterate over time range
    for (double t = t_start; t <= t_stop + deltamin - 1.0e-12; t += deltamin)
    {
      int retval = sat_propagate(&s, t, 10, 1.0e-12, &posteme, &velteme);

      if (retval != 0)
      {
        printf("[ERROR] Sat %5d (%.8e) code %2d at %8.f mfe\n",
               s.norad_number, s.period, retval, t);
        break;
      }
      else
      {
        fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f",
                t, posteme.x, posteme.y, posteme.z);

        if (argv[1][0] == 'v')
        {
          e = teme2coe(&posteme, &velteme);
          fprintf(outfile, " %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                  e.a, e.ecc, e.incl*RAD2DEG, e.omega*RAD2DEG, e.argp*RAD2DEG,
                  e.nu*RAD2DEG, e.m*RAD2DEG);
        }
        else
        {
          fprintf(outfile, "\n");
        }
      }
    }
  }

  fclose(tle_file);
  fclose(outfile);

  return 0;
}
