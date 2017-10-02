#include <stdio.h>
#include <string.h>
#include <time.h>

#include "libsgp4ansi.h"
#include "const.h"

int
main (int argc, char** argv)
{
  if (argc < 2)
    return 0;
  if ((argv[1][0] != 'c') && (argv[1][0] != 'v') && (argv[1][0] != 't'))
    return 0;

  char tlestr0[130];
  char tlestr1[130];
  char tlestr2[130];

  char timedef[36];

  FILE* tle_file, * outfile;

  vec3 posteme = {0};
  vec3 velteme = {0};

  double t_start = -1440, t_stop = 1440, deltamin = 20;

  sat s = {0};

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

    fprintf(outfile, "%ld (%.8e)\n", s.norad_number, s.period);

    // Iterate over time range
    for (double t = t_start; t <= t_stop + deltamin - 1.0e-12; t += deltamin)
    {
      int retval = sat_propagate(&s, t, 10, 1.0e-12, &posteme, &velteme);
      //int retval = 0;

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
          double ro[3] = {posteme.x, posteme.y, posteme.z};
          double vo[3] = {velteme.x, velteme.y, velteme.z};

          rv2coe(ro, vo, &p, &a, &ecc, &incl, &node, &argp, &nu, &m,
                 &arglat, &truelon, &lonper);
          fprintf(outfile, " %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                  a, ecc, incl*RAD2DEG, node*RAD2DEG, argp*RAD2DEG, nu*RAD2DEG, m*RAD2DEG);
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
