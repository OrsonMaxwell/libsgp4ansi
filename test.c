#include <stdio.h>
#include <string.h>
#include <time.h>

#include "libsgp4ansi.h"

int
main (int argc, char** argv)
{
  if (argc != 3)
    return 0;
  if ((argv[1][0] != 'c') && (argv[1][0] != 'v') && (argv[1][0] != 'h'))
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

  int error = 0;

  tle_file = fopen(argv[2], "r");

  if (!tle_file) {
    printf("File opening failed!\n");
    return 0;
  }

  outfile  = fopen("ansi.out", "w");

  while (feof(tle_file) == 0)
  {
    fgets(tlestr0, 130, tle_file);
    fgets(tlestr1, 130, tle_file);
    fgets(tlestr2, 130, tle_file);

    if (argv[1][0] == 'v'){
      strncpy(timedef, tlestr2 + 69, 35);
      timedef[35] = '\0';
      sscanf(timedef, "%lf %lf %lf", &t_start, &t_stop, &deltamin);
    }

    tle2orbit(tlestr0, tlestr1, tlestr2, &s);

    fprintf(outfile, "%ld (%12.9lf)\n", s.norad_number, 3.14159265358979323846 * 2 / s.mean_motion);

    // Iterate over time range
    for (double t = t_start; t <= t_stop; t += deltamin)
    {
      error = orbit_prop(&s, t, 10, 1.0e-12, &posteme, &velteme);

      if (error != 0)
        printf("[ERROR] at %f: %3d\n", t, error);

      if (error == 0)
      {
        fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f",
                t, posteme.x, posteme.y, posteme.z);

        if (argv[1][0] == 'v')
        {
          /*
          rv2coe(ro, vo, mu, p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper );
                          fprintf(outfile, " %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                                   a, ecc, incl*rad, node*rad, argp*rad, nu*rad, m*rad);
          */
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
