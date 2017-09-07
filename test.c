#include <stdio.h>
#include <string.h>
#include <time.h>

#include "libsgp4ansi.h"
#include "const.h"
#include "vector.h"

int
main (int argc, char** argv)
{

  char tlestr1[130];
  char tlestr2[130];

  FILE* tle_file, * outfile;

  vec3 posteme;
  vec3 velteme;

  double t_start = -1440, t_stop = 1440, deltamin = 20;//, tsince = 0;

  orbit sat = {0};

  char str[2];

  int error = 0;

  tle_file = fopen(argv[1], "r");

  if(!tle_file) {
    perror("File opening failed!");
    return 1;
  }

  outfile  = fopen("ansi.out", "w");

  while (feof(tle_file) == 0)
  {
    do
    {
      fgets(tlestr1, 130, tle_file);
      strncpy(str, &tlestr1[0], 1);
      str[1] = '\0';
    } while ((strcmp(str, "#") == 0) && (feof(tle_file) == 0));

    if (feof(tle_file) == 0)
    {
      fgets(tlestr2, 130, tle_file);

      tle2orbit(tlestr1, tlestr2, &sat);

      fprintf(outfile, "%ld (%12.9lf)\n", sat.number, PI * 2 / sat.no);

      // Iterate over time range
      for (double t = t_start; t <= t_stop; t += deltamin)
      {
        error = orbit_prop(&sat, t, 10, 1.0e-12, &posteme, &velteme);
        if (error != 0)
          printf("[ERROR] at %f: %3d\n", t, error);

        if (error == 0)
        {
          fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f\n",
                  t, posteme.x, posteme.y, posteme.z);
        }
      }
    }
  }

  fclose(tle_file);
  fclose(outfile);

  return 0;
}
