#include <stdio.h>
#include <string.h>
#include <time.h>

#include "libsgp4ansi.h"

int
main (int argc, char** argv)
{


  char tlestr1[130];
  char tlestr2[130];

  FILE* tle_file, * outfile;

  vec3 posteme = {0};
  vec3 velteme = {0};

  double t_start = -1440, t_stop = 1440, deltamin = 20;

  sat s = {0};

  char str[2];

  int error = 0;

  tle_file = fopen(argv[1], "r");
  //tle_file = fopen("full.tle", "r");

  if (!tle_file) {
    printf("File opening failed!\n");
    return 0;
  }

  outfile  = fopen("ansi.out", "w");

  while (feof(tle_file) == 0)
  {
    do
    {
      fgets(tlestr1, 130, tle_file);
      strncpy(str, &tlestr1[0], 1);
      str[1] = '\0';
    // TODO: Get read of #
    } while ((strcmp(str, "#") == 0) && (feof(tle_file) == 0));

    if (feof(tle_file) == 0)
    {
      fgets(tlestr2, 130, tle_file);

      // TODO: Read name from file!
      tle2orbit(tlestr1, tlestr2, &s);

      fprintf(outfile, "%ld (%12.9lf)\n", s.norad_number, 3.14159265358979323846 * 2 / s.mean_motion);

      // Iterate over time range
      for (double t = t_start; t <= t_stop; t += deltamin)
      {
        error = orbit_prop(&s, t, 10, 1.0e-12, &posteme, &velteme);

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
