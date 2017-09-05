#include <stdio.h>
#include <string.h>
#include <time.h>

#include "libsgp4ansi.h"
#include "const.h"
#include "transform.h"

int
main (int argc, char** argv)
{

  char tlestr1[130];
  char tlestr2[130];

  FILE* tle_file, * outfile;

  vect posteme;
  vect velteme;

  double startmfe = -1440, stopmfe = 1440, deltamin = 20, tsince = 0;

  orbit sat = {0};

  time_t prop_time;

  char str[2];

  int error = 0;

  tle_file = fopen("full.tle", "r");
  outfile  = fopen("ansi.out", "w");

  if(!tle_file) {
    perror("File opening failed!");
    return 1;
  }

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

      fprintf(outfile, "%ld (%12.9lf)\n", sat.number, pi * 2 / sat.no);

      // Iterate over time range
      tsince = startmfe;
      while (tsince < stopmfe)
      {
        tsince += deltamin;

        if(tsince > stopmfe)
          tsince = stopmfe;

        prop_time += tsince * 60;

        orbit_prop(&sat, &prop_time, sat.epoch_ms, 10, 1.0e-12, &posteme, &velteme);

        if (error > 0)
          printf("[ERROR] at %f: %3d\n", tsince, error);

        if (error == 0)
        {
          fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f\n",
                  tsince, posteme.x, posteme.y, posteme.z);
        }
      }
    }
  }

  fclose(tle_file);
  fclose(outfile);

  return 0;
}
