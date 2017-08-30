#ifndef TLE_H_
#define TLE_H_

#include <time.h>
#include <inttypes.h>

typedef struct tle {
  char         name[24];
  unsigned int number;
  char         class;
  char         designator[8];
  time_t       epoch;
  double       ndotdiv2;
  double       nddotdiv6;
  double       bstar;
  uint8_t      ephem_type;
  unsigned int element_set;
  double       inclination;
  double       right_asc;
  double       eccentricity;
  double       arg_of_perigee;
  double       mean_anomaly;
  double       mean_motion;
  unsigned int rev_number;
} tle;

#endif /* TLE_H_ */
