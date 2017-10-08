/*
 * transform.h - Time format routines for libsgp4ansi.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#include <time.h>
#include <math.h>

#include "libsgp4ansi.h"
#include "const.h"
#include "epoch.h"

/*
 * Convert year and fractional day to unix time
 *
 * Inputs:  year   - Year
 *          days   - Decimal day since year start
 * Outputs: unix   - Unix time, s since 00:00 Jan, 1 1970
 *          ms     - Fractional part of seconds, ms
 * Returns: 0      - Success
 *         -1      - Invalid input
 */
int
fractday2unix
(
  unsigned int year,
  double       days,
  time_t*      unix,
  float*       ms
)
{
  if (days > ((year % 4 == 0) ? 366 : 365))
  {
    return -1;
  }

  struct tm res_tm;

  int mon_len[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  int day_of_year = (int)floor(days);

  res_tm.tm_year = 100 + year;

  if (year >= 57) // Will be valid only until 2057!
    res_tm.tm_year -= 100;

  // Month and day of month
  if ((year % 4) == 0) // Leap year?
    mon_len[1] = 29;

  int i = 1, j = 0;
  while ((day_of_year > j + mon_len[i - 1]) && (i < 12))
  {
    j = j + mon_len[i-1];
    i++;
  }
  res_tm.tm_mon = i - 1;
  res_tm.tm_mday = day_of_year - j;

  // Hours, minutes, and seconds
  double temp, sec;
  int result;
  temp           = (days - day_of_year) * 24;
  res_tm.tm_hour = (int)floor(temp);
  temp           = (temp - res_tm.tm_hour) * 60;
  res_tm.tm_min  = (int)floor(temp);
  sec            = (temp - res_tm.tm_min) * 60;
  res_tm.tm_sec  = (int)floor(sec);
  *ms            = (sec - res_tm.tm_sec) * 1000;

  // Ignore DST
  res_tm.tm_isdst = 0;
  result = mktime(&res_tm) - TIMEZONE;

  if (result == -1)
  {
    return result;
  }

  *unix = (time_t)result;
  return 0;
}

/*
 * Convert unix time to Julian date
 *
 * Inputs:  time - Timestamp in unix time
 *          ms   - Fracitonal second part, ms
 * Returns: Julian date on success
 */
double
unix2jul
(
  const time_t* time,
        float   ms
)
{
  struct tm* t;
  t = gmtime(time);

  return 367 * (t->tm_year + 1900)
  - floor((7 * ((t->tm_year + 1900) + floor((t->tm_mon + 10) / 12))) * 0.25)
  + floor(275 * (t->tm_mon + 1) / 9)
  + t->tm_mday + 1721013.5
  + ((((double)t->tm_sec + ms / 1000) / 60 + t->tm_min) / 60 + t->tm_hour) / 24;
}

/*
 * Convert Julian date to Greenwich Siderial Time
 *
 * Inputs:  julian - Julian date
 * Returns: GST time, rad
 */
double
jul2gst
(
  double julian
)
{
  double result, tempUT1;

  tempUT1 = (julian - JAN1_2000_1200H) / 36525.0;
  result = -6.2e-6* tempUT1 * tempUT1 * tempUT1 + 0.093104 * tempUT1 * tempUT1 +
      (876600.0*3600 + 8640184.812866) * tempUT1 + 67310.54841;

  result = fmod(result * DEG2RAD / 240.0, TAU);

  // Check quadrants
  if (result < 0.0)
    result += TAU;

  return result;
}
