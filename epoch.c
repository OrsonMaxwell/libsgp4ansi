/*
 * epoch.c - Time routines for libsgp4ansi.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 * Fundamentals of Astrodynamics and Applications, D. Vallado, Second Edition
 * Astronomical Algorithms, Jean Meeus
 *
 * Copyright � 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#include <time.h>
#include <math.h>

#include "libsgp4ansi.h"
#include "const.h"
#include "epoch.h"

/*
 * Convert year and fractional day to unix time
 *
 * Inputs:  year - Year
 *          days - Decimal day since year start
 * Outputs: ms   - Millisecond portion of time
 * Returns: unix - Unix time, s since 00:00 Jan, 1 1970
 */
time_t
fractday2unix
(
  unsigned int year,
  double       days,
  float*       time_ms
)
{
  unsigned int days_in_yr = (year % 4 == 0) ? 366 : 365;
  days = (days > days_in_yr) ? (days_in_yr) : (days);

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

  if (time_ms != NULL)
  {
    *time_ms          = (sec - res_tm.tm_sec) * 1000;
  }

  // Ignore DST
  res_tm.tm_isdst = 0;
  result = mktime(&res_tm) - TIMEZONE;

  if (result == -1)
  {
    return result;
  }

  return (time_t)result;
}

/*
 * Convert unix time to Julian date
 *
 * Inputs:  time - Timestamp in unix time
 *          ms   - Millisecond portion of time
 * Returns: Julian date on success
 */
double
unix2jul
(
  time_t time,
  float  time_ms
)
{
  struct tm* t;
  t = gmtime(&time);

  time_ms = fmin(time_ms, nextafter(1000, 0));

  return 367 * (t->tm_year + 1900)
  - floor((7 * ((t->tm_year + 1900) + floor((t->tm_mon + 10) / 12))) * 0.25)
  + floor(275 * (t->tm_mon + 1) / 9)
  + t->tm_mday + 1721013.5
  + ((((double)t->tm_sec + time_ms / 1000) / 60 + t->tm_min) / 60 + t->tm_hour)
  / 24;
}

/*
 * Convert Julian date to Greenwich Siderial Time
 *
 * Inputs:  julian_date - Julian date
 * Returns: GST time, rad
 */
double
jul2gst
(
  double julian_date
)
{
  double gst, UT1;

  UT1 = (julian_date - J2000) / 36525.0;
  gst = -6.2e-6 * UT1 * UT1 * UT1 + 0.093104 * UT1 * UT1 +
      (876600.0 * 3600 + 8640184.812866) * UT1 + 67310.54841;

  gst = fmod(gst * DEG2RAD / 240.0, TAU);

  // Check quadrants
  if (gst < 0.0)
    gst += TAU;

  return gst;
}
