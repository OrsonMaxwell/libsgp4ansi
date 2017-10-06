/*
 * vectors.c - Vector math routines for libsgp4ansi.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#include <math.h>

#include "libsgp4ansi.h"
#include "vector.h"

/*
 * Find magnitude of a 3D vector
 *
 * Inputs:  arg    - Source vector
 * Outputs: None
 * Returns: vector arg magnitude
 */
double
vec3_mag
(
  const vec3* arg
)
{
  return sqrt(arg->x * arg->x + arg->y * arg->y + arg->z * arg->z);
}

/*
 * Copy a 3D vector
 *
 * Inputs:  from - Source vector
 * Outputs: to   - Destination vector
 * Returns: None
 */
void
vec3_copy
(
  const vec3* from,
        vec3* to
)
{
  to->x = from->x;
  to->y = from->y;
  to->z = from->z;
}

/*
 * Add two 3D vectors with coefficients
 *
 * Inputs:  c1     - 1st vector coefficient
 *          vect1  - 1st vector
 *          c2     - 2nd vector coefficient
 *          vect2  - 2nd vector
 * Outputs: result - Resulting vector
 * Returns: None
 */
void
vec3_add
(
        double c1,
  const vec3*  vect1,
        double c2,
  const vec3*  vect2,
        vec3*  result
)
{
  result->x = c1 * vect1->x + c2 * vect2->x;
  result->y = c1 * vect1->y + c2 * vect2->y;
  result->z = c1 * vect1->z + c2 * vect2->z;
}

/*
 * Dot product of two 3D vectors
 *
 * Inputs:  vect1  - 1st vector
 *          vect2  - 2nd vector
 * Outputs: None
 * Returns: dot product
 */
double
vec3_dot
(
  const vec3* vect1,
  const vec3* vect2
)
{
  return vect1->x * vect2->x + vect1->y * vect2->y + vect1->z * vect2->z;
}

/*
 * Crossing of two 3D vectors
 *
 * Inputs:  vect1  - 1st vector
 *          vect2  - 2nd vector
 * Outputs: result - Resulting vector
 * Returns: None
 */
void
vec3_cross
(
  const vec3* vect1,
  const vec3* vect2,
        vec3* result
)
{
  result->x= vect1->y * vect2->z - vect1->z * vect2->y;
  result->y= vect1->z * vect2->x - vect1->x * vect2->z;
  result->z= vect1->x * vect2->y - vect1->y * vect2->x;
}

/*
 * Planar angle between two 3D vectors
 *
 * Inputs:  vect1  - 1st vector
 *          vect2  - 2nd vector
 * Returns: angle, rad
 */
double
vec3_angle
(
  const vec3* vect1,
  const vec3* vect2
)
{
  double magv1 = vec3_mag(vect1);
  double magv2 = vec3_mag(vect2);

  if (magv1 * magv2 <= 1.0e-16)
  {
    return NAN;
  }

  double temp = vec3_dot(vect1, vect2) / (magv1 * magv2);

  if (fabs(temp) > 1)
  {
    temp = (temp >= 0)?(1):(-1) * 1;
  }

  return acos(temp);
}
