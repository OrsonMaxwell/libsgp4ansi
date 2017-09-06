/*
 * vectors.c - Vector math routines for libsgp4ansi.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 *
 * Copyright � 2017 Orson J. Maxwell. Please see LICENSE for details.
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
magvec3(vec3* arg)
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
copyvec3(vec3* from, vec3* to)
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
addvec3(double c1, vec3* vect1, double c2, vec3* vect2, vec3* result)
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
dotvec3(vec3* vect1, vec3* vect2)
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
crossvec3
(
    vec3* vect1, vec3* vect2, vec3* result
)
{
  result->x= vect1->y * vect2->z - vect1->z * vect2->y;
  result->y= vect1->z * vect2->x - vect1->x * vect2->z;
  result->z= vect1->x * vect2->y - vect1->y * vect2->x;
}

void rotate2(vec3* src, double angle, vec3* result)
{
  double tempz = src->z;
  double s     = sin(angle);
  double c     = cos(angle);

  result->x = src->x * c - tempz * s;
  result->y = src->y;
  result->z = src->z * c + src->x * s;
}

void rotate3(vec3* src, double angle, vec3* result)
{
  double tempy  = src->y;
  double sine   = sin(angle);
  double cosine = cos(angle);

  result->y = src->y * cosine - src->x * sine;
  result->x = src->x * cosine + tempy * sine;
  result->z = src->z;
}
