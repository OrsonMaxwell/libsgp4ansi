/*
 * vectors.h - Vector math routines for libsgp4ansi.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#ifndef VECTOR_H_
#define VECTOR_H_

// ************************************************************************* //
//                                 FUNCTIONS                                 //
// ************************************************************************* //

// Find magnitude of a 3D vector
double
vec3_mag(vec3*);

// Copy a 3D vector
void
vec3_copy(vec3*, vec3*);

// Add two 3D vectors with coefficients
void
vec3_add(double, vec3*, double, vec3*, vec3*);

// Dot product of two 3D vectors
double
vec3_dot(vec3*, vec3*);

// Crossing of two 3D vectors
void
vec3_cross(vec3*, vec3*, vec3*);

#endif /* VECTOR_H_ */
