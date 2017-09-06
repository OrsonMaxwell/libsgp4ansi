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

// Find magnitude of a 3D vector
double
magvec3(vec3*);

// Copy a 3D vector
void
copyvec3(vec3*, vec3*);

// Add two 3D vectors with coefficients
void
addvec3(double, vec3*, double, vec3*, vec3*);

// Dot product of two 3D vectors
double
dotvec3(vec3*, vec3*);

// Crossing of two 3D vectors
void
crossvec3(vec3*, vec3*, vec3*);

// TODO: Write decription or inline the procedures
void rotate2(vec3*, double, vec3*);
void rotate3(vec3*, double, vec3*);

#endif /* VECTOR_H_ */
