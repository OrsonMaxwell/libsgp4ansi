/*
 * const.h - precomputed mathematical and physical constants for libsgp4ansi.
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#ifndef CONST_H_
#define CONST_H_

// ************************************************************************* //
//                                    MATH                                   //
// ************************************************************************* //

#define pi          3.14159265358979323846
#define twopi       6.28318530717958647692
#define threepidiv2 4.71238898038468985769
#define pidiv2      1.57079632679489661923

#define deg2rad     0.01745329251994329576 // Degrees to radians coefficient
#define rad2deg     57.2957795130823208767 // Radians to degrees coefficient
#define rpd2radmin  229.183118052329283507 // Revolutions per day
                                           // to radians per minute coefficient

// ************************************************************************* //
//                                   PHYSICS                                 //
// ************************************************************************* //

#define G         6.67408e-11            // Newtonian gravitational constant
#define M         5.9721986e24           // Mass of the Earth, kg
#define ke        1.996470165874962430e7 // ??????

// ************************************************************************* //
//                               GEODETIC DATUM                              //
// ************************************************************************* //

#ifdef USE_WGS72  // WGS 72
#define mu        398600.8               // km3/s2
#define Re        6378.135               // Earth radius, km
#define xke       0.07436691613317341324
#define tumin     13.4468396969593099729
#define j2        0.001082616            // 2nd grav zonal harmonic of the Earth
#define j3       -0.00000253881          // 3rd grav zonal harmonic of the Earth
#define j4       -0.00000165597          // 4th grav zonal harmonic of the Earth
#define j3divj2  -0.00234506972001152763
#endif /* USE_WGS72 */

#ifdef USE_WGS84  // WGS 84
#define mu        398600.5               // km3 / s2
#define Re        6378.137               // Earth radius, km
#define xke       0.07436685316871384602
#define tumin     13.4468510820449809410
#define j2        0.00108262998905       // 2nd grav zonal harmonic of the Earth
#define j3       -0.00000253215306       // 3rd grav zonal harmonic of the Earth
#define j4       -0.00000161098761       // 4th grav zonal harmonic of the Earth
#define j3divj2  -0.00233889055874200014
#endif /* USE_WGS84 */

#endif /* CONST_H_ */
