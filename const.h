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

#define PI          3.141592653589793L
#define TWOPI       (PI * 2)
#define THREEPIDIV2 (3 * PI / 2)
#define PIDIV2      (PI / 2)

#define DEG2RAD     (PI / 180)             // Degrees to radians coefficient
#define RAD2DEG     (180 / PI)             // Radians to degrees coefficient
#define RPD2RADPM   (1440.0 / (2.0 * PI))  // Revolutions per day
                                           // to radians per minute coefficient

// ************************************************************************* //
//                               GEODETIC DATUM                              //
// ************************************************************************* //

#ifdef USE_WGS72_LEGACY  // Old WGS 72
#define GM        398600.79964           // mu, km3 / s2
#define RE        6378.135               // Equatorial Earth radius, km
#define ECC       0.0818188              // Eccentricity
#define OMEGAE    7.292115147e-5         // Angular velocity, rad/sec
#define RPTIM     (OMEGAE * 60)          // Angular velocity, rad/min
#define FLATT     0.00335277945416750486 // Flattening of the Earth ellipsoid
#define XKE       0.0743669161           // (60.0 / sqrt(Re*Re*Re/GM))
#define TUMIN     (1 / XKE)              // Time units in minute
#define J2        0.001082616            // 2nd grav zonal harmonic of the Earth
#define J3       -0.00000253881          // 3rd grav zonal harmonic of the Earth
#define J4       -0.00000165597          // 4th grav zonal harmonic of the Earth
#define J3DIVJ2   (J3 / J2)
#endif /* USE_WGS72_LEGACY */

#ifdef USE_WGS72  // WGS 72
#define GM        398600.8               // mu, km3 / s2
#define RE        6378.135               // Equatorial Earth radius, km
#define ECC       0.0818188              // Eccentricity
#define OMEGAE    7.292115147e-5         // Angular velocity, rad/sec
#define RPTIM     (OMEGAE * 60)          // Angular velocity, rad/min
#define FLATT     0.00335277945416750486 // Flattening of the Earth ellipsoid
#define XKE       0.07436691613317341324 // (60.0 / sqrt(Re*Re*Re/GM))
#define TUMIN     (1 / XKE)              // Time units in minute
#define J2        0.001082616            // 2nd grav zonal harmonic of the Earth
#define J3       -0.00000253881          // 3rd grav zonal harmonic of the Earth
#define J4       -0.00000165597          // 4th grav zonal harmonic of the Earth
#define J3DIVJ2   (J3 / J2)
#endif /* USE_WGS72 */

#ifdef USE_WGS84  // WGS 84
#define GM        398600.5               // mu, km3 / s2
#define RE        6378.137               // Earth radius, km
#define ECC       0.0818188              // Eccentricity
#define OMEGAE    7.292115e-5            // Angular velocity, rad/sec
#define RPTIM     (OMEGAE * 60)          // Angular velocity, rad/min
#define FLATT     0.00335281066478120473 // Flattening
#define XKE       0.07436685316871384602 // (60.0 / sqrt(Re*Re*Re/GM))
#define TUMIN     (1 / XKE)              // Time units in minute
#define J2        0.00108262998905       // 2nd grav zonal harmonic of the Earth
#define J3       -0.00000253215306       // 3rd grav zonal harmonic of the Earth
#define J4       -0.00000161098761       // 4th grav zonal harmonic of the Earth
#define J3DIVJ2   (J3 / J2)
#endif /* USE_WGS84 */

#endif /* CONST_H_ */
