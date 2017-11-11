Description:

  LibSGP4ANSI is an ANSI C11 shared library implementation of the SGP4 and SDP4
  algorithms with all the corrections according to the AIAA 2006-6753 paper
  "Revisiting Spacetrack Report #3" by Vallado, Crawford, Hujsak, and Kelso.

Motivation:

  Apparently the author was unable to locate a singular convenient project which
  could be used as a mathematical core for satellite propagation software of his
  own. Some of the implementations out there seemed inefficient and sometimes an
  overkill in their approach to scientific computing.
  This particular implementation was concieved with practical applications for
  radio amateurs and astro photography enthusiasts in mind as well as carrying
  an educational value for its author. Having said all this, this code base
  possesses the following notable features:
  
  - Cleaner code. Most of the historic routines that were migrated directly from
    FORTRAN sources by other authors were inlined into the main initialization
    and propagation routines saving up on call overheads using 40+ arguments.
    Also a light refactoring of the naming was in order as well as significantly
    reducing the number of variables used.
    Common constants were moved into macros, all relevant code subsets were 
    organized into respective translation units;
  - Condensed feature list. The author chose to concentrate on what was required
    for practical applications instead of trying to implement an exhaustive set
    of downstream astrodynamical mathematics;
    Thanks to this approach, this library has the facilities to not only 
    propagate NORAD TLE element sets in time, but also to predict satellite
    passes over given observer location on the geoid as well as to find
    satellite transits over solar and lunar discs;
  - Shared library format allows flexibility and compartmentalization of effort;
  - Optimization allowed to achieve temporal performance comparable to
    statically linked reference AIAA paper C++ code when propagating the orbits
    and significantly better performance when finding satellite passes when
    compared to popular amateur software (like Orbitron);

Usage:

  Basic usage consists of instantiating a 'sat' data structure and then passing
  it to either 'sat_load_tle()' function along with the TLE strings or to the 
  'sat_load_params()' function. In the latter case you will require to supply
  the mean elements encoded in TLE as arguments to that function.
  After that the 'sat' structure will contain all the mathematical terms
  required to propagate the orbit in time expanded according to SGP4/SDP4
  specification.
  
  To do so, you can use 'sat_propagate()' to obtain a TEME position and velocity
  vectors of the satellite at given time.
  
  It is also possible to use a number of downstream wrapper functions:
  - 'sat_observe()' - get observational data from a given geodetic point at a 
    given time including satellite azimuth, elevation, range, range rate,
    positions of the Sun and the Moon, skylight type etc;
  - 'sat_find_passes()' get information on all the passes the satellite makes in
    a given time frame over a given location on the geoid. Output includes times
    and observational vectors of all AOS, TCA and LOS events as well as
    information about satellite illumination by the sun etc. Be aware that this
    function requires a pointer to a preallocated array of 'pass' objects -
    plese refer to test.c for a useful heuristic that helps determine how many
    'pass' elements to allocate memory for;
  - 'sat_find_transits()' - get information about the transits which satellite
    makes over solar and lunar discs as viewed from a given geodetic location in
    a given time frame;
  - 'star_observe()' - get observational vector for a quazi-fixed celestial
    body (the one whose movement across the sky is insignificant for a single
    observation session) - stars, nebulae, and other catalogue objects as well
    as some planets. The scope of this library does not include the
    determination of true equatorial coordinates of such objects, therefore
    these coordinates must be supplied to this function explicitly.
    
   Please see the source code for more information - the author likes to think
   that it is well documented and does not require any further explanation.

Version history:

  - 0.9 Initial release. requires extensive observational fidelity testing.