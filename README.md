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
    Thanks to this approach, this library has the facilities to not only propa-
    gate NORAD TLE element sets in time, but also to predict satellite passes
    over given observer location on the geoid as well as to find satellite tran-
    sits over solar and lunar discs;
  - Shared library format allows flexibility and compartmentalization of effort;
  - Optimization allowed to achieve temporal performance comparable to statical-
    ly linked reference AIAA paper C++ code when propagating the orbits and
    significantly better performance when finding satellite passes when compa-
    red to popular amateur software (like Orbitron);

Version history:
  0.9 Initial release. requires extensive observational fidelity testing.