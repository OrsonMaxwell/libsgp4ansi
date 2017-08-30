#ifndef LIBTLE2TCP_H_
#define LIBTLE2TCP_H_

#include <time.h>
#include <inttypes.h>

// ************************************************************************* //
//                                CUSTOM TYPES                               //
// ************************************************************************* //

/*
 * NORAD TLE representation
 */
typedef struct tle {
  char         name[24];
  unsigned int number;
  char         class;
  char         designator[8];
  time_t       epoch;
  double       ndotdiv2;
  double       nddotdiv6;
  double       bstar;
  uint8_t      ephem_type;
  unsigned int element_set;
  double       inclination;
  double       right_asc;
  double       eccentricity;
  double       arg_of_perigee;
  double       mean_anomaly;
  double       mean_motion;
  unsigned int rev_number;
} tle;

/*
 * Orbital element set
 */
typedef struct orbit
{
  long int  satnum;
  int       epochyr, epochtynumrev;
  int       error;
  char      operationmode;
  char      init, method;

  /* Near Earth */
  int    isimp;
  double aycof  , con41  , cc1    , cc4      , cc5    , d2      , d3   , d4    ,
         delmo  , eta    , argpdot, omgcof   , sinmao , t       , t2cof, t3cof ,
         t4cof  , t5cof  , x1mth2 , x7thm1   , mdot   , nodedot, xlcof , xmcof ,
         nodecf;

  /* Deep Space */
  int    irez;
  double d2201  , d2211  , d3210  , d3222    , d4410  , d4422   , d5220 , d5232 ,
         d5421  , d5433  , dedt   , del1     , del2   , del3    , didt  , dmdt  ,
         dnodt  , domdt  , e3     , ee2      , peo    , pgho    , pho   , pinco ,
         plo    , se2    , se3    , sgh2     , sgh3   , sgh4    , sh2   , sh3   ,
         si2    , si3    , sl2    , sl3      , sl4    , gsto    , xfact , xgh2  ,
         xgh3   , xgh4   , xh2    , xh3      , xi2    , xi3     , xl2   , xl3   ,
         xl4    , xlamo  , zmol   , zmos     , atime  , xli     , xni;

  double a      , altp   , alta   , epochdays, jdsatepoch       , nddot , ndot  ,
         bstar  , rcse   , inclo  , nodeo    , ecco             , argpo , mo    ,
         no;
} orbit;

// ************************************************************************* //
//                                 FUNCTIONS                                 //
// ************************************************************************* //

extern int
sgp4_init(tle* t, orbit o);

extern int
azel(void);

#endif /* LIBTLE2TCP_H_ */
