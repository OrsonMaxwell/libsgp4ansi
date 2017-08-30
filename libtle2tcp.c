/*
 * libtle2tcp.c
 *
 *  Created on: 28 рту. 2017 у.
 *      Author: Orson
 */

#include <stdio.h>
#include <math.h>

#include "libtle2tcp.h"
#include "sgp4unit.h"
#include "sgp4ext.h"
#include "sgp4io.h"

/*
 * Initialize the SGP4/SDP4 math with orbital elements
 */
int
sgp4_init(tle* t, orbit o)
{
  /* --------------------- local variables ------------------------ */
  double ao, ainv,   con42, cosio, sinio, cosio2, eccsq,
       omeosq, posq,   rp,     rteosq,
       cnodm , snodm , cosim , sinim , cosomm, sinomm, cc1sq ,
       cc2   , cc3   , coef  , coef1 , cosio4, day   , dndt  ,
       em    , emsq  , eeta  , etasq , gam   , argpm , nodem ,
       inclm , mm    , nm    , perige, pinvsq, psisq , qzms24,
       rtemsq, s1    , s2    , s3    , s4    , s5    , s6    ,
       s7    , sfour , ss1   , ss2   , ss3   , ss4   , ss5   ,
       ss6   , ss7   , sz1   , sz2   , sz3   , sz11  , sz12  ,
       sz13  , sz21  , sz22  , sz23  , sz31  , sz32  , sz33  ,
       tc    , temp  , temp1 , temp2 , temp3 , tsi   , xpidot,
       xhdot1, z1    , z2    , z3    , z11   , z12   , z13   ,
       z21   , z22   , z23   , z31   , z32   , z33,
       qzms2t, ss, j2, j3oj2, j4, x2o3, r[3], v[3],
       tumin, mu, radiusearthkm, xke, j3, delmotemp, qzms2ttemp, qzms24temp;

  /* ----------- set all near earth variables to zero ------------ */
  o->isimp   = 0;   o->method = 'n'; o->aycof    = 0.0;
  o->con41   = 0.0; o->cc1    = 0.0; o->cc4      = 0.0;
  o->cc5     = 0.0; o->d2     = 0.0; o->d3       = 0.0;
  o->d4      = 0.0; o->delmo  = 0.0; o->eta      = 0.0;
  o->argpdot = 0.0; o->omgcof = 0.0; o->sinmao   = 0.0;
  o->t       = 0.0; o->t2cof  = 0.0; o->t3cof    = 0.0;
  o->t4cof   = 0.0; o->t5cof  = 0.0; o->x1mth2   = 0.0;
  o->x7thm1  = 0.0; o->mdot   = 0.0; o->nodedot  = 0.0;
  o->xlcof   = 0.0; o->xmcof  = 0.0; o->nodecf   = 0.0;

  /* ----------- set all deep space variables to zero ------------ */
  o->irez  = 0;   o->d2201 = 0.0; o->d2211 = 0.0;
  o->d3210 = 0.0; o->d3222 = 0.0; o->d4410 = 0.0;
  o->d4422 = 0.0; o->d5220 = 0.0; o->d5232 = 0.0;
  o->d5421 = 0.0; o->d5433 = 0.0; o->dedt  = 0.0;
  o->del1  = 0.0; o->del2  = 0.0; o->del3  = 0.0;
  o->didt  = 0.0; o->dmdt  = 0.0; o->dnodt = 0.0;
  o->domdt = 0.0; o->e3    = 0.0; o->ee2   = 0.0;
  o->peo   = 0.0; o->pgho  = 0.0; o->pho   = 0.0;
  o->pinco = 0.0; o->plo   = 0.0; o->se2   = 0.0;
  o->se3   = 0.0; o->sgh2  = 0.0; o->sgh3  = 0.0;
  o->sgh4  = 0.0; o->sh2   = 0.0; o->sh3   = 0.0;
  o->si2   = 0.0; o->si3   = 0.0; o->sl2   = 0.0;
  o->sl3   = 0.0; o->sl4   = 0.0; o->gsto  = 0.0;
  o->xfact = 0.0; o->xgh2  = 0.0; o->xgh3  = 0.0;
  o->xgh4  = 0.0; o->xh2   = 0.0; o->xh3   = 0.0;
  o->xi2   = 0.0; o->xi3   = 0.0; o->xl2   = 0.0;
  o->xl3   = 0.0; o->xl4   = 0.0; o->xlamo = 0.0;
  o->zmol  = 0.0; o->zmos  = 0.0; o->atime = 0.0;
  o->xli   = 0.0; o->xni   = 0.0;

  // sgp4fix - note the following variables are also passed directly via o->
  // it is possible to streamline the sgp4init call by deleting the "x"
  // variables, but the user would need to set the o->* values first. we
  // include the additional assignments in case twoline2rv is not used.
  o->bstar   = t->bstar;
  o->ecco    = t->eccentricity;
  o->argpo   = t->arg_of_perigee;
  o->inclo   = t->inclination;
  o->mo      = t->mean_anomaly;
  o->no      = t->mean_motion;
  o->nodeo   = t->right_asc;

  // sgp4fix add opsmode
  o->operationmode = 'i';

  /* ------------------------ earth constants ----------------------- */
  // sgp4fix identify constants and allow alternate values
  getgravconst( whichconst, &tumin, &mu, &radiusearthkm, &xke, &j2, &j3, &j4, &j3oj2 );
  ss     = 78.0 / radiusearthkm + 1.0;
  // sgp4fix use multiply for speed instead of pow
  qzms2ttemp = (120.0 - 78.0) / radiusearthkm;
  qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
  x2o3   =  2.0 / 3.0;

  o->init = 'y';
  o->t   = 0.0;

  initl
      (
        satn, whichconst, o->ecco, epoch, o->inclo, &o->no, &o->method,
        &ainv, &ao, &o->con41, &con42, &cosio, &cosio2, &eccsq, &omeosq,
        &posq, &rp, &rteosq, &sinio, &o->gsto, o->operationmode
      );
  o->error = 0;

  // sgp4fix remove this check as it is unnecessary
  // the mrt check in sgp4 handles decaying satellite cases even if the starting
  // condition is below the surface of te earth
//     if (rp < 1.0)
//       {
//         printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
//         o->error = 5;
//       }

  if ((omeosq >= 0.0 ) || ( o->no >= 0.0))
    {
      o->isimp = 0;
      if (rp < (220.0 / radiusearthkm + 1.0))
          o->isimp = 1;
      sfour  = ss;
      qzms24 = qzms2t;
      perige = (rp - 1.0) * radiusearthkm;

      /* - for perigees below 156 km, s and qoms2t are altered - */
      if (perige < 156.0)
        {
          sfour = perige - 78.0;
          if (perige < 98.0)
              sfour = 20.0;
          // sgp4fix use multiply for speed instead of pow
          qzms24temp =  (120.0 - sfour) / radiusearthkm;
          qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
          sfour  = sfour / radiusearthkm + 1.0;
        }
      pinvsq = 1.0 / posq;

      tsi  = 1.0 / (ao - sfour);
      o->eta  = ao * o->ecco * tsi;
      etasq = o->eta * o->eta;
      eeta  = o->ecco * o->eta;
      psisq = fabs(1.0 - etasq);
      coef  = qzms24 * pow(tsi, 4.0);
      coef1 = coef / pow(psisq, 3.5);
      cc2   = coef1 * o->no * (ao * (1.0 + 1.5 * etasq + eeta *
                     (4.0 + etasq)) + 0.375 * j2 * tsi / psisq * o->con41 *
                     (8.0 + 3.0 * etasq * (8.0 + etasq)));
      o->cc1   = o->bstar * cc2;
      cc3   = 0.0;
      if (o->ecco > 1.0e-4)
          cc3 = -2.0 * coef * tsi * j3oj2 * o->no * sinio / o->ecco;
      o->x1mth2 = 1.0 - cosio2;
      o->cc4    = 2.0* o->no * coef1 * ao * omeosq *
                        (o->eta * (2.0 + 0.5 * etasq) + o->ecco *
                        (0.5 + 2.0 * etasq) - j2 * tsi / (ao * psisq) *
                        (-3.0 * o->con41 * (1.0 - 2.0 * eeta + etasq *
                        (1.5 - 0.5 * eeta)) + 0.75 * o->x1mth2 *
                        (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * o->argpo)));
      o->cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
                     (etasq + eeta) + eeta * etasq);
      cosio4 = cosio2 * cosio2;
      temp1  = 1.5 * j2 * pinvsq * o->no;
      temp2  = 0.5 * temp1 * j2 * pinvsq;
      temp3  = -0.46875 * j4 * pinvsq * pinvsq * o->no;
      o->mdot     = o->no + 0.5 * temp1 * rteosq * o->con41 + 0.0625 *
                         temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
      o->argpdot  = -0.5 * temp1 * con42 + 0.0625 * temp2 *
                          (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                          temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
      xhdot1            = -temp1 * cosio;
      o->nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
                           2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
      xpidot            =  o->argpdot+ o->nodedot;
      o->omgcof   = o->bstar * cc3 * cos(o->argpo);
      o->xmcof    = 0.0;
      if (o->ecco > 1.0e-4)
          o->xmcof = -x2o3 * coef * o->bstar / eeta;
      o->nodecf = 3.5 * omeosq * xhdot1 * o->cc1;
      o->t2cof   = 1.5 * o->cc1;
      // sgp4fix for divide by zero with xinco = 180 deg
      // sgp4fix divisor for divide by zero check on inclination
      // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
      // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
      const double temp4    =   1.5e-12;
      if (fabs(cosio+1.0) > 1.5e-12)
          o->xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
        else
          o->xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
      o->aycof   = -0.5 * j3oj2 * sinio;
      // sgp4fix use multiply for speed instead of pow
      delmotemp = 1.0 + o->eta * cos(o->mo);
      o->delmo   = delmotemp * delmotemp * delmotemp;
      o->sinmao  = sin(o->mo);
      o->x7thm1  = 7.0 * cosio2 - 1.0;

      /* --------------- deep space initialization ------------- */
      if ((2*pi / o->no) >= 225.0)
        {
          o->method = 'd';
          o->isimp  = 1;
          tc    =  0.0;
          inclm = o->inclo;

          dscom
              (
                epoch, o->ecco, o->argpo, tc, o->inclo, o->nodeo,
                o->no, &snodm, &cnodm,  &sinim, &cosim,&sinomm,     &cosomm,
                &day, &o->e3, &o->ee2, &em,         &emsq, &gam,
                &o->peo,  &o->pgho,   &o->pho, &o->pinco,
                &o->plo,  &rtemsq,        &o->se2, &o->se3,
                &o->sgh2, &o->sgh3,   &o->sgh4,
                &o->sh2,  &o->sh3,    &o->si2, &o->si3,
                &o->sl2,  &o->sl3,    &o->sl4, &s1, &s2, &s3, &s4, &s5,
                &s6,   &s7,   &ss1,  &ss2,  &ss3,  &ss4,  &ss5,  &ss6,  &ss7, &sz1, &sz2, &sz3,
                &sz11, &sz12, &sz13, &sz21, &sz22, &sz23, &sz31, &sz32, &sz33,
                &o->xgh2, &o->xgh3,   &o->xgh4, &o->xh2,
                &o->xh3,  &o->xi2,    &o->xi3,  &o->xl2,
                &o->xl3,  &o->xl4,    &nm, &z1, &z2, &z3, &z11,
                &z12, &z13, &z21, &z22, &z23, &z31, &z32, &z33,
                &o->zmol, &o->zmos
              );
          dpper
              (
                o->e3, o->ee2, o->peo, o->pgho,
                o->pho, o->pinco, o->plo, o->se2,
                o->se3, o->sgh2, o->sgh3, o->sgh4,
                o->sh2, o->sh3, o->si2,  o->si3,
                o->sl2, o->sl3, o->sl4,  o->t,
                o->xgh2,o->xgh3,o->xgh4, o->xh2,
                o->xh3, o->xi2, o->xi3,  o->xl2,
                o->xl3, o->xl4, o->zmol, o->zmos, inclm, o->init,
                &o->ecco, &o->inclo, &o->nodeo, &o->argpo, &o->mo,
                o->operationmode
              );

          argpm  = 0.0;
          nodem  = 0.0;
          mm     = 0.0;

          dsinit
              (
                whichconst,
                cosim, emsq, o->argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
                ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, o->t, tc,
                o->gsto, o->mo, o->mdot, o->no, o->nodeo,
                o->nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
                o->ecco, eccsq, &em, &argpm, &inclm, &mm, &nm, &nodem,
                &o->irez,  &o->atime,
                &o->d2201, &o->d2211, &o->d3210, &o->d3222 ,
                &o->d4410, &o->d4422, &o->d5220, &o->d5232,
                &o->d5421, &o->d5433, &o->dedt,  &o->didt,
                &o->dmdt,  &dndt,         &o->dnodt, &o->domdt ,
                &o->del1,  &o->del2,  &o->del3,  &o->xfact,
                &o->xlamo, &o->xli,   &o->xni
              );
        }

    /* ----------- set variables if not deep space ----------- */
    if (o->isimp != 1)
      {
        cc1sq          = o->cc1 * o->cc1;
        o->d2    = 4.0 * ao * tsi * cc1sq;
        temp           = o->d2 * tsi * o->cc1 / 3.0;
        o->d3    = (17.0 * ao + sfour) * temp;
        o->d4    = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) *
                         o->cc1;
        o->t3cof = o->d2 + 2.0 * cc1sq;
        o->t4cof = 0.25 * (3.0 * o->d3 + o->cc1 *
                         (12.0 * o->d2 + 10.0 * cc1sq));
        o->t5cof = 0.2 * (3.0 * o->d4 +
                         12.0 * o->cc1 * o->d3 +
                         6.0 * o->d2 * o->d2 +
                         15.0 * cc1sq * (2.0 * o->d2 + cc1sq));
      }
    } // if omeosq = 0 ...

    /* finally propogate to zero epoch to initialize all others. */
    // sgp4fix take out check to let satellites process until they are actually below earth surface
//       if(o->error == 0)
    sgp4(whichconst, o, 0.0, r, v);

    o->init = 'n';

//#include "debug6.cpp"
    //sgp4fix return boolean. satrec.error contains any error codes
    return true;
  return 0;
}

/*
 * Get instantaneous azimuth and elevation of a given satellite from the
 * ground based observer.
 */
int
azel(void)
{
  return 0;
}
