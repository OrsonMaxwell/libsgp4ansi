/*
 * libsgp4ansi.c - an ANSI C-11 SGP4/SDP4 implementation library for sgp4ansid.
 *
 * References:
 * https://www.celestrak.com/NORAD/documentation/spacetrk.pdf
 * https://celestrak.com/publications/AIAA/2006-6753/
 * IERS Bulletin - A (Vol. XXVIII No. 030)
 *
 * Copyright © 2017 Orson J. Maxwell. Please see LICENSE for details.
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "libsgp4ansi.h"
#include "const.h"
#include "transform.h"

// ************************************************************************* //
//                            PRIVATE PROTOTYPES                             //
// ************************************************************************* //

// SGP4 propagation function implementation
int
orbit_sgp4(orbit*, double, unsigned int, double, vect*, vect*);

// SDP4 propagation function implementation
int
orbit_sdp4(orbit*, double, unsigned int, double, vect*, vect*);

// Deep space long period periodics contributions to the mean elements
void
orbit_dslongper(orbit*, double, double*, double*, double*, double*, double*);

// ************************************************************************* //
//                             PRIVATE FUNCTIONS                             //
// ************************************************************************* //

/*
 * SGP4 propagation function implementation
 *
 * Inputs:  sat       - orbit struct pointer with initialized orbital data
 *          tdelta    - Time since epoch, minutes
 *          maxiter   - Kepler's equation maximum iteration count
 *          tolerance - Kepler's equation desired precision tolerance
 * Outputs: pos       - 3D position vector in TEME frame in km
 *          vel       - 3D velocity vector in TEME frame in km/sec
 * Returns: 0         - Success
 *          1         - Mean motion zero or less
 *          2         - Nonsensical orbital eccentricity (e >= 1; e < -0.00001)
 *          3         - Long periodics result error
 *          4         - Short period preliminary quantities error
 *          5         - Decaying satellite
 */
int
orbit_sgp4
(
    orbit* sat,
    double tdelta,
    unsigned int maxiter,
    double tolerance,
    vect* pos,
    vect* vel
)
{
  // Update for secular gravity and atmospheric drag
  double xmdf    = sat->Mo + sat->mdot * tdelta;
  double argpdf  = sat->omega + sat->omegaprime * tdelta;
  double nodedf  = sat->alpha + sat->nodedot * tdelta;
  double argpm   = argpdf;
  double mm      = xmdf;
  double t2      = tdelta * tdelta;
  double nodem   = nodedf + sat->nodecf * t2;
  double tempa   = 1.0 - sat->C1 * tdelta;
  double tempe   = sat->Bstar * sat->C4 * tdelta;
  double templ   = sat->t2cof * t2;

  double delomg, delmtemp, delm, temp, t3, t4;

  if (sat->islowperigee != true)
  {
    delomg   = sat->omgcof * tdelta;
    delmtemp =  1.0 + sat->eta * cos(xmdf);
    delm     = sat->xmcof * (delmtemp * delmtemp * delmtemp - sat->delMo);
    temp     = delomg + delm;
    mm       = xmdf + temp;
    argpm    = argpdf - temp;
    t3       = t2 * tdelta;
    t4       = t3 * tdelta;
    tempa    = tempa - sat->d2 * t2 - sat->d3 * t3 - sat->d4 * t4;
    tempe    = tempe + sat->Bstar * sat->C5 * (sin(mm) - sat->sinMo);
    templ   = templ + sat->t3cof * t3 + t4 * (sat->t4cof + tdelta * sat->t5cof);
  }

  if (sat->no <= 0.0)
  {
    return 1;
  }

  double em = sat->e;
  double am = pow((xke / sat->no), 2.0L / 3.0L) * tempa * tempa;
  sat->no   = xke / pow(am, 1.5);
  em        = em - tempe;

  // fix tolerance for error recognition
  if ((em >= 1.0) || (em < -0.00001))
  {
    return 2;
  }

  // Tolerance to avoid a division by zero
  if (em < 1.0e-12)
  {
    em = 1.0e-12;
  }

         mm   = mm + sat->no * templ;
  double xlm  = mm + argpm + nodem;
  double emsq = em * em;
         temp = 1.0 - emsq;

  nodem  = fmod(nodem, twopi);
  argpm  = fmod(argpm, twopi);
  xlm    = fmod(xlm, twopi);
  mm     = fmod(xlm - argpm - nodem, twopi);

  // Compute extra mean quantities
  double sinim = sin(sat->i);
  double cosim = cos(sat->i);

  double axnl = em * cos(argpm);
         temp = 1.0 / (am * (1.0 - em * em));
  double aynl = em* sin(argpm) + temp * sat->aycof;
  double xl   = mm + argpm + nodem + temp * sat->xlcof * axnl;

  // Kepler's equation
  double       u     = fmod(xl - nodem, twopi);
  double       eo1   = u;
  double       temp5 = 9999.9;
  unsigned int kiter = 0;

  double sineo1, coseo1;

  while ((fabs(temp5) >= tolerance) && (kiter < maxiter))
  {
    sineo1  = sin(eo1);
    coseo1  = cos(eo1);
    temp5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
    temp5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / temp5;

    if (fabs(temp5) >= 0.95)
    {
      temp5 = temp5 > 0.0 ? 0.95 : -0.95;
    }

    eo1     = eo1 + temp5;
    kiter++;
  }

  // Short period preliminary quantities
  double ecose     = axnl * coseo1 + aynl * sineo1;
  double esine     = axnl * sineo1 - aynl * coseo1;
  double el2       = axnl * axnl + aynl * aynl;
  double pl        = am * (1.0 - el2);
  double mrt       = 0.0;
  double vkmpersec = Re * xke / 60.0;

  double betal, rl,    rdotl, rvdotl, sinu, cosu,  su,    sin2u, cos2u,
         temp1, temp2, xnode, xinc,   mvt,  rvdot, sinsu, cossu, snod,
         cnod,  sini,  cosi,  xmx,    xmy,  ux,    uy,    uz,    vx,  vy,  vz;

  if (pl < 0.0)
  {
    return 4;
  }
  else
  {
    rl     = am * (1.0 - ecose);
    rdotl  = sqrt(am) * esine/rl;
    rvdotl = sqrt(pl) / rl;
    betal  = sqrt(1.0 - el2);
    temp   = esine / (1.0 + betal);
    sinu   = am / rl * (sineo1 - aynl - axnl * temp);
    cosu   = am / rl * (coseo1 - axnl + aynl * temp);
    su     = atan2(sinu, cosu);
    sin2u  = (cosu + cosu) * sinu;
    cos2u  = 1.0 - 2.0 * sinu * sinu;
    temp   = 1.0 / pl;
    temp1  = 0.5 * j2 * temp;
    temp2  = temp1 * temp;
    mrt    = rl * (1.0 - 1.5 * temp2 * betal * sat->con41) + 0.5 * temp1
             * sat->x1mth2 * cos2u;
    su     = su - 0.25 * temp2 * sat->x7thm1 * sin2u;
    xnode  = nodem + 1.5 * temp2 * cosim * sin2u;
    xinc   = sat->i + 1.5 * temp2 * cosim * sinim * cos2u;
    mvt    = rdotl - sat->no * temp1 * sat->x1mth2 * sin2u / xke;
    rvdot  = rvdotl + sat->no * temp1 * (sat->x1mth2 * cos2u + 1.5
             * sat->con41) / xke;

    // Orientation vectors
    sinsu =  sin(su);
    cossu =  cos(su);
    snod  =  sin(xnode);
    cnod  =  cos(xnode);
    sini  =  sin(xinc);
    cosi  =  cos(xinc);
    xmx   = -snod * cosi;
    xmy   =  cnod * cosi;
    ux    =  xmx * sinsu + cnod * cossu;
    uy    =  xmy * sinsu + snod * cossu;
    uz    =  sini * sinsu;
    vx    =  xmx * cossu - cnod * sinsu;
    vy    =  xmy * cossu - snod * sinsu;
    vz    =  sini * cossu;

    // Position and velocity vectors
    pos->x = (mrt * ux) * Re;
    pos->y = (mrt * uy) * Re;
    pos->z = (mrt * uz) * Re;
    vel->x = (mvt * ux + rvdot * vx) * vkmpersec;
    vel->y = (mvt * uy + rvdot * vy) * vkmpersec;
    vel->z = (mvt * uz + rvdot * vz) * vkmpersec;
  }

  // Decaying satellites
  if (mrt < 1.0)
  {
    return 5;
  }

  return 0;
}

/*
 * SDP4 propagation function implementation
 *
 * Inputs:  sat       - orbit struct pointer with initialized orbital data
 *          tdelta    - Time since epoch, minutes
 *          maxiter   - Kepler's equation maximum iteration count
 *          tolerance - Kepler's equation desired precision threshold
 * Outputs: pos       - 3D position vector in TEME frame in km
 *          vel       - 3D velocity vector in TEME frame in km/sec
 * Returns: 0         - Success
 *          1         - Mean motion zero or less
 *          2         - Nonsensical orbital eccentricity (e >= 1; e < -0.00001)
 *          3         - Long periodics result error
 *          4         - Short period preliminary quantities error
 *          5         - Decaying satellite
 */
int
orbit_sdp4
(
    orbit* sat,
    double tdelta,
    unsigned int iter,
    double thresh,
    vect* pos,
    vect* vel
)
{
  // Update for secular gravity and atmospheric drag
  double xmdf    = sat->Mo + sat->mdot * tdelta;
  double argpdf  = sat->omega + sat->omegaprime * tdelta;
  double nodedf  = sat->alpha + sat->nodedot * tdelta;
  double argpm   = argpdf;
  double mm      = xmdf;
  double t2      = tdelta * tdelta;
  double nodem   = nodedf + sat->nodecf * t2;
  double tempa   = 1.0 - sat->C1 * tdelta;
  double tempe   = sat->Bstar * sat->C4 * tdelta;
  double templ   = sat->t2cof * t2;

  double delomg, delmtemp, delm, temp, t3, t4;

  if (sat->islowperigee != true)
  {
    delomg   = sat->omgcof * tdelta;
    delmtemp =  1.0 + sat->eta * cos(xmdf);
    delm     = sat->xmcof * (delmtemp * delmtemp * delmtemp - sat->delMo);
    temp     = delomg + delm;
    mm       = xmdf + temp;
    argpm    = argpdf - temp;
    t3       = t2 * tdelta;
    t4       = t3 * tdelta;
    tempa    = tempa - sat->d2 * t2 - sat->d3 * t3 - sat->d4 * t4;
    tempe    = tempe + sat->Bstar * sat->C5 * (sin(mm) - sat->sinMo);
    templ   = templ + sat->t3cof * t3 + t4 * (sat->t4cof + tdelta * sat->t5cof);
  }

  double nm    = sat->no;
  double em    = sat->e;
  double inclm = sat->i;
  if (sat->isdeepspace == true)
  {
    double tc = tdelta;
    /*dspace
    (
        sat->irez,
        sat->d2201, sat->d2211, sat->d3210,
        sat->d3222, sat->d4410, sat->d4422,
        sat->d5220, sat->d5232, sat->d5421,
        sat->d5433, sat->dedt,  sat->del1,
        sat->del2,  sat->del3,  sat->didt,
        sat->dmdt,  sat->dnodt, sat->domdt,
        sat->argpo, sat->argpdot, sat->t, tc,
        sat->gsto, sat->xfact, sat->xlamo,
        sat->no, &sat->atime,
        &em, &argpm, &inclm, &sat->xli, &mm, &sat->xni,
        &nodem, &dndt, &nm
    );*/
  }

  if (nm <= 0.0)
  {
    return 1;
  }

  double am = pow((xke / nm), 2.0L / 3.0L) * tempa * tempa;
  nm = xke / pow(am, 1.5);
  em = em - tempe;

  // fix tolerance for error recognition
  if ((em >= 1.0) || (em < -0.00001))
  {
    return 2;
  }

  // Tolerance to avoid a division by zero
  if (em < 1.0e-6)
    em  = 1.0e-6;

  mm     = mm + sat->no * templ;
  double xlm    = mm + argpm + nodem;
  double emsq   = em * em;
  temp   = 1.0 - emsq;

  nodem  = fmod(nodem, twopi);
  argpm  = fmod(argpm, twopi);
  xlm    = fmod(xlm, twopi);
  mm     = fmod(xlm - argpm - nodem, twopi);

  // Compute extra mean quantities
  double sinim = sin(inclm);
  double cosim = cos(inclm);

  // Add lunar-solar
  double ep     = em;
  double xincp  = inclm;
  double argpp  = argpm;
  double nodep  = nodem;
  double mp     = mm;
  double sinip  = sinim;
  double cosip  = cosim;
  if (sat->isdeepspace == true)
  {
    /*dpper
    (
        sat->e3,   sat->ee2,  sat->peo,
        sat->pgho, sat->pho,  sat->pinco,
        sat->plo,  sat->se2,  sat->se3,
        sat->sgh2, sat->sgh3, sat->sgh4,
        sat->sh2,  sat->sh3,  sat->si2,
        sat->si3,  sat->sl2,  sat->sl3,
        sat->sl4,  tdelta,    sat->xgh2,
        sat->xgh3, sat->xgh4, sat->xh2,
        sat->xh3,  sat->xi2,  sat->xi3,
        sat->xl2,  sat->xl3,  sat->xl4,
        sat->zmol, sat->zmos, sat->i,
        'n', &ep, &xincp, &nodep, &argpp, &mp, 'i'
    );*/
    if (xincp < 0.0)
    {
      xincp  = -xincp;
      nodep = nodep + pi;
      argpp  = argpp - pi;
    }

    if ((ep < 0.0 ) || ( ep > 1.0))
    {
      return 3;
    }
  }

  // Long period periodics
  if (sat->isdeepspace == true)
  {
    sinip =  sin(xincp);
    cosip =  cos(xincp);
    sat->aycof = -0.5*j3divj2*sinip;

    // Division by zero check for xincp = 180 deg
    if (fabs(cosip+1.0) > 1.5e-12)
      sat->xlcof = -0.25 * j3divj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
    else
      sat->xlcof = -0.25 * j3divj2 * sinip * (3.0 + 5.0 * cosip) / 1.5e-12;
  }

  double axnl = ep * cos(argpp);
  temp = 1.0 / (am * (1.0 - ep * ep));
  double aynl = ep* sin(argpp) + temp * sat->aycof;
  double xl   = mp + argpp + nodep + temp * sat->xlcof * axnl;

  // Kepler's equation
  double u    = fmod(xl - nodep, twopi);
  double eo1  = u;
  double temp5 = 9999.9;
  unsigned int ktr = 0;

  double sineo1, coseo1;

  while ((fabs(temp5) >= thresh) && (ktr < iter))
  {
    sineo1 = sin(eo1);
    coseo1 = cos(eo1);
    temp5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
    temp5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / temp5;
    if (fabs(temp5) >= 0.95)
      temp5 = temp5 > 0.0 ? 0.95 : -0.95;
    eo1    = eo1 + temp5;
    ktr++;
  }

  // Short period preliminary quantities
  double ecose = axnl * coseo1 + aynl * sineo1;
  double esine = axnl * sineo1 - aynl * coseo1;
  double el2   = axnl * axnl + aynl * aynl;
  double pl    = am * (1.0 - el2);
  double mrt = 0.0;
  double vkmpersec = Re * xke / 60.0;

  double betal, rl, rdotl, rvdotl, sinu, cosu, su, sin2u, cos2u, temp1, temp2,
  cosisq, xnode, xinc, mvt, rvdot, sinsu, cossu, snod, cnod, sini, cosi, xmx,
  xmy, ux, uy, uz, vx, vy, vz;

  if (pl < 0.0)
  {
    return 4;
  }
  else
  {
    rl     = am * (1.0 - ecose);
    rdotl  = sqrt(am) * esine/rl;
    rvdotl = sqrt(pl) / rl;
    betal  = sqrt(1.0 - el2);
    temp   = esine / (1.0 + betal);
    sinu   = am / rl * (sineo1 - aynl - axnl * temp);
    cosu   = am / rl * (coseo1 - axnl + aynl * temp);
    su     = atan2(sinu, cosu);
    sin2u  = (cosu + cosu) * sinu;
    cos2u  = 1.0 - 2.0 * sinu * sinu;
    temp   = 1.0 / pl;
    temp1  = 0.5 * j2 * temp;
    temp2  = temp1 * temp;

    // Update for short period periodics
    if (sat->isdeepspace == true)
    {
      cosisq      = cosip * cosip;
      sat->con41  = 3.0*cosisq - 1.0;
      sat->x1mth2 = 1.0 - cosisq;
      sat->x7thm1 = 7.0*cosisq - 1.0;
    }
    mrt   = rl * (1.0 - 1.5 * temp2 * betal * sat->con41) +
        0.5 * temp1 * sat->x1mth2 * cos2u;
    su    = su - 0.25 * temp2 * sat->x7thm1 * sin2u;
    xnode = nodep + 1.5 * temp2 * cosip * sin2u;
    xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
    mvt   = rdotl - nm * temp1 * sat->x1mth2 * sin2u / xke;
    rvdot = rvdotl + nm * temp1 * (sat->x1mth2 * cos2u +
        1.5 * sat->con41) / xke;

    // Orientation vectors
    sinsu =  sin(su);
    cossu =  cos(su);
    snod  =  sin(xnode);
    cnod  =  cos(xnode);
    sini  =  sin(xinc);
    cosi  =  cos(xinc);
    xmx   = -snod * cosi;
    xmy   =  cnod * cosi;
    ux    =  xmx * sinsu + cnod * cossu;
    uy    =  xmy * sinsu + snod * cossu;
    uz    =  sini * sinsu;
    vx    =  xmx * cossu - cnod * sinsu;
    vy    =  xmy * cossu - snod * sinsu;
    vz    =  sini * cossu;

    // Position and velocity vectors
    pos->x = (mrt * ux)* Re;
    pos->y = (mrt * uy)* Re;
    pos->z = (mrt * uz)* Re;
    vel->x = (mvt * ux + rvdot * vx) * vkmpersec;
    vel->y = (mvt * uy + rvdot * vy) * vkmpersec;
    vel->z = (mvt * uz + rvdot * vz) * vkmpersec;
  }

  // Decaying satellites
  if (mrt < 1.0)
  {
    return 5;
  }

  return 0;
}

/*
 * Deep space long period periodics contributions to the mean elements
 *
 * Inputs:  sat    - orbit struct pointer with full orbital data
 *          tdelta - Time since epoch, minutes
 * Outputs: sat    - orbit struct pointer with long period periodics applied
 *          e      - Eccentricity with lpp applied
 *          i      - Inclination with lpp applied
 *          alpha  - Right ascension of ascension node with lpp applied
 *          omega  - Argument of perigee with lpp applied
 *          mp     - Mean motion term with lpp applied
 * Returns: None
 *
 * If any of the output pointers are NULL, assume initial run and apply lpp to
 * satellite data directly
 */
void
orbit_dslongper
(
    orbit* sat,
    double tdelta,
    double* e,
    double* i,
    double* alpha,
    double* omega,
    double* mp
)
{
  bool init = false;
  if((e == NULL)||(i == NULL)||(alpha == NULL)||(omega == NULL)||(mp == NULL))
  {
    init = true;
    tdelta = 0;
  }

  double alfdp, betdp, cosip, cosop, dalf, dbet, dls,
  f2,    f3,    pe,    pgh,   ph,   pinc, pl ,
  sel,   ses,   sghl,  sghs,  shll, shs,  sil,
  sinip, sinop, sinzf, sis,   sll,  sls,  xls,
  xnoh,  zf,    zm,    zel,   zes,  znl,  zns;

  // Constants
  zns   = 1.19459e-5;
  zes   = 0.01675;
  znl   = 1.5835218e-4;
  zel   = 0.05490;

  zm = sat->zmos + zns * tdelta;

  zf    = zm + 2.0 * zes * sin(zm);
  sinzf = sin(zf);
  f2    =  0.5 * sinzf * sinzf - 0.25;
  f3    = -0.5 * sinzf * cos(zf);
  ses   = sat->se2* f2 + sat->se3 * f3;
  sis   = sat->si2 * f2 + sat->si3 * f3;
  sls   = sat->sl2 * f2 + sat->sl3 * f3 + sat->sl4 * sinzf;
  sghs  = sat->sgh2 * f2 + sat->sgh3 * f3 + sat->sgh4 * sinzf;
  shs   = sat->sh2 * f2 + sat->sh3 * f3;

  zm    = sat->zmol + znl * tdelta;

  zf    = zm + 2.0 * zel * sin(zm);
  sinzf = sin(zf);
  f2    =  0.5 * sinzf * sinzf - 0.25;
  f3    = -0.5 * sinzf * cos(zf);
  sel   = sat->ee2 * f2 + sat->e3 * f3;
  sil   = sat->xi2 * f2 + sat->xi3 * f3;
  sll   = sat->xl2 * f2 + sat->xl3 * f3 + sat->xl4 * sinzf;
  sghl  = sat->xgh2 * f2 + sat->xgh3 * f3 + sat->xgh4 * sinzf;
  shll  = sat->xh2 * f2 +sat-> xh3 * f3;
  pe    = ses + sel;
  pinc  = sis + sil;
  pl    = sls + sll;
  pgh   = sghs + sghl;
  ph    = shs + shll;

  if (init == false)
  {
    pe    = pe - sat->peo;
    pinc  = pinc - sat->pinco;
    pl    = pl - sat->plo;
    pgh   = pgh - sat->pgho;
    ph    = ph - sat->pho;
    *i    = *i + pinc;
    *e    = *e + pe;
    sinip = sin(*i);
    cosip = cos(*i);

    // Apply periodics directly accounting for Lyddane choice
    //  Readjust the 0.2 limit value and limit discontinuity
    //  0.2 rad = 11.45916 deg
    if (*i >= 0.2)
    {
      ph     = ph / sinip;
      pgh    = pgh - cosip * ph;
      *omega  = *omega + pgh;
      *alpha  = *alpha + ph;
      *mp     = *mp + pl;
    }
    else
    {
      /* ---- apply periodics with lyddane modification ---- */
      sinop  = sin(*alpha);
      cosop  = cos(*alpha);
      alfdp  = sinip * sinop;
      betdp  = sinip * cosop;
      dalf   =  ph * cosop + pinc * cosip * sinop;
      dbet   = -ph * sinop + pinc * cosip * cosop;
      alfdp  = alfdp + dalf;
      betdp  = betdp + dbet;
      *alpha  = fmod(*alpha, twopi);

      xls    = *mp + *omega + cosip * (*alpha);
      dls    = pl + pgh - pinc * (*alpha) * sinip;
      xls    = xls + dls;
      xnoh   = *alpha;
      *alpha  = atan2(alfdp, betdp);

      if (fabs(xnoh - *alpha) > pi)
      {
        if (*alpha < xnoh)
        {
          *alpha = *alpha + twopi;
        }
        else
        {
          *alpha = *alpha - twopi;
        }
      }

      *mp    = *mp + pl;
      *omega = xls - *mp - cosip * (*alpha);
    }
  }
}

// ************************************************************************* //
//                                 INTERFACE                                 //
// ************************************************************************* //

/*
 * Initialize SGP4/SDP4 orbit model from a raw NORAD TLE lines
 *
 * Inputs:  sat  - Empty orbit struct pointer with TLE data
 * Outputs: sat  - orbit struct pointer with full orbital data
 * Returns: None
 *
 * Calls: orbit_init
 */
int
tle2orbit(char* tlestr1, char* tlestr2, orbit* sat)
{
  // Pre-formatting the raw TLE input
  for (char j = 10; j <= 15; j++)
    if (tlestr1[j] == ' ')
      tlestr1[j] = '_';

  if (tlestr1[44] != ' ')
    tlestr1[43] = tlestr1[44];

  tlestr1[44] = '.';
  if (tlestr1[7] == ' ')
    tlestr1[7] = 'U';

  if (tlestr1[9] == ' ')
    tlestr1[9] = '.';

  for (char j = 45; j <= 49; j++)
    if (tlestr1[j] == ' ')
      tlestr1[j] = '0';

  if (tlestr1[51] == ' ')
    tlestr1[51] = '0';

  if (tlestr1[53] != ' ')
    tlestr1[52] = tlestr1[53];

  tlestr1[53] = '.';
  tlestr2[25] = '.';

  for (char j = 26; j <= 32; j++)
    if (tlestr2[j] == ' ')
      tlestr2[j] = '0';

  if (tlestr1[62] == ' ')
    tlestr1[62] = '0';

  if (tlestr1[68] == ' ')
    tlestr1[68] = '0';

  struct tm epoch_tm, prop_tm;

  int cardnum, epochyr, epochdays, nexp, Bexp;
  double nddot, Bstar;

  sscanf(tlestr1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",
         &cardnum, &sat->number, &sat->sec_class, sat->designator, &epochyr,
         &epochdays,&sat->nprimediv2, &nddot, &nexp, &Bstar,
         &Bexp, &sat->ephem_type, &sat->elset_number);

  fractday2unix(epochyr, epochdays, &sat->epoch);

  sat->ndprimediv6 = nddot * pow(10, nexp);
  sat->Bstar = Bstar * pow(10, Bexp);

  if (tlestr2[52] == ' ') // check for minus sign
    sscanf(tlestr2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %10lf %6ld \n",
           &cardnum,&sat->number, &sat->i, &sat->alpha, &sat->e, &sat->omega,
           &sat->Mo, &sat->no, &sat->rev_number);
  else
    sscanf(tlestr2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld \n",
           &cardnum,&sat->number, &sat->i, &sat->alpha, &sat->e, &sat->omega,
           &sat->Mo, &sat->no, &sat->rev_number);

  // Convert to SGP4 units and expand
  return orbit_init(sat);
}

/*
 * Expand SGP4/SDP4 orbit elements from an orbit containing NORAD TLE portion
 *
 * Inputs:  sat  - orbit struct pointer to an empty orbit
 * Outputs: sat  - orbit struct pointer with full orbital data
 * Returns: None
 */
int
orbit_init(orbit* sat)
{
  // Convert to SGP4 units
  sat->no          = sat->no / rpd2radmin;
  sat->a           = pow(sat->no * tumin, (-2.0L / 3.0L));
  sat->nprimediv2  = sat->nprimediv2  / (rpd2radmin * 1440);
  sat->ndprimediv6 = sat->ndprimediv6 / (rpd2radmin * 1440 * 1440);
  sat->julepoch    = unix2jul(&sat->epoch, sat->epoch_ms);
  sat->GSTo        = jul2gst(sat->julepoch);

  // Standard orbital elements
  sat->i       = sat->i * deg2rad;
  sat->alpha   = sat->alpha * deg2rad;
  sat->omega   = sat->omega * deg2rad;
  sat->Mo      = sat->Mo * deg2rad;
  sat->altapoR = sat->a * (1.0 + sat->e) - 1.0;
  sat->altperR = sat->a * (1.0 - sat->e) - 1.0;
  sat->sini    = sin(sat->i);
  sat->cosi    = cos(sat->i);

  // Aux epoch quantities
  double esq     = pow(sat->e, 2);
  double omegasq = 1.0 - esq;
  double rteosq  = sqrt(omegasq);
  double cosi2   = pow(sat->cosi, 2);

  // Reverse Kozai notation of the mean motion
  double ak    = pow(xke / sat->no, 2.0L / 3.0L);
  double d1    = 0.75 * j2 * (3.0 * cosi2 - 1.0) / (rteosq * omegasq);
  double del   = d1 / pow(ak, 2);
  double adel  = ak * (1.0 - del * del - del * (1.0L / 3.0L + 134.0 * del
                 * del / 81.0));
  del          = d1 / (adel * adel);

  sat->no      = sat->no / (1.0 + del);

  double po    = sat->a * omegasq;
  double con42 = 1.0 - 5.0 * cosi2;

  sat->con41   = -con42 -2 * cosi2;

  double ainv  = 1.0 / sat->a;
  double posq  = pow(po, 2);

  if ((omegasq >= 0.0 ) || (sat->no >= 0.0))
  {
    double s4th   = 78.0 / Re + 1.0;
    double qzms24 = pow((120.0 - 78.0) / Re, 4);
    double altperkm = sat->altperR * Re;

    sat->islowperigee    = false;
    if (altperkm < 220.0)
    {
      sat->islowperigee  = true;
    }

    // For perigees below 156km, recompute s and qoms2t
    if (altperkm < 156.0)
    {
      s4th = altperkm - 78.0;
      if (altperkm < 98.0)
      {
        s4th   = 20.0;
      }

      qzms24 = pow((120.0 - s4th) / Re, 4);
      s4th   = s4th / Re + 1.0;
    }

    double pinvsq = 1.0 / posq;

    double tsi    = 1.0 / (sat->a - s4th);
    sat->eta      = sat->a * sat->e * tsi;

    double etasq  = pow(sat->eta, 2);
    double eeta   = sat->e * sat->eta;
    double psisq  = fabs(1.0 - etasq);
    double coef   = qzms24 * pow(tsi, 4.0);
    double coef1  = coef / pow(psisq, 3.5);
    double C2     = coef1 * sat->no * (sat->a * (1.0 + 1.5 * etasq + eeta *
                    (4.0 + etasq)) + 0.375 * j2 * tsi / psisq *
                    sat->con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));

    sat->C1       = sat->Bstar * C2;

    double C3 = 0.0;
    // Low eccentricity case
    if (sat->e > 1.0e-4)
    {
      C3 = -2.0 * coef * tsi * j3divj2 * sat->no * sat->sini / sat->e;
    }

    sat->x1mth2   = 1.0 - cosi2;
    sat->C4       = 2.0* sat->no * coef1 * sat->a * omegasq * (sat->eta * (2.0 + 0.5
                    * etasq) + sat->e * (0.5 + 2.0 * etasq) - j2 * tsi / (sat->a *
                    psisq) * (-3.0 * sat->con41 * (1.0 - 2.0 * eeta + etasq *
                    (1.5 - 0.5 * eeta)) + 0.75 * sat->x1mth2 * (2.0 * etasq -
                    eeta * (1.0 + etasq)) * cos(2.0 * sat->omega)));

    sat->C5       = 2.0 * coef1 * sat->a * omegasq * (1.0 + 2.75 * (etasq + eeta) +
                   eeta * etasq);

    double cosi4  = pow(cosi2, 2);
    double temp1  = 1.5 * j2 * pinvsq * sat->no;
    double temp2  = 0.5 * temp1 * j2 * pinvsq;
    double temp3  = -0.46875 * j4 * pinvsq * pinvsq * sat->no;

    sat->mdot     = sat->no + 0.5 * temp1 * rteosq * sat->con41 + 0.0625 * temp2
                    * rteosq * (13.0 - 78.0 * cosi2 + 137.0 * cosi4);

    sat->omegaprime = -0.5 * temp1 * con42 + 0.0625 * temp2 * (7.0 -
                      114.0 * cosi2 + 395.0 * cosi4) + temp3 * (3.0 - 36.0
                      * cosi2 + 49.0 * cosi4);

    double xhdot1     = -temp1 * sat->cosi;

    sat->nodedot= xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosi2) + 2.0 *
                        temp3 * (3.0 - 7.0 * cosi2)) * sat->cosi;

    double xpidot     = sat->omegaprime + sat->nodedot;
    double sinomega   = sin(sat->omega);
    double cosomega   = cos(sat->omega);

    sat->omgcof  = sat->Bstar * C3 * cosomega;
    sat->xmcof   = 0.0;

    if (sat->e > 1.0e-4)
      sat->xmcof = -(2.0L / 3.0L) * coef * sat->Bstar / eeta;

    sat->nodecf  = 3.5 * omegasq * xhdot1 * sat->C1;
    sat->t2cof   = 1.5 * sat->C1;

    // Division by zero protection in case xinco = 180deg
    if (fabs(sat->cosi + 1.0) > 1.5e-12)
      sat->xlcof = -0.25 * j3divj2 * sat->sini * (3.0 + 5.0 * sat->cosi) / (1.0 + sat->cosi);
    else
      sat->xlcof = -0.25 * j3divj2 * sat->sini * (3.0 + 5.0 * sat->cosi) / 1.5e-12;

    sat->aycof   = -0.5 * j3divj2 * sat->sini;

    sat->delMo   = pow(1.0 + sat->eta * cos(sat->Mo), 3);
    sat->sinMo   = sin(sat->Mo);
    sat->x7thm1  = 7.0 * cosi2 - 1.0;

    //*************************** DEEP SPACE ***********************************
    // TODO: Deepspace routines give me headaches
    // Local variables for deep space and resonance perturbations
    double z1,  z2,  z3,  z11,  z12,  z13,  z21,  z22,  z23,  z31,  z32,  z33;
    double sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32,sz33;
    double s1,  s2,  s3,  s4,  s5,  s6,  s7;
    double ss1, ss2, ss3, ss4, ss5, ss6, ss7;

    // Deep space initializations
    //if ((twopi / sat->no) >= 225.0)
    if (false)
    {
      sat->isdeepspace = true;
      sat->islowperigee       = true;

      // Constants
      double zes     =  0.01675;
      double zel     =  0.05490;
      double C1ss    =  2.9864797e-6;
      double C1l     =  4.7968065e-7;
      double zsinis  =  0.39785416;
      double zcosis  =  0.91744867;
      double zcosgs  =  0.1945905;
      double zsings  = -0.98088458;

      double sinalpha  = sin(sat->alpha);
      double cosalpha  = cos(sat->alpha);
      double betasq    = 1.0 - esq;
      double rtemsq    = sqrt(betasq);

      // Initialize lunar and solar terms
      double day    = sat->julepoch + 18261.5;
      double xnodce = fmod(4.5236020 - 9.2422029e-4 * day, twopi);
      double stem   = sin(xnodce);
      double ctem   = cos(xnodce);
      double zcosil = 0.91375164 - 0.03568096 * ctem;
      double zsinil = sqrt(1.0 - zcosil * zcosil);
      double zsinhl = 0.089683511 * stem / zsinil;
      double zcoshl = sqrt(1.0 - zsinhl * zsinhl);
      double gamma  = 5.8351514 + 0.0019443680 * day;
      double zx     = 0.39785416 * stem / zsinil;
      double zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
      zx            = atan2(zx, zy);
      zx            = gamma + zx - xnodce;
      double zcosgl = cos(zx);
      double zsingl = sin(zx);

      // Get solar terms
      // Local vars for forked iteration
      double zcosg = zcosgs;
      double zsing = zsings;
      double zcosi = zcosis;
      double zsini = zsinis;
      double zcosh = cosalpha;
      double zsinh = sinalpha;
      double Cc    = C1ss;
      double xnoi  = 1.0 / sat->no;

      double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
      double x1, x2, x3, x4, x5, x6, x7, x8;

      for (uint8_t iter = 1; iter <= 2; iter++)
      {
        a1  =   zcosg * zcosh + zsing * zcosi * zsinh;
        a3  =  -zsing * zcosh + zcosg * zcosi * zsinh;
        a7  =  -zcosg * zsinh + zsing * zcosi * zcosh;
        a8  =   zsing * zsini;
        a9  =   zsing * zsinh + zcosg * zcosi * zcosh;
        a10 =   zcosg * zsini;
        a2  =   sat->cosi * a7 + sat->sini * a8;
        a4  =   sat->cosi * a9 + sat->sini * a10;
        a5  =  -sat->sini * a7 + sat->cosi * a8;
        a6  =  -sat->sini * a9 + sat->cosi * a10;

        x1  =  a1 * cosomega + a2 * sinomega;
        x2  =  a3 * cosomega + a4 * sinomega;
        x3  = -a1 * sinomega + a2 * cosomega;
        x4  = -a3 * sinomega + a4 * cosomega;
        x5  =  a5 * sinomega;
        x6  =  a6 * sinomega;
        x7  =  a5 * cosomega;
        x8  =  a6 * cosomega;

        z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
        z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
        z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
        z1  =  3.0 *  (a1 * a1 + a2 * a2) + z31 * esq;
        z2  =  6.0 *  (a1 * a3 + a2 * a4) + z32 * esq;
        z3  =  3.0 *  (a3 * a3 + a4 * a4) + z33 * esq;
        z11 = -6.0 * a1 * a5 + esq *  (-24.0 * x1 * x7-6.0 * x3 * x5);
        z12 = -6.0 *  (a1 * a6 + a3 * a5) + esq *
            (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
        z13 = -6.0 * a3 * a6 + esq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
        z21 =  6.0 * a2 * a5 + esq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
        z22 =  6.0 *  (a4 * a5 + a2 * a6) + esq *
            (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
        z23 =  6.0 * a4 * a6 + esq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
        z1  = z1 + z1 + betasq * z31;
        z2  = z2 + z2 + betasq * z32;
        z3  = z3 + z3 + betasq * z33;
        s3  = Cc * xnoi;
        s2  = -0.5 * s3 / rtemsq;
        s4  = s3 * rtemsq;
        s1  = -15.0 * sat->e * s4;
        s5  = x1 * x3 + x2 * x4;
        s6  = x2 * x3 + x1 * x4;
        s7  = x2 * x4 - x1 * x3;

        // Get lunar terms
        if (iter == 1)
        {
          ss1   = s1;
          ss2   = s2;
          ss3   = s3;
          ss4   = s4;
          ss5   = s5;
          ss6   = s6;
          ss7   = s7;
          sz1   = z1;
          sz2   = z2;
          sz3   = z3;
          sz11  = z11;
          sz12  = z12;
          sz13  = z13;
          sz21  = z21;
          sz22  = z22;
          sz23  = z23;
          sz31  = z31;
          sz32  = z32;
          sz33  = z33;
          zcosg = zcosgl;
          zsing = zsingl;
          zcosi = zcosil;
          zsini = zsinil;
          zcosh = zcoshl * cosalpha + zsinhl * sinalpha;
          zsinh = sinalpha * zcoshl - cosalpha * zsinhl;
          Cc    = C1l;
        }
      }

      sat->zmol = fmod(4.7199672 + 0.22997150  * day - gamma, twopi);
      sat->zmos = fmod(6.2565837 + 0.017201977 * day, twopi);

      // Apply solar terms
      sat->se2  =   2.0 * ss1 * ss6;
      sat->se3  =   2.0 * ss1 * ss7;
      sat->si2  =   2.0 * ss2 * sz12;
      sat->si3  =   2.0 * ss2 * (sz13 - sz11);
      sat->sl2  =  -2.0 * ss3 * sz2;
      sat->sl3  =  -2.0 * ss3 * (sz3 - sz1);
      sat->sl4  =  -2.0 * ss3 * (-21.0 - 9.0 * esq) * zes;
      sat->sgh2 =   2.0 * ss4 * sz32;
      sat->sgh3 =   2.0 * ss4 * (sz33 - sz31);
      sat->sgh4 = -18.0 * ss4 * zes;
      sat->sh2  =  -2.0 * ss2 * sz22;
      sat->sh3  =  -2.0 * ss2 * (sz23 - sz21);

      // Apply lunar terms
      sat->ee2  =   2.0 * s1 * s6;
      sat->e3   =   2.0 * s1 * s7;
      sat->xi2  =   2.0 * s2 * z12;
      sat->xi3  =   2.0 * s2 * (z13 - z11);
      sat->xl2  =  -2.0 * s3 * z2;
      sat->xl3  =  -2.0 * s3 * (z3 - z1);
      sat->xl4  =  -2.0 * s3 * (-21.0 - 9.0 * esq) * zel;
      sat->xgh2 =   2.0 * s4 * z32;
      sat->xgh3 =   2.0 * s4 * (z33 - z31);
      sat->xgh4 = -18.0 * s4 * zel;
      sat->xh2  =  -2.0 * s2 * z22;
      sat->xh3  =  -2.0 * s2 * (z23 - z21);

      orbit_dslongper(sat, 0L, NULL, NULL, NULL, NULL, NULL);
    }

    if (sat->islowperigee != true)
    {
      double C1sq  = pow(sat->C1, 2);

      sat->d2     = 4.0 * sat->a * tsi * C1sq;

      double temp = sat->d2 * tsi * sat->C1 / 3.0;

      sat->d3     = (17.0 * sat->a + s4th) * temp;
      sat->d4     = 0.5 * temp * sat->a * tsi * (221.0 * sat->a + 31.0 * s4th)
                          * sat->C1;
      sat->t3cof  = sat->d2 + 2.0 * C1sq;
      sat->t4cof  = 0.25 * (3.0 * sat->d3 + sat->C1 * (12.0 * sat->d2 + 10.0
          * C1sq));
      sat->t5cof  = 0.2 * (3.0 * sat->d4 + 12.0 * sat->C1 * sat->d3 + 6.0 *
          sat->d2 * sat->d2 + 15.0 * C1sq * (2.0 * sat->d2 + C1sq));
    }
  }

  vect pos = {0}, vel = {0};

  if (false) // TODO: Deepspace FFS
  //if (sat->isdeepspace == true)
  {
    orbit_sdp4(sat, 0.0L, 5, 1.0e-12, &pos, &vel);
  }
  else
  {
    orbit_sgp4(sat, 0.0L, 5, 1.0e-12, &pos, &vel);
  }

  return 0;
}

/*
 * Get position and velocity vectors in the TEME frame at given time
 *
 * Inputs:  sat       - orbit struct pointer with initialized orbital data
 *          time      - UTC time to find satellite data at
 *          msec      - Millisecond portion of time
 *          maxiter   - Kepler's equation maximum iteration count
 *          tolerance - Kepler's equation desired precision tolerance
 * Outputs: pos       - 3D position vector in TEME frame in km
 *          vel       - 3D velocity vector in TEME frame in km/sec
 * Returns: 0         - Success
 *         -1         - Invalid inputs or parametres
 *          1         - Mean motion zero or less
 *          2         - Nonsensical orbital eccentricity (e >= 1; e < -0.00001)
 *          3         - Long periodics result error
 *          4         - Short period preliminary quantities error
 *          5         - Decaying satellite
 */
int
orbit_prop
(
    orbit* sat,
    time_t* time,
    unsigned int msec,
    unsigned int maxiter,
    double tolerance,
    vect* pos,
    vect* vel
)
{
  if ((sat == NULL) ||
      (time == NULL) ||
      (msec >= 1000) ||
      (maxiter < 1) ||
      (tolerance <= 0.0) ||
      (pos == NULL) ||
      (vel == NULL))
  {
    return -1;
  }

  double tdelta = difftime(*time + msec / 1000,
                           sat->epoch + sat->epoch_ms / 1000) / 60.0L;

  if (sat->isdeepspace == true)
  {
    return orbit_sdp4(sat, tdelta, maxiter, tolerance, pos, vel);
  }
  else
  {
    return orbit_sgp4(sat, tdelta, maxiter, tolerance, pos, vel);
  }
}
