/*
 * libtle2tcp.c
 *
 *  Created on: 28 ���. 2017 �.
 *      Author: Orson
 */

#include "libtle2tcp.h"
#include <math.h>
#include <stdio.h>

// ************************************************************************* //
//                             PRIVATE FUNCTIONS                             //
// ************************************************************************* //



// ************************************************************************* //
//                                 INTERFACE                                 //
// ************************************************************************* //

/*
 * Initialize SGP4 orbit from a TLE representation
 */
int
orbit_init(orbit* sat)
{
  // Initializing SGP4 propagation model

  // Convert to SGP4 units
  sat->no          = sat->no / rpd2radmin;
  sat->a           = pow(sat->no * tumin, (-2.0L / 3.0L));
  sat->nprimediv2  = sat->nprimediv2  / (rpd2radmin * 1440);
  sat->ndprimediv6 = sat->ndprimediv6 / (rpd2radmin * 1440 * 1440);

  // Standard orbital elements
  sat->i     = sat->i * deg2rad;
  sat->alpha = sat->alpha * deg2rad;
  sat->omega = sat->omega * deg2rad;
  sat->Mo    = sat->Mo * deg2rad;
  sat->alta  = sat->a * (1.0 + sat->e) - 1.0;
  sat->altp  = sat->a * (1.0 - sat->e) - 1.0;

  // Aux epoch quantities
  double esq    = pow(sat->e, 2);
  double omeosq = 1.0 - esq;
  double rteosq = sqrt(omeosq);
  double cosi   = cos(sat->i);
  double cosi2  = pow(cosi, 2);

  // Un-Kozai the mean motion
  double ak    = pow(xke / sat->no, 2.0L / 3.0L);
  double d1    = 0.75 * j2 * (3.0 * cosi2 - 1.0) / (rteosq * omeosq);
  double del   = d1 / pow(ak, 2);
  double adel  = ak * (1.0 - del * del - del * (1.0L / 3.0L + 134.0 * del
                 * del / 81.0));
  del          = d1 / (adel * adel);

  sat->no      = sat->no / (1.0 + del);

  double ao    = pow(xke / sat->no, (2.0L / 3.0L));
  double sini  = sin(sat->i);

  double po    = ao * omeosq;
  double con42 = 1.0 - 5.0 * cosi2;

  sat->con41   = -con42 -2 * cosi2;

  double ainv  = 1.0 / ao;
  double posq  = pow(po, 2);
  double rp    = ao * (1.0 - sat->e);

  // Sidereal time at epoch
  sat->gsto = 0; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ((omeosq >= 0.0 ) || (sat->no >= 0.0))
  {
    sat->isimp    = false;
    if (rp < (220.0 / Re + 1.0))
      sat->isimp  = 1;

    double s4  = 78.0 / Re + 1.0;
    double qzms24 = pow((120.0 - 78.0) / Re, 4);
    double perige = (rp - 1.0) * Re;

    // For perigees below 156km, recompute s and qoms2t
    if (perige < 156.0)
    {
      s4 = perige - 78.0;
      if (perige < 98.0)
        s4   = 20.0;

      qzms24 = pow((120.0 - s4) / Re, 4);
      s4     = s4 / Re + 1.0;
    }

    double pinvsq = 1.0 / posq;

    double tsi    = 1.0 / (ao - s4);
    sat->eta      = ao * sat->e * tsi;

    double etasq  = pow(sat->eta, 2);
    double eeta   = sat->e * sat->eta;
    double psisq  = fabs(1.0 - etasq);
    double coef   = qzms24 * pow(tsi, 4.0);
    double coef1  = coef / pow(psisq, 3.5);
    double C2     = coef1 * sat->no * (ao * (1.0 + 1.5 * etasq + eeta *
                    (4.0 + etasq)) + 0.375 * j2 * tsi / psisq *
                    sat->con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));

    sat->C1       = sat->Bstar * C2;

    double C3 = 0.0;
    if (sat->e > 1.0e-4)
      C3 = -2.0 * coef * tsi * j3divj2 * sat->no * sini / sat->e;

    sat->x1mth2   = 1.0 - cosi2;
    sat->C4       = 2.0* sat->no * coef1 * ao * omeosq * (sat->eta * (2.0 + 0.5
                    * etasq) + sat->e * (0.5 + 2.0 * etasq) - j2 * tsi / (ao *
                    psisq) * (-3.0 * sat->con41 * (1.0 - 2.0 * eeta + etasq *
                    (1.5 - 0.5 * eeta)) + 0.75 * sat->x1mth2 * (2.0 * etasq -
                    eeta * (1.0 + etasq)) * cos(2.0 * sat->omega)));

    sat->C5       = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) +
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

    double xhdot1     = -temp1 * cosi;

    sat->nodedot= xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosi2) + 2.0 *
                        temp3 * (3.0 - 7.0 * cosi2)) * cosi;

    double xpidot     =  sat->omegaprime + sat->nodedot;
    double sinomega   = sin(sat->omega);
    double cosomega   = cos(sat->omega);

    sat->omgcof  = sat->Bstar * C3 * cosomega;
    sat->xmcof   = 0.0;

    if (sat->e > 1.0e-4)
      sat->xmcof = -(2.0L / 3.0L) * coef * sat->Bstar / eeta;

    sat->nodecf  = 3.5 * omeosq * xhdot1 * sat->C1;
    sat->t2cof   = 1.5 * sat->C1;

    // Division by zero protection in case xinco = 180deg
    if (fabs(cosi + 1.0) > 1.5e-12)
      sat->xlcof = -0.25 * j3divj2 * sini * (3.0 + 5.0 * cosi) / (1.0 + cosi);
    else
      sat->xlcof = -0.25 * j3divj2 * sini * (3.0 + 5.0 * cosi) / 1.5e-12;

    sat->aycof   = -0.5 * j3divj2 * sini;

    sat->delMo   = pow(1.0 + sat->eta * cos(sat->Mo), 3);
    sat->sinMo   = sin(sat->Mo);
    sat->x7thm1  = 7.0 * cosi2 - 1.0;

    printf("NEW: no=%f\t sinMo=%f\t x7thm1=%f\n", sat->no, sat->sinMo, sat->x7thm1);

    // Deep space initializations
    if ((twopi / sat->no) >= 225.0)
    {
      sat->isdeepspace = true;
      sat->isimp       = 1;

      double tc        =  0.0;
      double im        = sat->i;

      /*dscom
          (
            epoch, satrec->ecco, satrec->argpo, tc, satrec->inclo, satrec->nodeo,
            satrec->no, &snodm, &cnodm,  &sini, &cosim,&sinomm,     &cosomm,
            &day, &satrec->e3, &satrec->ee2, &em,         &emsq, &gam,
            &satrec->peo,  &satrec->pgho,   &satrec->pho, &satrec->pinco,
            &satrec->plo,  &rtemsq,        &satrec->se2, &satrec->se3,
            &satrec->sgh2, &satrec->sgh3,   &satrec->sgh4,
            &satrec->sh2,  &satrec->sh3,    &satrec->si2, &satrec->si3,
            &satrec->sl2,  &satrec->sl3,    &satrec->sl4, &s1, &s2, &s3, &s4, &s5,
            &s6,   &s7,   &ss1,  &ss2,  &ss3,  &ss4,  &ss5,  &ss6,  &ss7, &sz1, &sz2, &sz3,
            &sz11, &sz12, &sz13, &sz21, &sz22, &sz23, &sz31, &sz32, &sz33,
            &satrec->xgh2, &satrec->xgh3,   &satrec->xgh4, &satrec->xh2,
            &satrec->xh3,  &satrec->xi2,    &satrec->xi3,  &satrec->xl2,
            &satrec->xl3,  &satrec->xl4,    &nm, &z1, &z2, &z3, &z11,
            &z12, &z13, &z21, &z22, &z23, &z31, &z32, &z33,
            &satrec->zmol, &satrec->zmos
          );*/

      // Constants
      double zes     =  0.01675;
      double zel     =  0.05490;
      double c1ss    =  2.9864797e-6;
      double c1l     =  4.7968065e-7;
      double zsinis  =  0.39785416;
      double zcosis  =  0.91744867;
      double zcosgs  =  0.1945905;
      double zsings  = -0.98088458;

      /* --------------------- local variables ------------------------ */
      /*int lsflg;
      double a1    , a2    , a3    , a4    , a5    , a6    , a7    ,
         a8    , a9    , a10   , betasq, cc    , ctem  , stem  ,
         x1    , x2    , x3    , x4    , x5    , x6    , x7    ,
         x8    , xnodce, xnoi  , zcosg , zcosgl, zcosh , zcoshl,
         zcosi , zcosil, zsing , zsingl, zsinh , zsinhl, zsini ,
         zsinil, zx    , zy;*/

      double sinalpha  = sin(sat->alpha);
      double cisalpha  = cos(sat->alpha);
      double betasq    = 1.0 - esq;
      double rtemsq    = sqrt(betasq);

      // Initialize lunar solar terms
      double day    = sat->epoch + 18261.5 + tc / 1440.0;
      double xnodce = fmod(4.5236020 - 9.2422029e-4 * day, twopi);
      double stem   = sin(xnodce);
      double ctem   = cos(xnodce);
      double zcosil = 0.91375164 - 0.03568096 * ctem;
      double zsinil = sqrt(1.0 - zcosil * zcosil);
      double zsinhl = 0.089683511 * stem / zsinil;
      double zcoshl = sqrt(1.0 - zsinhl * zsinhl);
      double gam    = 5.8351514 + 0.0019443680 * day;
      double zx     = 0.39785416 * stem / zsinil;
      double zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
      zx            = atan2(zx, zy);
      zx            = gam + zx - xnodce;
      double zcosgl = cos(zx);
      double zsingl = sin(zx);

      // Apply solar terms
      zcosg = zcosgs;
      zsing = zsings;
      zcosi = zcosis;
      zsini = zsinis;
      zcosh = *cosalpha;
      zsinh = *sinalpha;
      cc    = c1ss;
      xnoi  = 1.0 / sat->no;

      for (lsflg = 1; lsflg <= 2; lsflg++)
        {
          a1  =   zcosg * zcosh + zsing * zcosi * zsinh;
          a3  =  -zsing * zcosh + zcosg * zcosi * zsinh;
          a7  =  -zcosg * zsinh + zsing * zcosi * zcosh;
          a8  =   zsing * zsini;
          a9  =   zsing * zsinh + zcosg * zcosi * zcosh;
          a10 =   zcosg * zsini;
          a2  =   *cosi * a7 + *sini * a8;
          a4  =   *cosi * a9 + *sini * a10;
          a5  =  -*sini * a7 + *cosi * a8;
          a6  =  -*sini * a9 + *cosi * a10;

          x1  =  a1 * *cosomega + a2 * *sinomega;
          x2  =  a3 * *cosomega + a4 * *sinomega;
          x3  = -a1 * *sinomega + a2 * *cosomega;
          x4  = -a3 * *sinomega + a4 * *cosomega;
          x5  =  a5 * *sinomega;
          x6  =  a6 * *sinomega;
          x7  =  a5 * *cosomega;
          x8  =  a6 * *cosomega;

          *z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
          *z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
          *z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
          *z1  =  3.0 *  (a1 * a1 + a2 * a2) + *z31 * *esq;
          *z2  =  6.0 *  (a1 * a3 + a2 * a4) + *z32 * *esq;
          *z3  =  3.0 *  (a3 * a3 + a4 * a4) + *z33 * *esq;
          *z11 = -6.0 * a1 * a5 + *esq *  (-24.0 * x1 * x7-6.0 * x3 * x5);
          *z12 = -6.0 *  (a1 * a6 + a3 * a5) + *esq *
                 (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
          *z13 = -6.0 * a3 * a6 + *esq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
          *z21 =  6.0 * a2 * a5 + *esq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
          *z22 =  6.0 *  (a4 * a5 + a2 * a6) + *esq *
                 (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
          *z23 =  6.0 * a4 * a6 + *esq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
          *z1  = *z1 + *z1 + betasq * *z31;
          *z2  = *z2 + *z2 + betasq * *z32;
          *z3  = *z3 + *z3 + betasq * *z33;
          *s3  = cc * xnoi;
          *s2  = -0.5 * *s3 / *rtesq;
          *s4  = *s3 * *rtesq;
          *s1  = -15.0 * sat->e * *s4;
          *s5  = x1 * x3 + x2 * x4;
          *s6  = x2 * x3 + x1 * x4;
          *s7  = x2 * x4 - x1 * x3;

          // Apply lunar terms
          if (lsflg == 1)
            {
              *ss1   = *s1;
              *ss2   = *s2;
              *ss3   = *s3;
              *ss4   = *s4;
              *ss5   = *s5;
              *ss6   = *s6;
              *ss7   = *s7;
              *sz1   = *z1;
              *sz2   = *z2;
              *sz3   = *z3;
              *sz11  = *z11;
              *sz12  = *z12;
              *sz13  = *z13;
              *sz21  = *z21;
              *sz22  = *z22;
              *sz23  = *z23;
              *sz31  = *z31;
              *sz32  = *z32;
              *sz33  = *z33;
              zcosg = zcosgl;
              zsing = zsingl;
              zcosi = zcosil;
              zsini = zsinil;
              zcosh = zcoshl * *cosalpha + zsinhl * *sinalpha;
              zsinh = *sinalpha * zcoshl - *cosalpha * zsinhl;
              cc    = c1l;
           }
        }

      *zmol = fmod(4.7199672 + 0.22997150  * *day - *gam, twopi);
      *zmos = fmod(6.2565837 + 0.017201977 * *day, twopi);

      // Apply solar terms
      *se2  =   2.0 * *ss1 * *ss6;
      *se3  =   2.0 * *ss1 * *ss7;
      *si2  =   2.0 * *ss2 * *sz12;
      *si3  =   2.0 * *ss2 * (*sz13 - *sz11);
      *sl2  =  -2.0 * *ss3 * *sz2;
      *sl3  =  -2.0 * *ss3 * (*sz3 - *sz1);
      *sl4  =  -2.0 * *ss3 * (-21.0 - 9.0 * *esq) * zes;
      *sgh2 =   2.0 * *ss4 * *sz32;
      *sgh3 =   2.0 * *ss4 * (*sz33 - *sz31);
      *sgh4 = -18.0 * *ss4 * zes;
      *sh2  =  -2.0 * *ss2 * *sz22;
      *sh3  =  -2.0 * *ss2 * (*sz23 - *sz21);

      // Apply lunar terms
      *ee2  =   2.0 * *s1 * *s6;
      *e3   =   2.0 * *s1 * *s7;
      *xi2  =   2.0 * *s2 * *z12;
      *xi3  =   2.0 * *s2 * (*z13 - *z11);
      *xl2  =  -2.0 * *s3 * *z2;
      *xl3  =  -2.0 * *s3 * (*z3 - *z1);
      *xl4  =  -2.0 * *s3 * (-21.0 - 9.0 * *esq) * zel;
      *xgh2 =   2.0 * *s4 * *z32;
      *xgh3 =   2.0 * *s4 * (*z33 - *z31);
      *xgh4 = -18.0 * *s4 * zel;
      *xh2  =  -2.0 * *s2 * *z22;
      *xh3  =  -2.0 * *s2 * (*z23 - *z21);

  }

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
