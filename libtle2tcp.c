/*
 * libtle2tcp.c
 *
 *  Created on: 28 рту. 2017 у.
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
orbit_init(sat* satellite)
{
  // Initializing SGP4 propagation model

  // Convert to SGP4 units
  satellite->no = satellite->no / rpd2radmin;
  satellite->a  = pow(satellite->no * tumin, (-2.0L / 3.0L));
  satellite->nprimediv2  = satellite->nprimediv2  / (rpd2radmin * 1440);
  satellite->ndprimediv6 = satellite->ndprimediv6 / (rpd2radmin * 1440 * 1440);

  // Standard orbital elements
  satellite->i     = satellite->i * deg2rad;
  satellite->alpha = satellite->alpha * deg2rad;
  satellite->omega = satellite->omega * deg2rad;
  satellite->Mo    = satellite->Mo * deg2rad;
  satellite->alta  = satellite->a * (1.0 + satellite->e) - 1.0;
  satellite->altp  = satellite->a * (1.0 - satellite->e) - 1.0;

  // Aux epoch quantities
  double eccsq  = satellite->e * satellite->e;
  double omeosq = 1.0 - eccsq;
  double rteosq = sqrt(omeosq);
  double cosio  = cos(satellite->i);
  double cosio2 = pow(cosio, 2);

  // Un-Kozai the mean motion
  double ak     = pow(xke / satellite->no, 2.0L / 3.0L);
  double d1     = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
  double del    = d1 / pow(ak, 2);
  double adel   = ak * (1.0 - del * del - del * (1.0L / 3.0L + 134.0 * del
                  * del / 81.0));
  del           = d1 / (adel * adel);

  satellite->no = satellite->no / (1.0 + del);

  double ao     = pow(xke / satellite->no, (2.0L / 3.0L));
  double sinio  = sin(satellite->i);

  double po     = ao * omeosq;
  double con42  = 1.0 - 5.0 * cosio2;

  // ???
  satellite->con41 = -con42 -2 * cosio2;

  double ainv  = 1.0 / ao;
  double posq  = pow(po, 2);
  double rp    = ao * (1.0 - satellite->e);

  // Sidereal time at epoch
  satellite->gsto = 0;

  if ((omeosq >= 0.0 ) || ( satellite->no >= 0.0))
  {
    satellite->isimp = false;
    if (rp < (220.0 / Re + 1.0))
      satellite->isimp = 1;

    double sfour  = 78.0 / Re + 1.0;
    double qzms24 = pow((120.0 - 78.0) / Re, 4);
    double perige = (rp - 1.0) * Re;

    // For perigees below 156km, recompute s and qoms2t
    if (perige < 156.0)
    {
      sfour = perige - 78.0;
      if (perige < 98.0)
        sfour = 20.0;

      qzms24 = pow((120.0 - sfour) / Re, 4);
      sfour = sfour / Re + 1.0;
    }

    double pinvsq  = 1.0 / posq;

    double tsi     = 1.0 / (ao - sfour);
    satellite->eta = ao * satellite->e * tsi;

    double etasq = satellite->eta * satellite->eta;
    double eeta  = satellite->e * satellite->eta;
    double psisq = fabs(1.0 - etasq);
    double coef  = qzms24 * pow(tsi, 4.0);
    double coef1 = coef / pow(psisq, 3.5);
    double C2   = coef1 * satellite->no * (ao * (1.0 + 1.5 * etasq + eeta *
                   (4.0 + etasq)) + 0.375 * j2 * tsi / psisq * satellite->con41
                   * (8.0 + 3.0 * etasq * (8.0 + etasq)));

    satellite->C1 = satellite->bstar * C2;

    double C3 = 0.0;
    if (satellite->e > 1.0e-4)
      C3 = -2.0 * coef * tsi * j3divj2 * satellite->no * sinio / satellite->e;

    satellite->x1mth2 = 1.0 - cosio2;
    satellite->C4     = 2.0* satellite->no * coef1 * ao * omeosq *
                        (satellite->eta * (2.0 + 0.5 * etasq) + satellite->e *
                        (0.5 + 2.0 * etasq) - j2 * tsi / (ao * psisq) *
                        (-3.0 * satellite->con41 * (1.0 - 2.0 * eeta + etasq *
                        (1.5 - 0.5 * eeta)) + 0.75 * satellite->x1mth2 *
                        (2.0 * etasq - eeta * (1.0 + etasq)) *
                        cos(2.0 * satellite->omega)));

    satellite->C5     = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
                        (etasq + eeta) + eeta * etasq);

    double cosio4     = pow(cosio2, 2);
    double temp1      = 1.5 * j2 * pinvsq * satellite->no;
    double temp2      = 0.5 * temp1 * j2 * pinvsq;
    double temp3      = -0.46875 * j4 * pinvsq * pinvsq * satellite->no;

    satellite->mdot   = satellite->no + 0.5 * temp1 * rteosq * satellite->con41
                        + 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 +
                        137.0 * cosio4);

    satellite->omegaprime = -0.5 * temp1 * con42 + 0.0625 * temp2 * (7.0 -
                        114.0 * cosio2 + 395.0 * cosio4) + temp3 * (3.0 - 36.0
                        * cosio2 + 49.0 * cosio4);

    double xhdot1     = -temp1 * cosio;

    satellite->nodedot= xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 *
                        temp3 * (3.0 - 7.0 * cosio2)) * cosio;

    double xpidot     =  satellite->omegaprime + satellite->nodedot;

    satellite->omgcof = satellite->bstar * C3 * cos(satellite->omega);
    satellite->xmcof  = 0.0;

    if (satellite->e > 1.0e-4)
      satellite->xmcof= -(2.0L / 3.0L) * coef * satellite->bstar / eeta;

    satellite->nodecf = 3.5 * omeosq * xhdot1 * satellite->C1;
    satellite->t2cof  = 1.5 * satellite->C1;

    // Division by zero protection in case xinco = 180deg
    double temp4      =   1.5e-12;
    if (fabs(cosio+1.0) > 1.5e-12)
      satellite->xlcof= -0.25 * j3divj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    else
      satellite->xlcof= -0.25 * j3divj2 * sinio * (3.0 + 5.0 * cosio) / temp4;

    satellite->aycof  = -0.5 * j3divj2 * sinio;

    satellite->delMo  = pow(1.0 + satellite->eta * cos(satellite->Mo), 3);
    satellite->sinMo = sin(satellite->Mo);
    satellite->x7thm1 = 7.0 * cosio2 - 1.0;

    printf("NEW: no=%f\t sinMo=%f\t x7thm1=%f\n", satellite->no, satellite->sinMo, satellite->x7thm1);
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
