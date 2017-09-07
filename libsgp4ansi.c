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

#include <ctype.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "libsgp4ansi.h"
#include "epoch.h"
#include "const.h"
#include "coord.h"
#include "vector.h"

// ************************************************************************* //
//                                VERSION                                    //
// ************************************************************************* //

const char library_version[] = VERSION;

// ************************************************************************* //
//                            PRIVATE PROTOTYPES                             //
// ************************************************************************* //

// Trim leading and trailing whitespaces from a string
size_t
str_trim(const char*, size_t, char*);

// Evaluate cubic polynomial
inline double
__attribute__((always_inline))
eval_cubic_poly(const double, const double, const double, const double,
                const double);

// Initialize deep space terms
inline void
sat_deep_init(sat*, const double, const double, const double, const double,
         const double, const double, const double, const double, const double);

// Run deep space integrator
inline void
sat_deep_dot(sat*, struct intv*);

// ************************************************************************* //
//                             PRIVATE FUNCTIONS                             //
// ************************************************************************* //

/*
 * Printout sat struct
 *
 * Inputs:  s       - sat struct to print
 *          caption - Caption to add to output
 * Outputs: None
 * Returns: None
 */
void
sat_print(sat* s, const char* caption)
{
  printf("------> %s\n", caption);
  printf("// TLE\n");
  printf("tle.name:\t\t%s\n", s->tle.name);
  printf("tle.int_designator:\t%s\n", s->tle.int_designator);
  printf("tle.sec_class:\t\t%c\n", s->tle.sec_class);
  printf("tle.epoch:\t\t%s\n", ctime(&s->tle.epoch));
  printf("tle.epoch_ms:\t\t%f\n", s->tle.epoch_ms);
  printf("tle.epoch_jul:\t\t%lf\n", s->tle.epoch_jul);
  printf("tle.mean_motion_dt2:\t%lf\n", s->tle.mean_motion_dt2);
  printf("tle.mean_motion_ddt6:\t%lf\n", s->tle.mean_motion_ddt6);
  printf("tle.bstar:\t\t%lf\n", s->tle.Bstar);
  printf("tle.inclination:\t%lf\n", s->tle.inclination);
  printf("tle.right_asc_node:\t%lf\n", s->tle.right_asc_node);
  printf("tle.eccentricity:\t%lf\n", s->tle.eccentricity);
  printf("tle.argument_perigee:\t%lf\n", s->tle.argument_perigee);
  printf("tle.mean_anomaly:\t%lf\n", s->tle.mean_anomaly);
  printf("tle.mean_motion:\t%lf\n", s->tle.mean_motion);
  printf("tle.norad_number:\t%d\n", s->tle.norad_number);
  printf("tle.orbit_number:\t%d\n", s->tle.orbit_number);
  printf("// Common constants\n");
  printf("comm.cosio:\t\t%lf\n", s->comm.cosio);
  printf("comm.sinio:\t\t%lf\n", s->comm.sinio);
  printf("comm.eta:\t\t%lf\n", s->comm.eta);
  printf("comm.t2cof:\t\t%lf\n", s->comm.t2cof);
  printf("comm.a3ovk2:\t\t%lf\n", s->comm.a3ovk2);
  printf("comm.x1mth2:\t\t%lf\n", s->comm.x1mth2);
  printf("comm.x3thm1:\t\t%lf\n", s->comm.x3thm1);
  printf("comm.x7thm1:\t\t%lf\n", s->comm.x7thm1);
  printf("comm.aycof:\t\t%lf\n", s->comm.aycof);
  printf("comm.xlcof:\t\t%lf\n", s->comm.xlcof);
  printf("comm.xnodcf:\t\t%lf\n", s->comm.xnodcf);
  printf("comm.c1:\t\t%lf\n", s->comm.c1);
  printf("comm.c4:\t\t%lf\n", s->comm.c4);
  printf("comm.omgdot:\t\t%lf\n", s->comm.omgdot);
  printf("comm.xnodot:\t\t%lf\n", s->comm.xnodot);
  printf("comm.xmdot:\t\t%lf\n", s->comm.xmdot);
  printf("comm.xnodp:\t\t%lf\n", s->comm.xnodp);
  printf("comm.aodp:\t\t%lf\n", s->comm.aodp);
  printf("comm.perigee:\t\t%lf\n", s->comm.perigee);
  printf("comm.period:\t\t%lf\n", s->comm.period);
  printf("comm.use_simple_model:\t%d\n", s->comm.use_simple_model);
  printf("comm.is_deep_space:\t%d\n", s->comm.is_deep_space);
  printf("// Near space constants\n");
  printf("near.c5:\t\t%lf\n", s->near.c5);
  printf("near.omgcof:\t\t%lf\n", s->near.omgcof);
  printf("near.xmcof:\t\t%lf\n", s->near.xmcof);
  printf("near.delmo:\t\t%lf\n", s->near.delmo);
  printf("near.sinmo:\t\t%lf\n", s->near.sinmo);
  printf("near.d2:\t\t%lf\n", s->near.d2);
  printf("near.d3:\t\t%lf\n", s->near.d3);
  printf("near.d4:\t\t%lf\n", s->near.d4);
  printf("near.t3cof:\t\t%lf\n", s->near.t3cof);
  printf("near.t4cof:\t\t%lf\n", s->near.t4cof);
  printf("near.t5cof:\t\t%lf\n", s->near.t5cof);
  printf("// Deep space constants\n");
  printf("deep.gsto:\t\t%lf\n", s->deep.gsto);
  printf("deep.zmol:\t\t%lf\n", s->deep.zmol);
  printf("deep.zmos:\t\t%lf\n", s->deep.zmos);
  printf("// Lunar and Solar terms for epoch\n");
  printf("deep.sse:\t\t%lf\n", s->deep.sse);
  printf("deep.ssi:\t\t%lf\n", s->deep.ssi);
  printf("deep.ssl:\t\t%lf\n", s->deep.ssl);
  printf("deep.ssg:\t\t%lf\n", s->deep.ssg);
  printf("deep.ssh:\t\t%lf\n", s->deep.ssh);
  printf("// Lunar and Solar secular terms\n");
  printf("deep.se2:\t\t%lf\n", s->deep.se2);
  printf("deep.si2:\t\t%lf\n", s->deep.si2);
  printf("deep.sl2:\t\t%lf\n", s->deep.sl2);
  printf("deep.sgh2:\t\t%lf\n", s->deep.sgh2);
  printf("deep.sh2:\t\t%lf\n", s->deep.sh2);
  printf("deep.se3:\t\t%lf\n", s->deep.se3);
  printf("deep.si3:\t\t%lf\n", s->deep.si3);
  printf("deep.sl3:\t\t%lf\n", s->deep.sl3);
  printf("deep.sgh3:\t\t%lf\n", s->deep.sgh3);
  printf("deep.sh3:\t\t%lf\n", s->deep.sh3);
  printf("deep.sl4:\t\t%lf\n", s->deep.sl4);
  printf("deep.sgh4:\t\t%lf\n", s->deep.sgh4);
  printf("deep.ee2:\t\t%lf\n", s->deep.ee2);
  printf("deep.e3:\t\t%lf\n", s->deep.e3);
  printf("deep.xi2:\t\t%lf\n", s->deep.xi2);
  printf("deep.xi3:\t\t%lf\n", s->deep.xi3);
  printf("deep.xl2:\t\t%lf\n", s->deep.xl2);
  printf("deep.xl3:\t\t%lf\n", s->deep.xl3);
  printf("deep.xl4:\t\t%lf\n", s->deep.xl4);
  printf("deep.xgh2:\t\t%lf\n", s->deep.xgh2);
  printf("deep.xgh3:\t\t%lf\n", s->deep.xgh3);
  printf("deep.xgh4:\t\t%lf\n", s->deep.xgh4);
  printf("deep.xh2:\t\t%lf\n", s->deep.xh2);
  printf("deep.xh3:\t\t%lf\n", s->deep.xh3);
  printf("// Lunar and Solar dot terms\n");
  printf("deep.d2201:\t\t%lf\n", s->deep.d2201);
  printf("deep.d2211:\t\t%lf\n", s->deep.d2211);
  printf("deep.d3210:\t\t%lf\n", s->deep.d3210);
  printf("deep.d3222:\t\t%lf\n", s->deep.d3222);
  printf("deep.d4410:\t\t%lf\n", s->deep.d4410);
  printf("deep.d4422:\t\t%lf\n", s->deep.d4422);
  printf("deep.d5220:\t\t%lf\n", s->deep.d5220);
  printf("deep.d5232:\t\t%lf\n", s->deep.d5232);
  printf("deep.d5421:\t\t%lf\n", s->deep.d5421);
  printf("deep.d5433:\t\t%lf\n", s->deep.d5433);
  printf("deep.del1:\t\t%lf\n", s->deep.del1);
  printf("deep.del2:\t\t%lf\n", s->deep.del2);
  printf("deep.del3:\t\t%lf\n", s->deep.del3);
  printf("// Geopotential resonance (12h)\n");
  printf("deep.is_resonant:\t%d\n", s->deep.is_resonant);
  printf("// Geosynchronous resonance (24h)\n");
  printf("deep.is_synchronous:\t%d\n", s->deep.is_synchronous);
  printf("// Integrator values for epoch\n");
  printf("intv0.xndot:\t\t%lf\n", s->intv0.xndot);
  printf("intv0.xnddt:\t\t%lf\n", s->intv0.xnddt);
  printf("intv0.xldot:\t\t%lf\n", s->intv0.xldot);
  printf("// Integrator values for current a_time\n");
  printf("intvt.xndot:\t\t%lf\n", s->intvt.xndot);
  printf("intvt.xnddt:\t\t%lf\n", s->intvt.xnddt);
  printf("intvt.xldot:\t\t%lf\n", s->intvt.xldot);
  printf("// Integrator constants\n");
  printf("intc.xfact:\t\t%lf\n", s->intc.xfact);
  printf("intc.xlamo:\t\t%lf\n", s->intc.xlamo);
  printf("// Integrator parametres\n");
  printf("intp.xli:\t\t%lf\n", s->intp.xli);
  printf("intp.xni:\t\t%lf\n", s->intp.xni);
  printf("intp.atime:\t\t%lf\n", s->intp.atime);
}
/*
 * Trim leading and trailing whitespaces from a string
 *
 * Inputs:  len   - Length of the input in interest
 *          str   - String to trim
 * Outputs: out   - Trimmed and null terminated string, truncated to size
 * Returns: Any   - Length of resulting string
 *         -1     - Failure to read TLE lines
 */
size_t str_trim(const char* str, size_t len, char* out)
{
  if(len == 0)
  {
    return 0;
  }

  const char *end;
  size_t out_size;

  // Trim leading space
  while(isspace((unsigned char)* str)) str++;

  if(*str == 0)  // Whitespace only
  {
    *out = 0;
    return 1;
  }

  // Trim trailing spaces
  end = str + strlen(str) - 1;
  while(end > str && isspace((unsigned char)* end)) end--;
  end++;

  // Set output size to minimum of trimmed string length and buffer size minus 1
  out_size = ((end - str) < len - 1) ? (end - str) : len - 1;

  // Copy trimmed string to out
  memcpy(out, str, out_size);
  out[out_size] = 0; // Terminate

  return out_size;
}

/*
 * Evaluate cubic polynomial
 *
 * Inputs:  x  - Argument
 *          c  - Constant term
 *          a1 - 1st power coefficient
 *          a2 - 2nd power coefficient
 *          a3 - 3rd power coefficient
 * Outputs: None
 * Returns: Evaluated polynomial
 */
double
eval_cubic_poly
(
  const double x,
  const double c,
  const double a1,
  const double a2,
  const double a3
)
{
  return c + x * (a1 + x * (a2 + x * a3));
}

/*
 * Initialize deep space terms
 *
 * Inputs:  many - selected terms from near space initialization
 * Outputs: s    - sat struct with deep space terms initialized
 * Returns: None
 */
void
sat_deep_init
(
    sat* s,
    const double eosq,
    const double sinio,
    const double cosio,
    const double betao,
    const double theta2,
    const double betao2,
    const double xmdot,
    const double omgdot,
    const double xnodot
)
//at_deep_init(eosq, s->comm.sinio, s->comm.cosio, betao, theta2, betao2,
//s->comm.xmdot, s->comm.omgdot, s->comm.xnodot);
{
  double se    = 0.0;
  double si    = 0.0;
  double sl    = 0.0;
  double sgh   = 0.0;
  double shdq  = 0.0;
  double bfact = 0.0;

  static const double ZNS    =  1.19459E-5;
  static const double C1SS   =  2.9864797E-6;
  static const double ZES    =  0.01675;
  static const double ZNL    =  1.5835218E-4;
  static const double C1L    =  4.7968065E-7;
  static const double ZEL    =  0.05490;
  static const double ZCOSIS =  0.91744867;
  static const double ZSINI  =  0.39785416;
  static const double ZSINGS = -0.98088458;
  static const double ZCOSGS =  0.1945905;
  static const double Q22    =  1.7891679E-6;
  static const double Q31    =  2.1460748E-6;
  static const double Q33    =  2.2123015E-7;
  static const double ROOT22 =  1.7891679E-6;
  static const double ROOT32 =  3.7393792E-7;
  static const double ROOT44 =  7.3636953E-9;
  static const double ROOT52 =  1.1428639E-7;
  static const double ROOT54 =  2.1765803E-9;

  const double aqnv   = 1.0 / s->comm.aodp;
  const double xpidot = omgdot + xnodot;
  const double sinq   = sin(s->tle.right_asc_node);
  const double cosq   = cos(s->tle.right_asc_node);
  const double sing   = sin(s->tle.argument_perigee);
  const double cosg   = cos(s->tle.argument_perigee);

  // Initialuze lunar and solar terms
  const double jday   = s->tle.epoch_jul - JULIAN_JAN1_12H_2000;

  const double xnodce = 4.5236020 - 9.2422029e-4 * jday;
  const double xnodce_temp = fmod(xnodce, TWOPI);
  const double stem   = sin(xnodce_temp);
  const double ctem   = cos(xnodce_temp);
  const double zcosil = 0.91375164 - 0.03568096 * ctem;
  const double zsinil = sqrt(1.0 - zcosil * zcosil);
  const double zsinhl = 0.089683511 * stem / zsinil;
  const double zcoshl = sqrt(1.0 - zsinhl * zsinhl);
  const double c      = 4.7199672 + 0.22997150 * jday;
  const double gam    = 5.8351514 + 0.0019443680 * jday;

  s->deep.zmol        = c - gam;
  // Wrap around to [0, 2pi]
  s->deep.zmol       -= TWOPI * floor((c - gam) / TWOPI);

  double zx           = 0.39785416 * stem / zsinil;
  double zy           = zcoshl * ctem + 0.91744867 * zsinhl * stem;
  zx                  = atan2(zx, zy);
  zx                  = fmod(gam + zx - xnodce, TWOPI);
  const double zcosgl = cos(zx);
  const double zsingl = sin(zx);

  // Wrap around to [0, 2pi]
  s->deep.zmos        = 6.2565837 + 0.017201977 * jday;
  s->deep.zmos       -= TWOPI * floor((6.2565837 + 0.017201977 * jday)
                                      / TWOPI);

  // Calculate solar terms
  double zcosg = ZCOSGS;
  double zsing = ZSINGS;
  double zcosi = ZCOSIS;
  double zsini = ZSINI;
  double zcosh = cosq;
  double zsinh = sinq;
  double cc = C1SS;
  double zn = ZNS;
  double ze = ZES;
  const double xnoi = 1.0 / s->comm.xnodp;

  for (int cnt = 0; cnt < 2; cnt++)
  {
      // Repeat solar term calculation a second time after linart terms
      const double a1  = zcosg * zcosh + zsing * zcosi * zsinh;
      const double a3  = -zsing * zcosh + zcosg * zcosi * zsinh;
      const double a7  = -zcosg * zsinh + zsing * zcosi * zcosh;
      const double a8  = zsing * zsini;
      const double a9  = zsing * zsinh + zcosg * zcosi*zcosh;
      const double a10 = zcosg * zsini;
      const double a2  = cosio * a7 + sinio * a8;
      const double a4  = cosio * a9 + sinio * a10;
      const double a5  = -sinio * a7 + cosio * a8;
      const double a6  = -sinio * a9 + cosio * a10;
      const double x1  = a1 * cosg + a2 * sing;
      const double x2  = a3 * cosg + a4 * sing;
      const double x3  = -a1 * sing + a2 * cosg;
      const double x4  = -a3 * sing + a4 * cosg;
      const double x5  = a5 * sing;
      const double x6  = a6 * sing;
      const double x7  = a5 * cosg;
      const double x8  = a6 * cosg;
      const double z31 = 12.0 * x1 * x1 - 3. * x3 * x3;
      const double z32 = 24.0 * x1 * x2 - 6. * x3 * x4;
      const double z33 = 12.0 * x2 * x2 - 3. * x4 * x4;
      double z1        = 3.0 * (a1 * a1 + a2 * a2) + z31 * eosq;
      double z2        = 6.0 * (a1 * a3 + a2 * a4) + z32 * eosq;
      double z3        = 3.0 * (a3 * a3 + a4 * a4) + z33 * eosq;

      const double z11 = -6.0 * a1 * a5 + eosq
                       * (-24. * x1 * x7 - 6. * x3 * x5);
      const double z12 = -6.0 * (a1 * a6 + a3 * a5) + eosq
                       *(-24. * (x2 * x7 + x1 * x8) - 6. * (x3 * x6 + x4 * x5));
      const double z13 = -6.0 * a3 * a6 + eosq
                       * (-24. * x2 * x8 - 6. * x4 * x6);
      const double z21 = 6.0 * a2 * a5 + eosq
                       * (24. * x1 * x5 - 6. * x3 * x7);
      const double z22 = 6.0 * (a4 * a5 + a2 * a6) + eosq
                       * (24. * (x2 * x5 + x1 * x6) - 6. * (x4 * x7 + x3 * x8));
      const double z23 = 6.0 * a4 * a6 + eosq
                       * (24. * x2 * x6 - 6. * x4 * x8);

      z1 = z1 + z1 + betao2 * z31;
      z2 = z2 + z2 + betao2 * z32;
      z3 = z3 + z3 + betao2 * z33;

      const double s3 = cc * xnoi;
      const double s2 = -0.5 * s3 / betao;
      const double s4 = s3 * betao;
      const double s1 = -15.0 * s->tle.eccentricity * s4;
      const double s5 = x1 * x3 + x2 * x4;
      const double s6 = x2 * x3 + x1 * x4;
      const double s7 = x2 * x4 - x1 * x3;

      se  = s1 * zn * s5;
      si  = s2 * zn * (z11 + z13);
      sl  = -zn * s3 * (z1 + z3 - 14.0 - 6.0 * eosq);
      sgh = s4 * zn * (z31 + z33 - 6.0);

      // Fix for certain inclination
      if (s->tle.inclination < 5.2359877e-2 ||
        s->tle.inclination > PI - 5.2359877e-2)
      {
        shdq = 0.0;
      }
      else
      {
        shdq = (-zn * s2 * (z21 + z23)) / sinio;
      }

      s->deep.ee2  =  2.0  * s1 * s6;
      s->deep.e3   =  2.0  * s1 * s7;
      s->deep.xi2  =  2.0  * s2 * z12;
      s->deep.xi3  =  2.0  * s2 * (z13 - z11);
      s->deep.xl2  = -2.0  * s3 * z2;
      s->deep.xl3  = -2.0  * s3 * (z3 - z1);
      s->deep.xl4  = -2.0  * s3 * (-21.0 - 9.0 * eosq) * ze;
      s->deep.xgh2 =  2.0  * s4 * z32;
      s->deep.xgh3 =  2.0  * s4 * (z33 - z31);
      s->deep.xgh4 = -18.0 * s4 * ze;
      s->deep.xh2  = -2.0  * s2 * z22;
      s->deep.xh3  = -2.0  * s2 * (z23 - z21);

      if (cnt == 1)
      {
          break;
      }
      // Calculate lunar terms
      s->deep.sse  = se;
      s->deep.ssi  = si;
      s->deep.ssl  = sl;
      s->deep.ssh  = shdq;
      s->deep.ssg  = sgh - cosio * s->deep.ssh;
      s->deep.se2  = s->deep.ee2;
      s->deep.si2  = s->deep.xi2;
      s->deep.sl2  = s->deep.xl2;
      s->deep.sgh2 = s->deep.xgh2;
      s->deep.sh2  = s->deep.xh2;
      s->deep.se3  = s->deep.e3;
      s->deep.si3  = s->deep.xi3;
      s->deep.sl3  = s->deep.xl3;
      s->deep.sgh3 = s->deep.xgh3;
      s->deep.sh3  = s->deep.xh3;
      s->deep.sl4  = s->deep.xl4;
      s->deep.sgh4 = s->deep.xgh4;

      zcosg = zcosgl;
      zsing = zsingl;
      zcosi = zcosil;
      zsini = zsinil;
      zcosh = zcoshl * cosq + zsinhl * sinq;
      zsinh = sinq * zcoshl - cosq * zsinhl;
      zn    = ZNL;
      cc    = C1L;
      ze    = ZEL;
  }

  s->deep.sse += se;
  s->deep.ssi += si;
  s->deep.ssl += sl;
  s->deep.ssg += sgh - cosio * shdq;
  s->deep.ssh += shdq;

  s->deep.is_resonant     = false;
  s->deep.is_synchronous   = false;
  bool do_init_integrator    = true;

  if (s->comm.xnodp < 0.0052359877
   && s->comm.xnodp > 0.0034906585)
  {
    /*
     * 24h synchronous resonance terms initialisation
     */
    s->deep.is_resonant = true;
    s->deep.is_synchronous = true;

    const double g200 = 1.0 + eosq * (-2.5 + 0.8125 * eosq);
    const double g310 = 1.0 + 2.0 * eosq;
    const double g300 = 1.0 + eosq * (-6.0 + 6.60937 * eosq);
    const double f220 = 0.75 * (1.0 + cosio) * (1.0 + cosio);
    const double f311 = 0.9375 * sinio * sinio * (1.0 + 3.0 * cosio) - 0.75
                      * (1.0 + cosio);
    double f330       = 1.0 + cosio;
    f330              = 1.875 * f330 * f330 * f330;
    s->deep.del1      = 3.0 * s->comm.xnodp * s->comm.xnodp * aqnv * aqnv;
    s->deep.del2      = 2.0 * s->deep.del1 * f220 * g200 * Q22;
    s->deep.del3      = 3.0 * s->deep.del1 * f330 * g300 * Q33 * aqnv;
    s->deep.del1      = s->deep.del1 * f311 * g310 * Q31 * aqnv;
    s->intc.xlamo     = s->tle.mean_anomaly + s->tle.right_asc_node
                      + s->tle.argument_perigee - s->deep.gsto;
    bfact             = xmdot + xpidot - THDT;
    bfact            += s->deep.ssl + s->deep.ssg + s->deep.ssh;
  }
  else if (s->comm.xnodp       < 8.26e-3
        || s->comm.xnodp       > 9.24e-3
        || s->tle.eccentricity < 0.5)
  {
    do_init_integrator = false;
  }
  else
  {
    // Geopotential resonant terms initialization
    s->deep.is_resonant = true;

    double g211;
    double g310;
    double g322;
    double g410;
    double g422;
    double g520;
    double g201 = -0.306 - (s->tle.eccentricity - 0.64) * 0.440;

    if (s->tle.eccentricity <= 0.65)
    {
      g211 = eval_cubic_poly(s->tle.eccentricity,
                             3.616, -13.247, +16.290, 0.0);
      g310 = eval_cubic_poly(s->tle.eccentricity,
                             -19.302, 117.390, -228.419, 156.591);
      g322 = eval_cubic_poly(s->tle.eccentricity,
                             -18.9068, 109.7927, -214.6334, 146.5816);
      g410 = eval_cubic_poly(s->tle.eccentricity,
                             -41.122, 242.694, -471.094, 313.953);
      g422 = eval_cubic_poly(s->tle.eccentricity,
                             -146.407, 841.880, -1629.014, 1083.435);
      g520 = eval_cubic_poly(s->tle.eccentricity,
                             -532.114, 3017.977, -5740.032, 3708.276);
    }
    else
    {
      g211 = eval_cubic_poly(s->tle.eccentricity,
                             -72.099, 331.819, -508.738, 266.724);
      g310 = eval_cubic_poly(s->tle.eccentricity,
                             -346.844, 1582.851, -2415.925, 1246.113);
      g322 = eval_cubic_poly(s->tle.eccentricity,
                             -342.585, 1554.908, -2366.899, 1215.972);
      g410 = eval_cubic_poly(s->tle.eccentricity,
                             -1052.797, 4758.686, -7193.992, 3651.957);
      g422 = eval_cubic_poly(s->tle.eccentricity,
                             -3581.69, 16178.11, -24462.77, 12422.52);

      if (s->tle.eccentricity <= 0.715)
      {
        g520 = eval_cubic_poly(s->tle.eccentricity,
                               1464.74, -4664.75, 3763.64, 0.0);
      }
      else
      {
        g520 = eval_cubic_poly(s->tle.eccentricity,
                               -5149.66, 29936.92, -54087.36, 31324.56);
      }
    }

    double g533;
    double g521;
    double g532;

    if (s->tle.eccentricity < 0.7)
    {
      g533 = eval_cubic_poly(s->tle.eccentricity,
                             -919.2277, 4988.61, -9064.77, 5542.21);
      g521 = eval_cubic_poly(s->tle.eccentricity,
                             -822.71072, 4568.6173, -8491.4146, 5337.524);
      g532 = eval_cubic_poly(s->tle.eccentricity,
                             -853.666, 4690.25, -8624.77, 5341.4);
    }
    else
    {
      g533 = eval_cubic_poly(s->tle.eccentricity,
                             -37995.78, 161616.52, -229838.2, 109377.94);
      g521 = eval_cubic_poly(s->tle.eccentricity,
                             -51752.104, 218913.95, -309468.16, 146349.42);
      g532 = eval_cubic_poly(s->tle.eccentricity,
                             -40023.88, 170470.89, -242699.48, 115605.82);
    }

    const double sini2 = sinio * sinio;
    const double f220  = 0.75 * (1.0 + 2.0 * cosio + theta2);
    const double f221  = 1.5 * sini2;
    const double f321  = 1.875 * sinio * (1.0 - 2.0 * cosio - 3.0 * theta2);
    const double f322  = -1.875 * sinio * (1.0 + 2.0 * cosio - 3.0 * theta2);
    const double f441  = 35.0 * sini2 * f220;
    const double f442  = 39.3750 * sini2 * sini2;
    const double f522  = 9.84375 * sinio * (sini2 * (1.0 - 2.0 * cosio - 5.0
                       * theta2) + 0.33333333 * (-2.0 + 4.0 * cosio + 6.0
                       * theta2));
    const double f523  = sinio * (4.92187512 * sini2 * (-2.0 - 4.0 * cosio
                       + 10.0 * theta2) + 6.56250012 * (1.0 + 2.0 * cosio
                       - 3.0 * theta2));
    const double f542  = 29.53125 * sinio * (2.0 - 8.0 * cosio + theta2
                       * (-12.0 + 8.0 * cosio + 10.0 * theta2));
    const double f543  = 29.53125 * sinio * (-2.0 - 8.0 * cosio + theta2
                       * (12.0 + 8.0 * cosio - 10.0 * theta2));

    const double xno2 = s->comm.xnodp * s->comm.xnodp;
    const double ainv2 = aqnv * aqnv;

    double temp1  = 3.0   * xno2  * ainv2;
    double temp   = temp1 * ROOT22;
    s->deep.d2201 = temp  * f220  * g201;
    s->deep.d2211 = temp  * f221  * g211;
    temp1         = temp1 * aqnv;
    temp          = temp1 * ROOT32;
    s->deep.d3210 = temp  * f321  * g310;
    s->deep.d3222 = temp  * f322  * g322;
    temp1         = temp1 * aqnv;
    temp          = 2.0   * temp1 * ROOT44;
    s->deep.d4410 = temp  * f441  * g410;
    s->deep.d4422 = temp  * f442  * g422;
    temp1         = temp1 * aqnv;
    temp          = temp1 * ROOT52;
    s->deep.d5220 = temp  * f522  * g520;
    s->deep.d5232 = temp  * f523  * g532;
    temp          = 2.0   * temp1 * ROOT54;
    s->deep.d5421 = temp  * f542  * g521;
    s->deep.d5433 = temp  * f543  * g533;

    s->intc.xlamo = s->tle.mean_anomaly
                  + 2 * s->tle.right_asc_node
                  - 2 * s->deep.gsto;
    bfact         = xmdot + 2 * xnodot - 2 * THDT;
    bfact         = bfact + s->deep.ssl + 2 * s->deep.ssh;
  }

  if (do_init_integrator)
  {
    // Initialize integrator
    s->intc.xfact = bfact - s->comm.xnodp;
    s->intp.atime = 0.0;
    s->intp.xni   = s->comm.xnodp;
    s->intp.xli   = s->intc.xlamo;

    // Prepare dot terms at epoch
    sat_deep_dot(s, &s->intv0);
  }
}

/*
 * Run deep space integrator
 *
 * Inputs:  s     - sat struct with deep space terms initialized
 * Outputs: intv0 - intv struct to hold calculated terms
 * Returns: None
 */
void
sat_deep_dot(sat* s, struct intv* intv0)
{
  static const double G22   = 5.7686396;
  static const double G32   = 0.95240898;
  static const double G44   = 1.8014998;
  static const double G52   = 1.0508330;
  static const double G54   = 4.4108898;
  static const double FASX2 = 0.13130908;
  static const double FASX4 = 2.8843198;
  static const double FASX6 = 0.37448087;

  if (s->deep.is_synchronous)
  {

    intv0->xndot = s->deep.del1 * sin(s->intp.xli - FASX2)
                  + s->deep.del2 * sin(2.0 * (s->intp.xli - FASX4))
                  + s->deep.del3 * sin(3.0 * (s->intp.xli - FASX6));
    intv0->xnddt = s->deep.del1 * cos(s->intp.xli - FASX2)
                  + 2.0 * s->deep.del2 * cos(2.0 * (s->intp.xli - FASX4))
                  + 3.0 * s->deep.del3 * cos(3.0 * (s->intp.xli - FASX6));
  }
  else
  {
    const double xomi  = s->tle.argument_perigee + s->comm.omgdot
                       * s->intp.atime;
    const double x2omi = xomi + xomi;
    const double x2li  = s->intp.xli + s->intp.xli;

    intv0->xndot = s->deep.d2201 * sin(x2omi + s->intp.xli - G22)
                  + s->deep.d2211 * sin(s->intp.xli - G22)
                  + s->deep.d3210 * sin(xomi + s->intp.xli - G32)
                  + s->deep.d3222 * sin(-xomi + s->intp.xli - G32)
                  + s->deep.d4410 * sin(x2omi + x2li - G44)
                  + s->deep.d4422 * sin(x2li - G44)
                  + s->deep.d5220 * sin(xomi + s->intp.xli - G52)
                  + s->deep.d5232 * sin(-xomi + s->intp.xli - G52)
                  + s->deep.d5421 * sin(xomi + x2li - G54)
                  + s->deep.d5433 * sin(-xomi + x2li - G54);
    intv0->xnddt = s->deep.d2201 * cos(x2omi + s->intp.xli - G22)
                  + s->deep.d2211 * cos(s->intp.xli - G22)
                  + s->deep.d3210 * cos(xomi + s->intp.xli - G32)
                  + s->deep.d3222 * cos(-xomi + s->intp.xli - G32)
                  + s->deep.d5220 * cos(xomi + s->intp.xli - G52)
                  + s->deep.d5232 * cos(-xomi + s->intp.xli - G52)
                  + 2.0 * (s->deep.d4410 * cos(x2omi + x2li - G44)
                  + s->deep.d4422 * cos(x2li - G44)
                  + s->deep.d5421 * cos(xomi + x2li - G54)
                  + s->deep.d5433 * cos(-xomi + x2li - G54));
  }

  intv0->xldot  = s->intp.xni + s->intc.xfact;
  intv0->xnddt *= intv0->xldot;
}

// ************************************************************************* //
//                                    API                                    //
// ************************************************************************* //

/*
 * Initialize SGP4/SDP4 orbit model from a raw NORAD TLE lines
 *
 * Inputs:  name_str - String containing satellite name
 *          line1    - TLE line 1
 *          line2    - TLE line 2
 * Outputs: s        - sat struct pointer with full orbital data expanded
 * Returns: 0        - Success
 *         -1        - Corrup TLE
 *
 * Calls: orbit_init
 *
 * TODO: Add checksum check
 */
int
sat_load_tle(char* name_str, char* line1, char* line2, sat* s)
{
  // Pre-formatting the raw TLE input
  for (char j = 10; j <= 15; j++)
    if (line1[j] == ' ')
      line1[j] = '_';

  if (line1[44] != ' ')
    line1[43] = line1[44];

  line1[44] = '.';
  if (line1[7] == ' ')
    line1[7] = 'U';

  if (line1[9] == ' ')
    line1[9] = '.';

  for (char j = 45; j <= 49; j++)
    if (line1[j] == ' ')
      line1[j] = '0';

  if (line1[51] == ' ')
    line1[51] = '0';

  if (line1[53] != ' ')
    line1[52] = line1[53];

  line1[53] = '.';
  line2[25] = '.';

  for (char j = 26; j <= 32; j++)
    if (line2[j] == ' ')
      line2[j] = '0';

  if (line1[62] == ' ')
    line1[62] = '0';

  if (line1[68] == ' ')
    line1[68] = '0';

  str_trim(name_str, 24, s->tle.name);

  struct tm epoch_tm, prop_tm;

  unsigned int cardnum;
  unsigned int norad_number1;
  unsigned int norad_number2;
  char int_designator[9];
  int epochyr;
  double epochdays;
  double mean_motion_dt2;
  double mean_motion_ddt6;
  int nexp;
  double bstar;
  int bexp;
  int ephem_type;
  int elset_number;
  double inclination;
  double right_asc_node;
  double argument_perigee;
  double mean_anomaly;
  double mean_motion;

  int retval = sscanf(line1,"%2d %5ld %1c %9s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %3d ",
         &cardnum, &norad_number1, &s->tle.sec_class, &s->tle.int_designator, &epochyr,
         &epochdays,&mean_motion_dt2, &mean_motion_ddt6, &nexp, &bstar,
         &bexp, &ephem_type, &elset_number);

  if (retval != 13 || cardnum != 1)
  {
    return -1;
  }

  if (line2[52] == ' ') // check for minus sign
    retval = sscanf(line2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %10lf %6ld \n",
           &cardnum, &norad_number2, &inclination, &right_asc_node,
           &s->tle.eccentricity, &argument_perigee, &mean_anomaly,
           &mean_motion, &s->tle.orbit_number);
  else
    retval = sscanf(line2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld \n",
           &cardnum, &norad_number2, &inclination, &right_asc_node,
           &s->tle.eccentricity, &argument_perigee, &mean_anomaly,
           &mean_motion, &s->tle.orbit_number);

  if (retval != 9 || cardnum != 2 || norad_number1 != norad_number2)
  {
    return -1;
  }

  s->tle.norad_number = norad_number2;

  // Convert to SGP4 units
  if (fractday2unix(epochyr, epochdays, &s->tle.epoch, &s->tle.epoch_ms) != 0)
  {
    return -1;
  }
  s->tle.epoch_jul        = unix2jul(&s->tle.epoch, s->tle.epoch_ms);
  s->tle.mean_motion_dt2  = mean_motion_dt2  / (RPD2RADPM * 1440.0);
  s->tle.mean_motion_ddt6 = mean_motion_ddt6 * pow(10, nexp)
                          / (RPD2RADPM * 1440.0 * 1440);
  s->tle.Bstar            = bstar * pow(10, bexp);
  s->tle.inclination      = inclination * DEG2RAD;
  s->tle.right_asc_node   = right_asc_node * DEG2RAD;
  s->tle.argument_perigee = argument_perigee * DEG2RAD;
  s->tle.mean_anomaly     = mean_anomaly * DEG2RAD;
  s->tle.mean_motion      = mean_motion * TWOPI / 1440;

  // Recover original mean motion and semimajor axis
  const double a1     = pow(XKE / s->tle.mean_motion, TWOTHIRD);
  const double cosio  = cos(s->tle.inclination);
  const double theta2 = cosio * cosio;
  const double x3thm1 = 3.0 * theta2 - 1.0;
  const double eosq   = s->tle.eccentricity * s->tle.eccentricity;
  const double betao2 = 1.0 - eosq;
  const double betao  = sqrt(betao2);
  const double temp   = (1.5 * J2 / 2.0) * x3thm1 / (betao * betao2);
  const double del1   = temp / (a1 * a1);
  const double a0     = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + del1
                      * 134.0 / 81.0)));
  const double del0   = temp / (a0 * a0);

  s->comm.xnodp       = s->tle.mean_motion / (1.0 + del0);
  s->comm.aodp        = a0 / (1.0 - del0);

  // Find perigee and period
  s->comm.perigee     = (s->comm.aodp * (1.0 - s->tle.eccentricity) - 1) * RE;
  s->comm.period      = TWOPI / s->comm.xnodp;

  // Expand constants
  return sat_init(s);
}

/*
 * Expand SGP4/SDP4 orbit elements from an orbit containing NORAD TLE portion
 *
 * Inputs:  s  - sat struct pointer to an unexpanded elements set
 * Outputs: s  - sat struct pointer with full orbital data
 * Returns: 0  - Success
 *         -1  - Eccenricity out of range [0, 1)
 *         -2  - Inclination out of range [0, pi)
 *
 * TODO: add return values
 */
int
sat_init(sat* s)
{

  /*
   * error checks
   */
  if (s->tle.eccentricity < 0.0 || s->tle.eccentricity > 0.999999)
  {
    return -1;
  }

  if (s->tle.inclination < 0.0 || s->tle.inclination > PI)
  {
    return -2;
  }

  s->comm.cosio       = cos(s->tle.inclination);
  s->comm.sinio       = sin(s->tle.inclination);
  const double theta2 = s->comm.cosio * s->comm.cosio;
  s->comm.x3thm1      = 3.0 * theta2 - 1.0;
  const double eosq   = s->tle.eccentricity * s->tle.eccentricity;
  const double betao2 = 1.0 - eosq;
  const double betao  = sqrt(betao2);

  if (s->comm.period >= 225.0)
  {
    s->comm.is_deep_space    = true;
  }
  else
  {
    s->comm.is_deep_space    = false;
    s->comm.use_simple_model = false;

    // Use simple model for perigee less than 220 kilometers
    if (s->comm.perigee < 220.0)
    {
      s->comm.use_simple_model = true;
    }
  }

  // find new s4 and qoms2t for perigees below 156km
  double s4 = 1.0 + 78.0 / RE;
  double qoms24 = pow(((120.0 - 78.0) / RE), 4.0);
  if (s->comm.perigee < 156.0)
  {
    s4 = s->comm.perigee - 78.0;
    if (s->comm.perigee < 98.0)
    {
      s4 = 20.0;
    }
    qoms24 = pow((120.0 - s4) * 1.0 / RE, 4.0);
    s4 = s4 / RE + 1;
  }

   // Expand constants
   const double pinvsq = 1.0 / (s->comm.aodp * s->comm.aodp * betao2 * betao2);
   const double tsi    = 1.0 / (s->comm.aodp - s4);
   s->comm.eta         = s->comm.aodp * s->tle.eccentricity * tsi;
   const double etasq  = s->comm.eta * s->comm.eta;
   const double eeta   = s->tle.eccentricity * s->comm.eta;
   const double psisq  = fabs(1.0 - etasq);
   const double coef   = qoms24 * pow(tsi, 4.0);
   const double coef1  = coef / pow(psisq, 3.5);
   const double c2     = coef1 * s->comm.xnodp * (s->comm.aodp * (1.0 + 1.5
                       * etasq + eeta * (4.0 + etasq)) + 0.75 * J2DIV2 * tsi
                       / psisq * s->comm.x3thm1 * (8.0 + 3.0 * etasq
                       * (8.0 + etasq)));
   s->comm.c1          = s->tle.Bstar * c2;
   s->comm.a3ovk2      = -J3 / J2DIV2;
   s->comm.x1mth2      = 1.0 - theta2;
   s->comm.c4          = 2.0 * s->comm.xnodp * coef1 * s->comm.aodp * betao2
                       * (s->comm.eta * (2.0 + 0.5 * etasq)
                       + s->tle.eccentricity * (0.5 + 2.0 * etasq) - J2 * tsi
                       / (s->comm.aodp * psisq) * (-3.0 * s->comm.x3thm1 * (1.0
                       - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) + 0.75
                       * s->comm.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq))
                       * cos(2.0 * s->tle.argument_perigee)));
   const double theta4 = theta2 * theta2;
   const double temp1  = 3.0 * J2DIV2 * pinvsq * s->comm.xnodp;
   const double temp2  = temp1 * J2DIV2 * pinvsq;
   const double temp3  = 1.25 * J2DIV2 * pinvsq * pinvsq * s->comm.xnodp;
   s->comm.xmdot       = s->comm.xnodp + 0.5 * temp1 * betao * s->comm.x3thm1
                       + 0.0625 * temp2 * betao * (13.0 - 78.0 * theta2 + 137.0
                       * theta4);
   const double x1m5th = 1.0 - 5.0 * theta2;
   s->comm.omgdot      = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7.0 - 114.0
                       * theta2 + 395.0 * theta4) + temp3 * (3.0 - 36.0
                       * theta2 + 49.0 * theta4);
   const double xhdot1 = -temp1 * s->comm.cosio;
   s->comm.xnodot      = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0
                       * temp3 * (3.0 - 7.0 * theta2)) * s->comm.cosio;
   s->comm.xnodcf      = 3.5 * betao2 * xhdot1 * s->comm.c1;
   s->comm.t2cof       = 1.5 * s->comm.c1;

   if (fabs(s->comm.cosio + 1.0) > 1.5e-12)
   {
     s->comm.xlcof = 0.125 * s->comm.a3ovk2 * s->comm.sinio * (3.0 + 5.0 * s->comm.cosio) / (1.0 + s->comm.cosio);
   }
   else
   {
     s->comm.xlcof = 0.125 * s->comm.a3ovk2 * s->comm.sinio * (3.0 + 5.0 * s->comm.cosio) / 1.5e-12;
   }

   s->comm.aycof = 0.25 * s->comm.a3ovk2 * s->comm.sinio;
   s->comm.x7thm1 = 7.0 * theta2 - 1.0;

   if (s->comm.is_deep_space)
   {
     s->deep.gsto = jul2gst(s->tle.epoch_jul);

     // Initialize deep space constants
     sat_deep_init(s, eosq, s->comm.sinio, s->comm.cosio, betao, theta2, betao2,
                   s->comm.xmdot, s->comm.omgdot, s->comm.xnodot);
   }
   else
   {
     // Shortcut for round orbits
     double c3 = 0.0;
     if (s->tle.eccentricity > 1.0e-4)
     {
       c3 = coef * tsi * s->comm.a3ovk2 * s->comm.xnodp * s->comm.sinio
           / s->tle.eccentricity;
     }

     s->near.c5      = 2.0 * coef1 * s->comm.aodp * betao2 * (1.0 + 2.75
         * (etasq + eeta) + eeta * etasq);
     s->near.omgcof  = s->tle.Bstar * c3 * cos(s->tle.argument_perigee);

     // Shortcut for round orbits
     s->near.xmcof   = 0.0;
     if (s->tle.eccentricity > 1.0e-4)
     {
       s->near.xmcof = -TWOTHIRD * coef * s->tle.Bstar / eeta;
     }

     s->near.delmo = pow(1.0 + s->comm.eta * (cos(s->tle.mean_anomaly)), 3.0);
     s->near.sinmo = sin(s->tle.mean_anomaly);

     if (!s->comm.use_simple_model)
     {
       const double c1sq = s->comm.c1 * s->comm.c1;
       s->near.d2        = 4.0 * s->comm.aodp * tsi * c1sq;
       const double temp = s->near.d2 * tsi * s->comm.c1 / 3.0;
       s->near.d3        = (17.0 * s->comm.aodp + s4) * temp;
       s->near.d4        = 0.5 * temp * s->comm.aodp * tsi * (221.0
           * s->comm.aodp + 31.0 * s4) * s->comm.c1;
       s->near.t3cof     = s->near.d2 + 2.0 * c1sq;
       s->near.t4cof     = 0.25 * (3.0 * s->near.d3 + s->comm.c1
           * (12.0 * s->near.d2 + 10.0 * c1sq));
       s->near.t5cof     = 0.2 * (3.0 * s->near.d4 + 12.0 * s->comm.c1 *
           s->near.d3 + 6.0 * s->near.d2 * s->near.d2 + 15.0 *
           c1sq * (2.0 * s->near.d2 + c1sq));
     }
   }
  return 0;
}
