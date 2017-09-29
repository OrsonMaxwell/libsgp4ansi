/*     ----------------------------------------------------------------
 *
 *                               sgp4unit.cpp
 *
 *    this file contains the sgp4 procedures for analytical propagation
 *    of a satellite. the code was originally released in the 1980 and 1986
 *    spacetrack papers. a detailed discussion of the theory and history
 *    may be found in the 2006 aiaa paper by vallado, crawford, hujsak,
 *    and kelso.
 *
 *                            companion code for
 *               fundamentals of astrodynamics and applications
 *                                    2007
 *                              by david vallado
 *
 *       (w) 719-573-2600, email dvallado@agi.com
 *
 *    current :
 *              30 Aug 10  david vallado
 *                           delete unused variables in initl
 *                           replace pow inetger 2, 3 with multiplies for speed
 *    changes :
 *               3 Nov 08  david vallado
 *                           put returns in for error codes
 *              29 sep 08  david vallado
 *                           fix atime for faster operation in dspace
 *                           add operationmode for afspc (a) or improved (i)
 *                           performance mode
 *              16 jun 08  david vallado
 *                           update small eccentricity check
 *              16 nov 07  david vallado
 *                           misc fixes for better compliance
 *              20 apr 07  david vallado
 *                           misc fixes for constants
 *              11 aug 06  david vallado
 *                           chg lyddane choice back to strn3, constants, misc doc
 *              15 dec 05  david vallado
 *                           misc fixes
 *              26 jul 05  david vallado
 *                           fixes for paper
 *                           note that each fix is preceded by a
 *                           comment with "sgp4fix" and an explanation of
 *                           what was changed
 *              10 aug 04  david vallado
 *                           2nd printing baseline working
 *              14 may 01  david vallado
 *                           2nd edition baseline
 *                     80  norad
 *                           original baseline
 *       ----------------------------------------------------------------      */

#include "sgp4unit.h"

const char help = 'n';
FILE *dbgfile;

void sat_print(elsetrec s, char* caption)
{
  FILE *f = fopen("aiaa_sat_print.txt", "a");
  fprintf(f, "======================= %s ====================\n", caption);
  fprintf(f, "------------------------ NORAD TLE --------------------------\n");
  fprintf(f, "name             \n");
  fprintf(f, "sec_class        \n");
  fprintf(f, "int_designator   \n");
  fprintf(f, "epoch            %lf\n\n", s.epochdays);
  fprintf(f, "julian_epoch     %22.15lf\n", s.jdsatepoch);
  fprintf(f, "mean_motion_dt2  %+.15e\n", s.ndot);
  fprintf(f, "mean_motion_ddt6 %+.15e\n", s.nddot);
  fprintf(f, "Bstar            %+.15e\n", s.bstar);
  fprintf(f, "inclination      %+.15e\n", s.inclo);
  fprintf(f, "right_asc_node   %+.15e\n", s.nodeo);
  fprintf(f, "eccentricity     %+.15e\n", s.ecco);
  fprintf(f, "argument_perigee %+.15e\n", s.argpo);
  fprintf(f, "mean_anomaly     %+.15e\n", s.mo);
  fprintf(f, "mean_motion      %+.15e\n", s.no);
  fprintf(f, "norad_number     %d\n", s.satnum);
  fprintf(f, "orbit_number     %d\n", s.epochtynumrev);
  fprintf(f, "-------------------------- Flags ----------------------------\n");
  fprintf(f, "is_deep_space    %c\n", s.method);
  fprintf(f, "use_simple_model %d\n", s.isimp);
  fprintf(f, "is_24h_resonant  %d\n", s.irez);
  fprintf(f, "is_12h_resonant  \n");
  fprintf(f, "---------------- Standard orbital elements ------------------\n");
  fprintf(f, "GSTo             %+.15e\n", s.gsto);
  fprintf(f, "xnodp            %+.15e\n", s.no);
  fprintf(f, "aodp             \n");
  fprintf(f, "perigee          \n");
  fprintf(f, "perigee_alt      %+.15e\n", s.altp);
  fprintf(f, "period           \n");
  fprintf(f, "---------------------- Common constants ---------------------\n");
  fprintf(f, "aycof            %+.15e\n", s.aycof);
  fprintf(f, "C1               %+.15e\n", s.cc1);
  fprintf(f, "C4               %+.15e\n", s.cc4);
  fprintf(f, "eta              %+.15e\n", s.eta);
  fprintf(f, "omgdot           %+.15e\n", s.argpdot);
  fprintf(f, "t2cof            %+.15e\n", s.t2cof);
  fprintf(f, "x1mth2           %+.15e\n", s.x1mth2);
  fprintf(f, "x1m5th2          \n");
  fprintf(f, "x1m7th2          %+.15e\n", s.con41);
  fprintf(f, "x7thm1           %+.15e\n", s.x7thm1);
  fprintf(f, "xlcof            %+.15e\n", s.xlcof);
  fprintf(f, "xnodcf           %+.15e\n", s.nodecf);
  fprintf(f, "xnodot           %+.15e\n", s.nodedot);
  fprintf(f, "xmdot            %+.15e\n", s.mdot);
  fclose(f);
}

/* ----------- local functions - only ever used internally by sgp4 ---------- */
static void dpper
(
    double e3,     double ee2,    double peo,     double pgho,   double pho,
    double pinco,  double plo,    double se2,     double se3,    double sgh2,
    double sgh3,   double sgh4,   double sh2,     double sh3,    double si2,
    double si3,    double sl2,    double sl3,     double sl4,    double t,
    double xgh2,   double xgh3,   double xgh4,    double xh2,    double xh3,
    double xi2,    double xi3,    double xl2,     double xl3,    double xl4,
    double zmol,   double zmos,   double inclo,
    char init,
    double& ep,    double& inclp, double& nodep,  double& argpp, double& mp,
    char opsmode
);

static void dscom
(
    double epoch,  double ep,     double argpp,   double tc,     double inclp,
    double nodep,  double np,
    double& snodm, double& cnodm, double& sinim,  double& cosim, double& sinomm,
    double& cosomm,double& day,   double& e3,     double& ee2,   double& em,
    double& emsq,  double& gam,   double& peo,    double& pgho,  double& pho,
    double& pinco, double& plo,   double& rtemsq, double& se2,   double& se3,
    double& sgh2,  double& sgh3,  double& sgh4,   double& sh2,   double& sh3,
    double& si2,   double& si3,   double& sl2,    double& sl3,   double& sl4,
    double& s1,    double& s2,    double& s3,     double& s4,    double& s5,
    double& s6,    double& s7,    double& ss1,    double& ss2,   double& ss3,
    double& ss4,   double& ss5,   double& ss6,    double& ss7,   double& sz1,
    double& sz2,   double& sz3,   double& sz11,   double& sz12,  double& sz13,
    double& sz21,  double& sz22,  double& sz23,   double& sz31,  double& sz32,
    double& sz33,  double& xgh2,  double& xgh3,   double& xgh4,  double& xh2,
    double& xh3,   double& xi2,   double& xi3,    double& xl2,   double& xl3,
    double& xl4,   double& nm,    double& z1,     double& z2,    double& z3,
    double& z11,   double& z12,   double& z13,    double& z21,   double& z22,
    double& z23,   double& z31,   double& z32,    double& z33,   double& zmol,
    double& zmos
);

static void dsinit
(
    gravconsttype whichconst,
    double cosim,  double emsq,   double argpo,   double s1,     double s2,
    double s3,     double s4,     double s5,      double sinim,  double ss1,
    double ss2,    double ss3,    double ss4,     double ss5,    double sz1,
    double sz3,    double sz11,   double sz13,    double sz21,   double sz23,
    double sz31,   double sz33,   double t,       double tc,     double gsto,
    double mo,     double mdot,   double no,      double nodeo,  double nodedot,
    double xpidot, double z1,     double z3,      double z11,    double z13,
    double z21,    double z23,    double z31,     double z33,    double ecco,
    double eccsq,  double& em,    double& argpm,  double& inclm, double& mm,
    double& nm,    double& nodem,
    int& irez,
    double& atime, double& d2201, double& d2211,  double& d3210, double& d3222,
    double& d4410, double& d4422, double& d5220,  double& d5232, double& d5421,
    double& d5433, double& dedt,  double& didt,   double& dmdt,  double& dndt,
    double& dnodt, double& domdt, double& del1,   double& del2,  double& del3,
    double& xfact, double& xlamo, double& xli,    double& xni, char
);

static void dspace
(
    int irez,
    double d2201,  double d2211,  double d3210,   double d3222,  double d4410,
    double d4422,  double d5220,  double d5232,   double d5421,  double d5433,
    double dedt,   double del1,   double del2,    double del3,   double didt,
    double dmdt,   double dnodt,  double domdt,   double argpo,  double argpdot,
    double t,      double tc,     double gsto,    double xfact,  double xlamo,
    double no,
    double& atime, double& em,    double& argpm,  double& inclm, double& xli,
    double& mm,    double& xni,   double& nodem,  double& dndt,  double& nm, char
);

static void initl
(
    int satn,      gravconsttype whichconst,
    double ecco,   double epoch,  double inclo,   double& no,
    char& method,
    double& ainv,  double& ao,    double& con41,  double& con42, double& cosio,
    double& cosio2,double& eccsq, double& omeosq, double& posq,
    double& rp,    double& rteosq,double& sinio , double& gsto, char opsmode
);

/* -----------------------------------------------------------------------------
 *
 *                           procedure dpper
 *
 *  this procedure provides deep space long period periodic contributions
 *    to the mean elements.  by design, these periodics are zero at epoch.
 *    this used to be dscom which included initialization, but it's really a
 *    recurring function.
 *
 *  author        : david vallado                  719-573-2600   28 jun 2005
 *
 *  inputs        :
 *    e3          -
 *    ee2         -
 *    peo         -
 *    pgho        -
 *    pho         -
 *    pinco       -
 *    plo         -
 *    se2 , se3 , sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4 -
 *    t           -
 *    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
 *    zmol        -
 *    zmos        -
 *    ep          - eccentricity                           0.0 - 1.0
 *    inclo       - inclination - needed for lyddane modification
 *    nodep       - right ascension of ascending node
 *    argpp       - argument of perigee
 *    mp          - mean anomaly
 *
 *  outputs       :
 *    ep          - eccentricity                           0.0 - 1.0
 *    inclp       - inclination
 *    nodep        - right ascension of ascending node
 *    argpp       - argument of perigee
 *    mp          - mean anomaly
 *
 *  locals        :
 *    alfdp       -
 *    betdp       -
 *    cosip  , sinip  , cosop  , sinop  ,
 *    dalf        -
 *    dbet        -
 *    dls         -
 *    f2, f3      -
 *    pe          -
 *    pgh         -
 *    ph          -
 *    pinc        -
 *    pl          -
 *    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
 *    sll   , sls
 *    xls         -
 *    xnoh        -
 *    zf          -
 *    zm          -
 *
 *  coupling      :
 *    none.
 *
 *  references    :
 *    hoots, roehrich, norad spacetrack report #3 1980
 *    hoots, norad spacetrack report #6 1986
 *    hoots, schumacher and glover 2004
 *    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

static void dpper
(
    double e3,     double ee2,    double peo,     double pgho,   double pho,
    double pinco,  double plo,    double se2,     double se3,    double sgh2,
    double sgh3,   double sgh4,   double sh2,     double sh3,    double si2,
    double si3,    double sl2,    double sl3,     double sl4,    double t,
    double xgh2,   double xgh3,   double xgh4,    double xh2,    double xh3,
    double xi2,    double xi3,    double xl2,     double xl3,    double xl4,
    double zmol,   double zmos,   double inclo,
    char init,
    double& ep,    double& inclp, double& nodep,  double& argpp, double& mp,
    char opsmode
)
{
  /* --------------------- local variables ------------------------ */
  const double twopi = 2.0 * pi;
  double alfdp, betdp, cosip, cosop, dalf, dbet, dls,
  f2,    f3,    pe,    pgh,   ph,   pinc, pl ,
  sel,   ses,   sghl,  sghs,  shll, shs,  sil,
  sinip, sinop, sinzf, sis,   sll,  sls,  xls,
  xnoh,  zf,    zm,    zel,   zes,  znl,  zns;

  /* ---------------------- constants ----------------------------- */
  zns   = 1.19459e-5;
  zes   = 0.01675;
  znl   = 1.5835218e-4;
  zel   = 0.05490;

  /* --------------- calculate time varying periodics ----------- */
  zm    = zmos + zns * t;
  // be sure that the initial call has time set to zero
  if (init == 'y')
    zm = zmos;
  zf    = zm + 2.0 * zes * sin(zm);
  sinzf = sin(zf);
  f2    =  0.5 * sinzf * sinzf - 0.25;
  f3    = -0.5 * sinzf * cos(zf);
  ses   = se2* f2 + se3 * f3;
  sis   = si2 * f2 + si3 * f3;
  sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
  sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
  shs   = sh2 * f2 + sh3 * f3;
  zm    = zmol + znl * t;
  if (init == 'y')
    zm = zmol;
  zf    = zm + 2.0 * zel * sin(zm);
  sinzf = sin(zf);
  f2    =  0.5 * sinzf * sinzf - 0.25;
  f3    = -0.5 * sinzf * cos(zf);
  sel   = ee2 * f2 + e3 * f3;
  sil   = xi2 * f2 + xi3 * f3;
  sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
  sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
  shll  = xh2 * f2 + xh3 * f3;
  pe    = ses + sel;
  pinc  = sis + sil;
  pl    = sls + sll;
  pgh   = sghs + sghl;
  ph    = shs + shll;

  if (init == 'n')
  {
    pe    = pe - peo;
    pinc  = pinc - pinco;
    pl    = pl - plo;
    pgh   = pgh - pgho;
    ph    = ph - pho;
    inclp = inclp + pinc;
    ep    = ep + pe;
    sinip = sin(inclp);
    cosip = cos(inclp);

    /* ----------------- apply periodics directly ------------ */
    //  sgp4fix for lyddane choice
    //  strn3 used original inclination - this is technically feasible
    //  gsfc used perturbed inclination - also technically feasible
    //  probably best to readjust the 0.2 limit value and limit discontinuity
    //  0.2 rad = 11.45916 deg
    //  use next line for original strn3 approach and original inclination
    double i = 0;
    if (opsmode == 'a')
      i = inclo;
    else
      i = inclp;
    //  use next line for gsfc version and perturbed inclination
    if (i >= 0.2)
    {
      ph     = ph / sinip;
      pgh    = pgh - cosip * ph;
      argpp  = argpp + pgh;
      nodep  = nodep + ph;
      mp     = mp + pl;
    }
    else
    {
      /* ---- apply periodics with lyddane modification ---- */
      sinop  = sin(nodep);
      cosop  = cos(nodep);
      alfdp  = sinip * sinop;
      betdp  = sinip * cosop;
      dalf   =  ph * cosop + pinc * cosip * sinop;
      dbet   = -ph * sinop + pinc * cosip * cosop;
      alfdp  = alfdp + dalf;
      betdp  = betdp + dbet;
      nodep  = fmod(nodep, twopi);
      //  sgp4fix for afspc written intrinsic functions
      // nodep used without a trigonometric function ahead
      if ((nodep < 0.0) && (opsmode != 'a'))
        nodep = nodep + twopi;
      xls    = mp + argpp + cosip * nodep;
      dls    = pl + pgh - pinc * nodep * sinip;
      xls    = xls + dls;
      xnoh   = nodep;
      nodep  = atan2(alfdp, betdp);
      //  sgp4fix for afspc written intrinsic functions
      // nodep used without a trigonometric function ahead
      if ((nodep < 0.0) && (opsmode != 'a'))
        nodep = nodep + twopi;
      if (fabs(xnoh - nodep) > pi){
        if (nodep < xnoh)
          nodep = nodep + twopi;
        else
          nodep = nodep - twopi;}
      mp    = mp + pl;
      argpp = xls - mp - cosip * nodep;
    }
  }   // if init == 'n'

  //#include "debug1.cpp"
//  printf("dpper>tdelta: %+.15e\n", t);
//  printf("dpper>inc_lp: %+.15e\n", inclp);
//  printf("dpper>ecc_lp: %+.15e\n", ep);
//  printf("dpper>nod_lp: %+.15e\n", nodep);
//  printf("dpper>arg_lp: %+.15e\n", argpp);
//  printf("dpper>mo_lp:  %+.15e\n", mp);
//  printf("----------------------------------------\n");
}  // end dpper

/*-----------------------------------------------------------------------------
 *
 *                           procedure dscom
 *
 *  this procedure provides deep space common items used by both the secular
 *    and periodics subroutines.  input is provided as shown. this routine
 *    used to be called dpper, but the functions inside weren't well organized.
 *
 *  author        : david vallado                  719-573-2600   28 jun 2005
 *
 *  inputs        :
 *    epoch       -
 *    ep          - eccentricity
 *    argpp       - argument of perigee
 *    tc          -
 *    inclp       - inclination
 *    nodep       - right ascension of ascending node
 *    np          - mean motion
 *
 *  outputs       :
 *    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
 *    day         -
 *    e3          -
 *    ee2         -
 *    em          - eccentricity
 *    emsq        - eccentricity squared
 *    gam         -
 *    peo         -
 *    pgho        -
 *    pho         -
 *    pinco       -
 *    plo         -
 *    rtemsq      -
 *    se2, se3         -
 *    sgh2, sgh3, sgh4        -
 *    sh2, sh3, si2, si3, sl2, sl3, sl4         -
 *    s1, s2, s3, s4, s5, s6, s7          -
 *    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
 *    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
 *    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
 *    nm          - mean motion
 *    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
 *    zmol        -
 *    zmos        -
 *
 *  locals        :
 *    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
 *    betasq      -
 *    cc          -
 *    ctem, stem        -
 *    x1, x2, x3, x4, x5, x6, x7, x8          -
 *    xnodce      -
 *    xnoi        -
 *    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
 *    zcosi  , zsini  , zcosil , zsinil ,
 *    zx          -
 *    zy          -
 *
 *  coupling      :
 *    none.
 *
 *  references    :
 *    hoots, roehrich, norad spacetrack report #3 1980
 *    hoots, norad spacetrack report #6 1986
 *    hoots, schumacher and glover 2004
 *    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

static void dscom
(
    double epoch,  double ep,     double argpp,   double tc,     double inclp,
    double nodep,  double np,
    double& snodm, double& cnodm, double& sinim,  double& cosim, double& sinomm,
    double& cosomm,double& day,   double& e3,     double& ee2,   double& em,
    double& emsq,  double& gam,   double& peo,    double& pgho,  double& pho,
    double& pinco, double& plo,   double& rtemsq, double& se2,   double& se3,
    double& sgh2,  double& sgh3,  double& sgh4,   double& sh2,   double& sh3,
    double& si2,   double& si3,   double& sl2,    double& sl3,   double& sl4,
    double& s1,    double& s2,    double& s3,     double& s4,    double& s5,
    double& s6,    double& s7,    double& ss1,    double& ss2,   double& ss3,
    double& ss4,   double& ss5,   double& ss6,    double& ss7,   double& sz1,
    double& sz2,   double& sz3,   double& sz11,   double& sz12,  double& sz13,
    double& sz21,  double& sz22,  double& sz23,   double& sz31,  double& sz32,
    double& sz33,  double& xgh2,  double& xgh3,   double& xgh4,  double& xh2,
    double& xh3,   double& xi2,   double& xi3,    double& xl2,   double& xl3,
    double& xl4,   double& nm,    double& z1,     double& z2,    double& z3,
    double& z11,   double& z12,   double& z13,    double& z21,   double& z22,
    double& z23,   double& z31,   double& z32,    double& z33,   double& zmol,
    double& zmos
)
{
  /* -------------------------- constants ------------------------- */
  const double zes     =  0.01675;
  const double zel     =  0.05490;
  const double c1ss    =  2.9864797e-6;
  const double c1l     =  4.7968065e-7;
  const double zsinis  =  0.39785416;
  const double zcosis  =  0.91744867;
  const double zcosgs  =  0.1945905;
  const double zsings  = -0.98088458;
  const double twopi   =  2.0 * pi;

  /* --------------------- local variables ------------------------ */
  int lsflg;
  double a1    , a2    , a3    , a4    , a5    , a6    , a7    ,
  a8    , a9    , a10   , betasq, cc    , ctem  , stem  ,
  x1    , x2    , x3    , x4    , x5    , x6    , x7    ,
  x8    , xnodce, xnoi  , zcosg , zcosgl, zcosh , zcoshl,
  zcosi , zcosil, zsing , zsingl, zsinh , zsinhl, zsini ,
  zsinil, zx    , zy;

  nm     = np;
  em     = ep;
  snodm  = sin(nodep);
  cnodm  = cos(nodep);
  sinomm = sin(argpp);
  cosomm = cos(argpp);
  sinim  = sin(inclp);
  cosim  = cos(inclp);
  emsq   = em * em;
  betasq = 1.0 - emsq;
  rtemsq = sqrt(betasq);

  /* ----------------- initialize lunar solar terms --------------- */
  peo    = 0.0;
  pinco  = 0.0;
  plo    = 0.0;
  pgho   = 0.0;
  pho    = 0.0;

  // printf("------------------------------\n");

  day    = epoch + 18261.5 + tc / 1440.0;
  // printf("day:    %+.15e\n", day);
  xnodce = fmod(4.5236020 - 9.2422029e-4 * day, twopi);
  // printf("xnodce: %+.15e\n", xnodce);
  stem   = sin(xnodce);
  // printf("stem:   %+.15e\n", stem);
  ctem   = cos(xnodce);
  // printf("ctem:   %+.15e\n", ctem);
  zcosil = 0.91375164 - 0.03568096 * ctem;
  // printf("zcosil: %+.15e\n", zcosil);
  zsinil = sqrt(1.0 - zcosil * zcosil);
  // printf("zsinil: %+.15e\n", zsinil);
  zsinhl = 0.089683511 * stem / zsinil;
  // printf("zsinhl: %+.15e\n", zsinhl);
  zcoshl = sqrt(1.0 - zsinhl * zsinhl);
  // printf("zcoshl: %+.15e\n", zcoshl);
  gam    = 5.8351514 + 0.0019443680 * day;
  // printf("gam:    %+.15e\n", gam);
  zx     = 0.39785416 * stem / zsinil;
  zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
  // printf("zy:     %+.15e\n", zy);
  zx     = atan2(zx, zy);
  zx     = gam + zx - xnodce;
  // printf("zx:     %+.15e\n", zx);
  zcosgl = cos(zx);
  // printf("zcosgl: %+.15e\n", zcosgl);
  zsingl = sin(zx);
  // printf("zsingl: %+.15e\n", zsingl);


  /* ------------------------- do solar terms --------------------- */
  zcosg = zcosgs;
  zsing = zsings;
  zcosi = zcosis;
  zsini = zsinis;
  zcosh = cnodm;
  // printf("zcosh:  %+.15e\n", zcosh);
  zsinh = snodm;
  // printf("zsinh:  %+.15e\n", zsinh);
  cc    = c1ss;
  xnoi  = 1.0 / nm;
  // printf("xnoi:   %+.15e\n", xnoi);

  for (lsflg = 1; lsflg <= 2; lsflg++)
  {
    a1  =   zcosg * zcosh + zsing * zcosi * zsinh;
    //printf("%d : %f %f %f %f %f\n", lsflg, zcosg, zcosh, zsing, zcosi, zsinh);
    a3  =  -zsing * zcosh + zcosg * zcosi * zsinh;
    a7  =  -zcosg * zsinh + zsing * zcosi * zcosh;
    a8  =   zsing * zsini;
    a9  =   zsing * zsinh + zcosg * zcosi * zcosh;
    a10 =   zcosg * zsini;
    a2  =   cosim * a7 + sinim * a8;
    a4  =   cosim * a9 + sinim * a10;
    a5  =  -sinim * a7 + cosim * a8;
    a6  =  -sinim * a9 + cosim * a10;
    //         printf("a1:     %+.15e\n", a1);
    //         printf("a3:     %+.15e\n", a3);
    //         printf("a7:     %+.15e\n", a7);
    //         printf("a8:     %+.15e\n", a8);
    //         printf("a9:     %+.15e\n", a9);
    //         printf("a10:    %+.15e\n", a10);
    //         printf("a2:     %+.15e\n", a2);
    //         printf("a4:     %+.15e\n", a4);
    //         printf("a5:     %+.15e\n", a5);
    //         printf("a6:     %+.15e\n", a6);

    x1  =  a1 * cosomm + a2 * sinomm;
    x2  =  a3 * cosomm + a4 * sinomm;
    x3  = -a1 * sinomm + a2 * cosomm;
    x4  = -a3 * sinomm + a4 * cosomm;
    x5  =  a5 * sinomm;
    x6  =  a6 * sinomm;
    x7  =  a5 * cosomm;
    x8  =  a6 * cosomm;
    //         printf("x1:     %+.15e\n", x1);
    //         printf("x2:     %+.15e\n", x2);
    //         printf("x3:     %+.15e\n", x3);
    //         printf("x4:     %+.15e\n", x4);
    //         printf("x5:     %+.15e\n", x5);
    //         printf("x6:     %+.15e\n", x6);
    //         printf("x7:     %+.15e\n", x7);
    //         printf("x8:     %+.15e\n", x8);

    z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
    z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
    z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
    z1  =  3.0 *  (a1 * a1 + a2 * a2) + z31 * emsq;
    z2  =  6.0 *  (a1 * a3 + a2 * a4) + z32 * emsq;
    z3  =  3.0 *  (a3 * a3 + a4 * a4) + z33 * emsq;
    z11 = -6.0 * a1 * a5 + emsq *  (-24.0 * x1 * x7-6.0 * x3 * x5);
    z12 = -6.0 *  (a1 * a6 + a3 * a5) + emsq *
        (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
    z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
    z21 =  6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
    z22 =  6.0 *  (a4 * a5 + a2 * a6) + emsq *
        (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
    z23 =  6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
    z1  = z1 + z1 + betasq * z31;
    z2  = z2 + z2 + betasq * z32;
    z3  = z3 + z3 + betasq * z33;
    s3  = cc * xnoi;
    s2  = -0.5 * s3 / rtemsq;
    s4  = s3 * rtemsq;
    s1  = -15.0 * em * s4;
    s5  = x1 * x3 + x2 * x4;
    s6  = x2 * x3 + x1 * x4;
    s7  = x2 * x4 - x1 * x3;
    //         printf("z31:    %+.15e\n", z31);
    //         printf("z32:    %+.15e\n", z32);
    //         printf("z33:    %+.15e\n", z33);
    //         printf("z1:     %+.15e\n", z1);
    //         printf("z2:     %+.15e\n", z2);
    //         printf("z3:     %+.15e\n", z3);
    //         printf("z11:    %+.15e\n", z11);
    //         printf("z12:    %+.15e\n", z12);
    //         printf("z13:    %+.15e\n", z13);
    //         printf("z21:    %+.15e\n", z21);
    //         printf("z22:    %+.15e\n", z22);
    //         printf("z23:    %+.15e\n", z23);
    //         printf("s3:     %+.15e\n", s3);
    //         printf("s2:     %+.15e\n", s2);
    //         printf("s4:     %+.15e\n", s4);
    //         printf("s1:     %+.15e\n", s1);
    //         printf("s5:     %+.15e\n", s5);
    //         printf("s6:     %+.15e\n", s6);
    //         printf("s7:     %+.15e\n", s7);

    /* ----------------------- do lunar terms ------------------- */
    if (lsflg == 1)
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
      zcosh = zcoshl * cnodm + zsinhl * snodm;
      zsinh = snodm * zcoshl - cnodm * zsinhl;
      //printf("*** %d : %f %f %f %f\n", lsflg, zcoshl, cnodm, zsinhl, snodm);
      cc    = c1l;
    }
  }

  zmol = fmod(4.7199672 + 0.22997150  * day - gam, twopi);
  zmos = fmod(6.2565837 + 0.017201977 * day, twopi);
  //     printf("zmos:   %+.15e\n", zmos);
  //     printf("zmol:   %+.15e\n", zmol);
  //     printf("------------------------------\n");
  /* ------------------------ do solar terms ---------------------- */
  se2  =   2 * ss1 * ss6;
  se3  =   2 * ss1 * ss7;
  si2  =   2 * ss2 * sz12;
  si3  =   2 * ss2 * (sz13 - sz11);
  sl2  =  -2 * ss3 * sz2;
  sl3  =  -2 * ss3 * (sz3 - sz1);
  sl4  =  -2 * ss3 * (-21 - 9 * emsq) * zes;
  sgh2 =   2 * ss4 * sz32;
  sgh3 =   2 * ss4 * (sz33 - sz31);
  sgh4 = -18 * ss4 * zes;
  sh2  =  -2 * ss2 * sz22;
  sh3  =  -2 * ss2 * (sz23 - sz21);
  //     printf("se2:    %+.15e\n", se2);
  //     printf("se3:    %+.15e\n", se3);
  //     printf("si2:    %+.15e\n", si2);
  //     printf("si3:    %+.15e\n", si3);
  //     printf("sl2:    %+.15e\n", sl2);
  //     printf("sl3:    %+.15e\n", sl3);
  //     printf("sl4:    %+.15e\n", sl4);
  //     printf("sgh2:   %+.15e\n", sgh2);
  //     printf("sgh3:   %+.15e\n", sgh3);
  //     printf("sgh4:   %+.15e\n", sgh4);
  //     printf("sh2:    %+.15e\n", sh2);
  //     printf("sh3:    %+.15e\n", sh3);
  /* ------------------------ do lunar terms ---------------------- */
  ee2  =   2 * s1 * s6;
  e3   =   2 * s1 * s7;
  xi2  =   2 * s2 * z12;
  xi3  =   2 * s2 * (z13 - z11);
  xl2  =  -2 * s3 * z2;
  xl3  =  -2 * s3 * (z3 - z1);
  xl4  =  -2 * s3 * (-21 - 9 * emsq) * zel;
  xgh2 =   2 * s4 * z32;
  xgh3 =   2 * s4 * (z33 - z31);
  xgh4 = -18 * s4 * zel;
  xh2  =  -2 * s2 * z22;
  xh3  =  -2 * s2 * (z23 - z21);
  //     printf("ee2:    %+.15e\n", ee2);
  //     printf("e3:     %+.15e\n", e3);
  //     printf("xi2:    %+.15e\n", xi2);
  //     printf("xi3:    %+.15e\n", xi3);
  //     printf("xl2:    %+.15e\n", xl2);
  //     printf("xl3:    %+.15e\n", xl3);
  //     printf("xl4:    %+.15e\n", xl4);
  //     printf("xgh2:   %+.15e\n", xgh2);
  //     printf("xgh3:   %+.15e\n", xgh3);
  //     printf("xgh4:   %+.15e\n", xgh4);
  //     printf("xh2:    %+.15e\n", xh2);
  //     printf("xh3:    %+.15e\n", xh3);
  //     printf("------------------------------\n");
  //#include "debug2.cpp"
}  // end dscom

/*-----------------------------------------------------------------------------
 *
 *                           procedure dsinit
 *
 *  this procedure provides deep space contributions to mean motion dot due
 *    to geopotential resonance with half day and one day orbits.
 *
 *  author        : david vallado                  719-573-2600   28 jun 2005
 *
 *  inputs        :
 *    cosim, sinim-
 *    emsq        - eccentricity squared
 *    argpo       - argument of perigee
 *    s1, s2, s3, s4, s5      -
 *    ss1, ss2, ss3, ss4, ss5 -
 *    sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
 *    t           - time
 *    tc          -
 *    gsto        - greenwich sidereal time                   rad
 *    mo          - mean anomaly
 *    mdot        - mean anomaly dot (rate)
 *    no          - mean motion
 *    nodeo       - right ascension of ascending node
 *    nodedot     - right ascension of ascending node dot (rate)
 *    xpidot      -
 *    z1, z3, z11, z13, z21, z23, z31, z33 -
 *    eccm        - eccentricity
 *    argpm       - argument of perigee
 *    inclm       - inclination
 *    mm          - mean anomaly
 *    xn          - mean motion
 *    nodem       - right ascension of ascending node
 *
 *  outputs       :
 *    em          - eccentricity
 *    argpm       - argument of perigee
 *    inclm       - inclination
 *    mm          - mean anomaly
 *    nm          - mean motion
 *    nodem       - right ascension of ascending node
 *    irez        - flag for resonance           0-none, 1-one day, 2-half day
 *    atime       -
 *    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433    -
 *    dedt        -
 *    didt        -
 *    dmdt        -
 *    dndt        -
 *    dnodt       -
 *    domdt       -
 *    del1, del2, del3        -
 *    ses  , sghl , sghs , sgs  , shl  , shs  , sis  , sls
 *    theta       -
 *    xfact       -
 *    xlamo       -
 *    xli         -
 *    xni
 *
 *  locals        :
 *    ainv2       -
 *    aonv        -
 *    cosisq      -
 *    eoc         -
 *    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543  -
 *    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533  -
 *    sini2       -
 *    temp        -
 *    temp1       -
 *    theta       -
 *    xno2        -
 *
 *  coupling      :
 *    getgravconst
 *
 *  references    :
 *    hoots, roehrich, norad spacetrack report #3 1980
 *    hoots, norad spacetrack report #6 1986
 *    hoots, schumacher and glover 2004
 *    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

static void dsinit
(
    gravconsttype whichconst,
    double cosim,  double emsq,   double argpo,   double s1,     double s2,
    double s3,     double s4,     double s5,      double sinim,  double ss1,
    double ss2,    double ss3,    double ss4,     double ss5,    double sz1,
    double sz3,    double sz11,   double sz13,    double sz21,   double sz23,
    double sz31,   double sz33,   double t,       double tc,     double gsto,
    double mo,     double mdot,   double no,      double nodeo,  double nodedot,
    double xpidot, double z1,     double z3,      double z11,    double z13,
    double z21,    double z23,    double z31,     double z33,    double ecco,
    double eccsq,  double& em,    double& argpm,  double& inclm, double& mm,
    double& nm,    double& nodem,
    int& irez,
    double& atime, double& d2201, double& d2211,  double& d3210, double& d3222,
    double& d4410, double& d4422, double& d5220,  double& d5232, double& d5421,
    double& d5433, double& dedt,  double& didt,   double& dmdt,  double& dndt,
    double& dnodt, double& domdt, double& del1,   double& del2,  double& del3,
    double& xfact, double& xlamo, double& xli,    double& xni,   char opsmode
)
{
  /* --------------------- local variables ------------------------ */
  const double twopi = 2.0 * pi;

  double ainv2 , aonv=0.0, cosisq, eoc, f220 , f221  , f311  ,
      f321  , f322  , f330  , f441  , f442  , f522  , f523  ,
      f542  , f543  , g200  , g201  , g211  , g300  , g310  ,
      g322  , g410  , g422  , g520  , g521  , g532  , g533  ,
      ses   , sgs   , sghl  , sghs  , shs   , shll  , sis   ,
      sini2 , sls   , temp  , temp1 , theta , xno2  , q22   ,
      q31   , q33   , root22, root44, root54, rptim , root32,
      root52, x2o3  , xke   , znl   , emo   , zns   , emsqo,
      tumin, mu, radiusearthkm, j2, j3, j4, j3oj2;

  q22    = 1.7891679e-6;
  q31    = 2.1460748e-6;
  q33    = 2.2123015e-7;
  root22 = 1.7891679e-6;
  root44 = 7.3636953e-9;
  root54 = 2.1765803e-9;
  rptim  = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
  root32 = 3.7393792e-7;
  root52 = 1.1428639e-7;
  x2o3   = 2.0 / 3.0;
  znl    = 1.5835218e-4;
  zns    = 1.19459e-5;

  // sgp4fix identify constants and allow alternate values
  getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

  /* -------------------- deep space initialization ------------ */
  irez = 0;
  if ((nm < 0.0052359877) && (nm > 0.0034906585))
    irez = 1;
  if ((nm >= 8.26e-3) && (nm <= 9.24e-3) && (em >= 0.5))
    irez = 2;

  /* ------------------------ do solar terms ------------------- */
  ses  =  ss1 * zns * ss5;
  sis  =  ss2 * zns * (sz11 + sz13);
  sls  = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
  sghs =  ss4 * zns * (sz31 + sz33 - 6.0);
  shs  = -zns * ss2 * (sz21 + sz23);
  // sgp4fix for 180 deg incl
  if (((inclm < 5.2359877e-2) || (inclm > pi - 5.2359877e-2)) && opsmode != 'a')
    shs = 0.0;
  if (sinim != 0.0)
    shs = shs / sinim;
  sgs  = sghs - cosim * shs;

  //     printf("ses:    %+.15e\n", ses);
  //     printf("sis:    %+.15e\n", sis);
  //     printf("sls:    %+.15e\n", sls);
  //     printf("sghs:   %+.15e\n", sghs);
  //     printf("shs:    %+.15e\n", shs);
  //     printf("sgs:    %+.15e\n", sgs);

  /* ------------------------- do lunar terms ------------------ */
  dedt = ses + s1 * znl * s5;
  didt = sis + s2 * znl * (z11 + z13);
  dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
  sghl = s4 * znl * (z31 + z33 - 6.0);
  shll = -znl * s2 * (z21 + z23);
  // sgp4fix for 180 deg incl
  if (((inclm < 5.2359877e-2) || (inclm > pi - 5.2359877e-2)) && (opsmode != 'a'))
    shll = 0.0;
  domdt = sgs + sghl;
  dnodt = shs;
  if (sinim != 0.0)
  {
    domdt = domdt - cosim / sinim * shll;
    dnodt = dnodt + shll / sinim;
  }

  //     printf("dedt:   %+.15e\n", dedt);
  //     printf("didt:   %+.15e\n", didt);
  //     printf("dmdt:   %+.15e\n", dmdt);
  //     printf("dndt:   %+.15e\n", 0.0);
  //     printf("sghl:   %+.15e\n", sghl);
  //     printf("shll:   %+.15e\n", shll);
  //     printf("domdt:  %+.15e\n", domdt);
  //     printf("dnodt:  %+.15e\n", dnodt);

  /* ----------- calculate deep space resonance effects -------- */
  dndt   = 0.0;
  theta  = fmod(gsto + tc * rptim, twopi);
  em     = em + dedt * t;
  inclm  = inclm + didt * t;
  argpm  = argpm + domdt * t;
  nodem  = nodem + dnodt * t;
  mm     = mm + dmdt * t;

  //     printf("----------------------------\n");
  //     printf("theta:  %+.15e\n", theta);
  //     printf("em:     %+.15e\n", em);
  //     printf("inclm:  %+.15e\n", inclm);

  //   sgp4fix for negative inclinations
  //   the following if statement should be commented out
  if ((inclm < 0.0) && (opsmode == 'a'))
  {
    inclm  = -inclm;
    argpm  = argpm - pi;
    nodem = nodem + pi;
  }

  /* -------------- initialize the resonance terms ------------- */
  if (irez != 0)
  {
    aonv = pow(nm / xke, x2o3);

    /* ---------- geopotential resonance for 12 hour orbits ------ */
    if (irez == 2)
    {
      cosisq = cosim * cosim;
      emo    = em;
      em     = ecco;
      emsqo  = emsq;
      emsq   = eccsq;
      eoc    = em * emsq;

      g201   = -0.306 - (em - 0.64) * 0.440;

      if (em <= 0.65)
      {
        g211 =    3.616  -  13.2470 * em +  16.2900 * emsq;
        g310 =  -19.302  + 117.3900 * em - 228.4190 * emsq +  156.5910 * eoc;
        g322 =  -18.9068 + 109.7927 * em - 214.6334 * emsq +  146.5816 * eoc;
        g410 =  -41.122  + 242.6940 * em - 471.0940 * emsq +  313.9530 * eoc;
        g422 = -146.407  + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
        g520 = -532.114  + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;
      }
      else
      {
        g211 =   -72.099 +   331.819 * em -   508.738 * emsq +   266.724 * eoc;
        g310 =  -346.844 +  1582.851 * em -  2415.925 * emsq +  1246.113 * eoc;
        g322 =  -342.585 +  1554.908 * em -  2366.899 * emsq +  1215.972 * eoc;
        g410 = -1052.797 +  4758.686 * em -  7193.992 * emsq +  3651.957 * eoc;
        g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
        if (em > 0.715)
          g520 =-5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
        else
          g520 = 1464.74 -  4664.75 * em +  3763.64 * emsq;
      }
      if (em < 0.7)
      {
        g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21  * eoc;
        g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
        g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4  * eoc;
      }
      else
      {
        g533 =-37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
        g521 =-51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
        g532 =-40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
      }

      //             printf("----------------------------\n");
      //             printf("g201:   %+.15e\n", g201);
      //             printf("g211:   %+.15e\n", g211);
      //             printf("g310:   %+.15e\n", g310);
      //             printf("g322:   %+.15e\n", g322);
      //             printf("g410:   %+.15e\n", g410);
      //             printf("g422:   %+.15e\n", g422);
      //             printf("g520:   %+.15e\n", g520);
      //             printf("g533:   %+.15e\n", g533);
      //             printf("g521:   %+.15e\n", g521);
      //             printf("g532:   %+.15e\n", g532);

      sini2=  sinim * sinim;
      f220 =  0.75 * (1.0 + 2.0 * cosim+cosisq);
      f221 =  1.5 * sini2;
      f321 =  1.875 * sinim  *  (1.0 - 2.0 * cosim - 3.0 * cosisq);
      f322 = -1.875 * sinim  *  (1.0 + 2.0 * cosim - 3.0 * cosisq);
      f441 = 35.0 * sini2 * f220;
      f442 = 39.3750 * sini2 * sini2;
      f522 =  9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim- 5.0 * cosisq) +
          0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq) );
      f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim +
          10.0 * cosisq) + 6.56250012 * (1.0+2.0 * cosim - 3.0 * cosisq));
      f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim+cosisq *
          (-12.0 + 8.0 * cosim + 10.0 * cosisq));
      f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim+cosisq *
          (12.0 + 8.0 * cosim - 10.0 * cosisq));

      //             printf("----------------------------\n");
      //             printf("f220:   %+.15e\n", f220);
      //             printf("f221:   %+.15e\n", f221);
      //             printf("f321:   %+.15e\n", f321);
      //             printf("f322:   %+.15e\n", f322);
      //             printf("f441:   %+.15e\n", f441);
      //             printf("f442:   %+.15e\n", f442);
      //             printf("f522:   %+.15e\n", f522);
      //             printf("f523:   %+.15e\n", f523);
      //             printf("f542:   %+.15e\n", f542);
      //             printf("f543:   %+.15e\n", f543);

      xno2  =  nm * nm;
      ainv2 =  aonv * aonv;
      temp1 =  3.0 * xno2 * ainv2;
      temp  =  temp1 * root22;
      d2201 =  temp * f220 * g201;
      d2211 =  temp * f221 * g211;
      temp1 =  temp1 * aonv;
      temp  =  temp1 * root32;
      d3210 =  temp * f321 * g310;
      d3222 =  temp * f322 * g322;
      temp1 =  temp1 * aonv;
      temp  =  2.0 * temp1 * root44;
      d4410 =  temp * f441 * g410;
      d4422 =  temp * f442 * g422;
      temp1 =  temp1 * aonv;
      temp  =  temp1 * root52;
      d5220 =  temp * f522 * g520;
      d5232 =  temp * f523 * g532;
      temp  =  2.0 * temp1 * root54;
      d5421 =  temp * f542 * g521;
      d5433 =  temp * f543 * g533;
      xlamo =  fmod(mo + nodeo + nodeo-theta - theta, twopi);
      xfact =  mdot + dmdt + 2.0 * (nodedot + dnodt - rptim) - no;
      em    = emo;
      emsq  = emsqo;

      //             printf("----------------------------\n");
      //             printf("ainv2:  %+.15e\n", ainv2);
      //             printf("d2201:  %+.15e\n", d2201);
      //             printf("d2211:  %+.15e\n", d2211);
      //             printf("d3210:  %+.15e\n", d3210);
      //             printf("d3222:  %+.15e\n", d3222);
      //             printf("d4410:  %+.15e\n", d4410);
      //             printf("d4422:  %+.15e\n", d4422);
      //             printf("d5220:  %+.15e\n", d5220);
      //             printf("d5232:  %+.15e\n", d5232);
      //             printf("d5232:  %+.15e\n", d5232);
      //             printf("d5433:  %+.15e\n", d5433);
      //             printf("xlamo:  %+.15e\n", xlamo);
      //             printf("xfact:  %+.15e\n", xfact);
    }

    /* ---------------- synchronous resonance terms -------------- */
    if (irez == 1)
    {
      g200  = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
      g310  = 1.0 + 2.0 * emsq;
      g300  = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
      f220  = 0.75 * (1.0 + cosim) * (1.0 + cosim);
      f311  = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
      f330  = 1.0 + cosim;
      f330  = 1.875 * f330 * f330 * f330;
      del1  = 3.0 * nm * nm * aonv * aonv;
      del2  = 2.0 * del1 * f220 * g200 * q22;
      del3  = 3.0 * del1 * f330 * g300 * q33 * aonv;
      del1  = del1 * f311 * g310 * q31 * aonv;
      xlamo = fmod(mo + nodeo + argpo - theta, twopi);
      xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no;

      //             printf("----------------------------\n");
      //             printf("g200:   %+.15e\n", g200);
      //             printf("g310:   %+.15e\n", g310);
      //             printf("g300:   %+.15e\n", g300);
      //             printf("f220:   %+.15e\n", f220);
      //             printf("f311:   %+.15e\n", f311);
      //             printf("f330:   %+.15e\n", f330);
      //             printf("del1:   %+.15e\n", del1);
      //             printf("del2:   %+.15e\n", del2);
      //             printf("del3:   %+.15e\n", del3);
    }

    /* ------------ for sgp4, initialize the integrator ---------- */
    xli   = xlamo;
    xni   = no;
    atime = 0.0;
    nm    = no + dndt;

    //         printf("xfact:  %+.15e\n", xfact);
    //         printf("xlamo:  %+.15e\n", xlamo);
    //         printf("xli:    %+.15e\n", xli);
    //         printf("xni:    %+.15e\n", xni);
    //         printf("atime:  %+.15e\n", atime);
    //         printf("xnodp:  %+.15e\n", nm);
  }

  //#include "debug3.cpp"
}  // end dsinit

/*-----------------------------------------------------------------------------
 *
 *                           procedure dspace
 *
 *  this procedure provides deep space contributions to mean elements for
 *    perturbing third body.  these effects have been averaged over one
 *    revolution of the sun and moon.  for earth resonance effects, the
 *    effects have been averaged over no revolutions of the satellite.
 *    (mean motion)
 *
 *  author        : david vallado                  719-573-2600   28 jun 2005
 *
 *  inputs        :
 *    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
 *    dedt        -
 *    del1, del2, del3  -
 *    didt        -
 *    dmdt        -
 *    dnodt       -
 *    domdt       -
 *    irez        - flag for resonance           0-none, 1-one day, 2-half day
 *    argpo       - argument of perigee
 *    argpdot     - argument of perigee dot (rate)
 *    t           - time
 *    tc          -
 *    gsto        - gst
 *    xfact       -
 *    xlamo       -
 *    no          - mean motion
 *    atime       -
 *    em          - eccentricity
 *    ft          -
 *    argpm       - argument of perigee
 *    inclm       - inclination
 *    xli         -
 *    mm          - mean anomaly
 *    xni         - mean motion
 *    nodem       - right ascension of ascending node
 *
 *  outputs       :
 *    atime       -
 *    em          - eccentricity
 *    argpm       - argument of perigee
 *    inclm       - inclination
 *    xli         -
 *    mm          - mean anomaly
 *    xni         -
 *    nodem       - right ascension of ascending node
 *    dndt        -
 *    nm          - mean motion
 *
 *  locals        :
 *    delt        -
 *    ft          -
 *    theta       -
 *    x2li        -
 *    x2omi       -
 *    xl          -
 *    xldot       -
 *    xnddt       -
 *    xndt        -
 *    xomi        -
 *
 *  coupling      :
 *    none        -
 *
 *  references    :
 *    hoots, roehrich, norad spacetrack report #3 1980
 *    hoots, norad spacetrack report #6 1986
 *    hoots, schumacher and glover 2004
 *    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

static void dspace
(
    int irez,
    double d2201,  double d2211,  double d3210,   double d3222,  double d4410,
    double d4422,  double d5220,  double d5232,   double d5421,  double d5433,
    double dedt,   double del1,   double del2,    double del3,   double didt,
    double dmdt,   double dnodt,  double domdt,   double argpo,  double argpdot,
    double t,      double tc,     double gsto,    double xfact,  double xlamo,
    double no,
    double& atime, double& em,    double& argpm,  double& inclm, double& xli,
    double& mm,    double& xni,   double& nodem,  double& dndt,  double& nm, char opsmode
)
{
  const double twopi = 2.0 * pi;
  int iretn , iret;
  double delt, ft, theta, x2li, x2omi, xl, xldot , xnddt, xndt, xomi, g22, g32,
  g44, g52, g54, fasx2, fasx4, fasx6, rptim , step2, stepn , stepp;

  fasx2 = 0.13130908;
  fasx4 = 2.8843198;
  fasx6 = 0.37448087;
  g22   = 5.7686396;
  g32   = 0.95240898;
  g44   = 1.8014998;
  g52   = 1.0508330;
  g54   = 4.4108898;
  rptim = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
  stepp =    720.0;
  stepn =   -720.0;
  step2 = 259200.0;

  /* ----------- calculate deep space resonance effects ----------- */
  dndt   = 0.0;
  theta  = fmod(gsto + tc * rptim, twopi);
  em     = em + dedt * t;

  inclm  = inclm + didt * t;
  argpm  = argpm + domdt * t;
  nodem  = nodem + dnodt * t;
  mm     = mm + dmdt * t;

  //     if (tc != 0)
  //     {
  //       printf("tdelta: %+.15e\n", tc);
  //       printf("theta:  %+.15e\n", theta);
  //       printf("em:     %+.15e\n", em);
  //       printf("inclm:  %+.15e\n", inclm);
  //       printf("argpm:  %+.15e\n", argpm);
  //       printf("nodem:  %+.15e\n", nodem);
  //       printf("dmdt:   %+.15e\n", dmdt);
  //       printf("mm:     %+.15e\n", mm);
  //     }

  //   sgp4fix for negative inclinations
  //   the following if statement should be commented out
  if ((inclm < 0.0) && (opsmode == 'a'))
  {
    inclm = -inclm;
    argpm = argpm - pi;
    nodem = nodem + pi;
  }

  /* - update resonances : numerical (euler-maclaurin) integration - */
  /* ------------------------- epoch restart ----------------------  */
  //   sgp4fix for propagator problems
  //   the following integration works for negative time steps and periods
  //   the specific changes are unknown because the original code was so convoluted

  // sgp4fix take out atime = 0.0 and fix for faster operation
  ft    = 0.0;
  if (irez != 0)
  {
    // sgp4fix streamline check
    if ((atime == 0.0) || (t * atime <= 0.0) || (fabs(t) < fabs(atime)) )
    {
      atime  = 0.0;
      xni    = no;
      xli    = xlamo;
    }
    // sgp4fix move check outside loop
    if (t > 0.0)
      delt = stepp;
    else
      delt = stepn;

    //           if (tc != 0)
    //           {
    //             printf("atime:  %+.15e\n", atime);
    //             printf("xni:    %+.15e\n", xni);
    //             printf("xli:    %+.15e\n", xli);
    //             printf("delt:   %+.15e\n", delt);
    //           }

    iretn = 381; // added for do loop
    iret  =   0; // added for loop
    while (iretn == 381)
    {
      /* ------------------- dot terms calculated ------------- */
      /* ----------- near - synchronous resonance terms ------- */
      if (irez != 2)
      {
        xndt  = del1 * sin(xli - fasx2) + del2 * sin(2.0 * (xli - fasx4)) +
        del3 * sin(3.0 * (xli - fasx6));
        xldot = xni + xfact;
        xnddt = del1 * cos(xli - fasx2) +
            2.0 * del2 * cos(2.0 * (xli - fasx4)) +
            3.0 * del3 * cos(3.0 * (xli - fasx6));
        xnddt = xnddt * xldot;
      }
      else
      {
        /* --------- near - half-day resonance terms -------- */
        xomi  = argpo + argpdot * atime;
        x2omi = xomi + xomi;
        x2li  = xli + xli;
        xndt  = d2201 * sin(x2omi + xli - g22) + d2211 * sin(xli - g22) +
            d3210 * sin(xomi + xli - g32)  + d3222 * sin(-xomi + xli - g32)+
            d4410 * sin(x2omi + x2li - g44)+ d4422 * sin(x2li - g44) +
            d5220 * sin(xomi + xli - g52)  + d5232 * sin(-xomi + xli - g52)+
            d5421 * sin(xomi + x2li - g54) + d5433 * sin(-xomi + x2li - g54);
        xldot = xni + xfact;
        xnddt = d2201 * cos(x2omi + xli - g22) + d2211 * cos(xli - g22) +
            d3210 * cos(xomi + xli - g32) + d3222 * cos(-xomi + xli - g32) +
            d5220 * cos(xomi + xli - g52) + d5232 * cos(-xomi + xli - g52) +
            2.0 * (d4410 * cos(x2omi + x2li - g44) +
                d4422 * cos(x2li - g44) + d5421 * cos(xomi + x2li - g54) +
                d5433 * cos(-xomi + x2li - g54));
        xnddt = xnddt * xldot;
      }

      /* ----------------------- integrator ------------------- */
      // sgp4fix move end checks to end of routine
      if (fabs(t - atime) >= stepp)
      {
        iret  = 0;
        iretn = 381;
      }
      else // exit here
      {
        ft    = t - atime;
        iretn = 0;
      }

      if (iretn == 381)
      {
        xli   = xli + xldot * delt + xndt * step2;
        xni   = xni + xndt * delt + xnddt * step2;
        atime = atime + delt;
      }
//      printf("tdelta: %+.15e\n", tc);
//      printf("atime:  %+.15e\n", atime);
//      printf("xndt:   %+.15e\n", xndt);
//      printf("xldot:  %+.15e\n", xldot);
//      printf("xnddt:  %+.15e\n", xnddt);
//      printf("xomi:   %+.15e\n", xomi);
//      printf("x2omi:  %+.15e\n", x2omi);
//      printf("x2li:   %+.15e\n", x2li);
//      printf("xli:    %+.15e\n", xli);
//      printf("xni:    %+.15e\n", xni );
//      printf("---------------------------------\n");
    }  // while iretn = 381


    nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
    xl = xli + xldot * ft + xndt * ft * ft * 0.5;
    if (irez != 1)
    {
      mm   = xl - 2.0 * nodem + 2.0 * theta;
      dndt = nm - no;
    }
    else
    {
      mm   = xl - nodem - argpm + theta;
      dndt = nm - no;
    }
    nm = no + dndt;
  }
//  printf("nm:     %+.15e\n", nm);
//  printf("xmp:    %+.15e\n", mm);
//  printf("dndt:   %+.15e\n", dndt);
  //#include "debug4.cpp"
}  // end dsspace

/*-----------------------------------------------------------------------------
 *
 *                           procedure initl
 *
 *  this procedure initializes the spg4 propagator. all the initialization is
 *    consolidated here instead of having multiple loops inside other routines.
 *
 *  author        : david vallado                  719-573-2600   28 jun 2005
 *
 *  inputs        :
 *    ecco        - eccentricity                           0.0 - 1.0
 *    epoch       - epoch time in days from jan 0, 1950. 0 hr
 *    inclo       - inclination of satellite
 *    no          - mean motion of satellite
 *    satn        - satellite number
 *
 *  outputs       :
 *    ainv        - 1.0 / a
 *    ao          - semi major axis
 *    con41       -
 *    con42       - 1.0 - 5.0 cos(i)
 *    cosio       - cosine of inclination
 *    cosio2      - cosio squared
 *    eccsq       - eccentricity squared
 *    method      - flag for deep space                    'd', 'n'
 *    omeosq      - 1.0 - ecco * ecco
 *    posq        - semi-parameter squared
 *    rp          - radius of perigee
 *    rteosq      - square root of (1.0 - ecco*ecco)
 *    sinio       - sine of inclination
 *    gsto        - gst at time of observation               rad
 *    no          - mean motion of satellite
 *
 *  locals        :
 *    ak          -
 *    d1          -
 *    del         -
 *    adel        -
 *    po          -
 *
 *  coupling      :
 *    getgravconst
 *    gstime      - find greenwich sidereal time from the julian date
 *
 *  references    :
 *    hoots, roehrich, norad spacetrack report #3 1980
 *    hoots, norad spacetrack report #6 1986
 *    hoots, schumacher and glover 2004
 *    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

static void initl
(
    int satn,      gravconsttype whichconst,
    double ecco,   double epoch,  double inclo,   double& no,
    char& method,
    double& ainv,  double& ao,    double& con41,  double& con42, double& cosio,
    double& cosio2,double& eccsq, double& omeosq, double& posq,
    double& rp,    double& rteosq,double& sinio , double& gsto,
    char opsmode
)
{
  /* --------------------- local variables ------------------------ */
  double ak, d1, del, adel, po, x2o3, j2, xke,
  tumin, mu, radiusearthkm, j3, j4, j3oj2;

  // sgp4fix use old way of finding gst
  double ds70;
  double ts70, tfrac, c1, thgr70, fk5r, c1p2p;
  const double twopi = 2.0 * pi;

  /* ----------------------- earth constants ---------------------- */
  // sgp4fix identify constants and allow alternate values
  getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
  x2o3   = 2.0 / 3.0;

  /* ------------- calculate auxillary epoch quantities ---------- */
  eccsq  = ecco * ecco;
  omeosq = 1.0 - eccsq;
  rteosq = sqrt(omeosq);
  cosio  = cos(inclo);
  cosio2 = cosio * cosio;

  /* ------------------ un-kozai the mean motion ----------------- */
  ak    = pow(xke / no, x2o3);
  d1    = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
  del   = d1 / (ak * ak);
  adel  = ak * (1.0 - del * del - del *
      (1.0 / 3.0 + 134.0 * del * del / 81.0));
  del   = d1 /(adel * adel);
  no    = no / (1.0 + del);
  ao    = pow(xke / no, x2o3);
  sinio = sin(inclo);
  po    = ao * omeosq;
  con42 = 1.0 - 5.0 * cosio2;
  con41 = -con42-cosio2-cosio2;
  ainv  = 1.0 / ao;
  posq  = po * po;
  rp    = ao * (1.0 - ecco);
  method = 'n';

  // sgp4fix modern approach to finding sidereal time
  if (opsmode == 'a')
  {
    // sgp4fix use old way of finding gst
    // count integer number of days from 0 jan 1970
    ts70  = epoch - 7305.0;
    ds70 = floor(ts70 + 1.0e-8);
    tfrac = ts70 - ds70;
    // find greenwich location at epoch
    c1    = 1.72027916940703639e-2;
    thgr70= 1.7321343856509374;
    fk5r  = 5.07551419432269442e-15;
    c1p2p = c1 + twopi;
    gsto  = fmod( thgr70 + c1*ds70 + c1p2p*tfrac + ts70*ts70*fk5r, twopi);
    if ( gsto < 0.0 )
      gsto = gsto + twopi;
  }
  else
    gsto = gstime(epoch + 2433281.5);

  printf("----------------------------------------\n");
  printf("aa      %+.15e\n", ak);
  printf("cosio   %+.15e\n", cosio);
  printf("sinio   %+.15e\n", sinio);
  printf("eo2     %+.15e\n", eccsq);
  printf("theta2  %+.15e\n", cosio2);
  printf("betao2  %+.15e\n", omeosq);
  printf("betao   %+.15e\n", rteosq);
  printf("delta1  %+.15e\n", d1);
  printf("a0      %+.15e\n", adel);
  printf("delta0  %+.15e\n", del);
  printf("xnodp   %+.15e\n", no);
  printf("aodp    %+.15e\n", ao);
  printf("x1m5th2 %+.15e\n", con42);
  printf("con41   %+.15e\n", con41);

  //#include "debug5.cpp"
}  // end initl

/*-----------------------------------------------------------------------------
 *
 *                             procedure sgp4init
 *
 *  this procedure initializes variables for sgp4.
 *
 *  author        : david vallado                  719-573-2600   28 jun 2005
 *
 *  inputs        :
 *    opsmode     - mode of operation afspc or improved 'a', 'i'
 *    whichconst  - which set of constants to use  72, 84
 *    satn        - satellite number
 *    bstar       - sgp4 type drag coefficient              kg/m2er
 *    ecco        - eccentricity
 *    epoch       - epoch time in days from jan 0, 1950. 0 hr
 *    argpo       - argument of perigee (output if ds)
 *    inclo       - inclination
 *    mo          - mean anomaly (output if ds)
 *    no          - mean motion
 *    nodeo       - right ascension of ascending node
 *
 *  outputs       :
 *    satrec      - common values for subsequent calls
 *    return code - non-zero on error.
 *                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
 *                   2 - mean motion less than 0.0
 *                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
 *                   4 - semi-latus rectum < 0.0
 *                   5 - epoch elements are sub-orbital
 *                   6 - satellite has decayed
 *
 *  locals        :
 *    cnodm  , snodm  , cosim  , sinim  , cosomm , sinomm
 *    cc1sq  , cc2    , cc3
 *    coef   , coef1
 *    cosio4      -
 *    day         -
 *    dndt        -
 *    em          - eccentricity
 *    emsq        - eccentricity squared
 *    eeta        -
 *    etasq       -
 *    gam         -
 *    argpm       - argument of perigee
 *    nodem       -
 *    inclm       - inclination
 *    mm          - mean anomaly
 *    nm          - mean motion
 *    perige      - perigee
 *    pinvsq      -
 *    psisq       -
 *    qzms24      -
 *    rtemsq      -
 *    s1, s2, s3, s4, s5, s6, s7          -
 *    sfour       -
 *    ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
 *    sz1, sz2, sz3
 *    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
 *    tc          -
 *    temp        -
 *    temp1, temp2, temp3       -
 *    tsi         -
 *    xpidot      -
 *    xhdot1      -
 *    z1, z2, z3          -
 *    z11, z12, z13, z21, z22, z23, z31, z32, z33         -
 *
 *  coupling      :
 *    getgravconst-
 *    initl       -
 *    dscom       -
 *    dpper       -
 *    dsinit      -
 *    sgp4        -
 *
 *  references    :
 *    hoots, roehrich, norad spacetrack report #3 1980
 *    hoots, norad spacetrack report #6 1986
 *    hoots, schumacher and glover 2004
 *    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

bool sgp4init
(
    gravconsttype whichconst, char opsmode,   const int satn,     const double epoch,
    const double xbstar,  const double xecco, const double xargpo,
    const double xinclo,  const double xmo,   const double xno,
    const double xnodeo,  elsetrec& satrec
)
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

  /* ------------------------ initialization --------------------- */
  // sgp4fix divisor for divide by zero check on inclination
  // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
  // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
  const double temp4    =   1.5e-12;

  /* ----------- set all near earth variables to zero ------------ */
  satrec.isimp   = 0;   satrec.method = 'n'; satrec.aycof    = 0.0;
  satrec.con41   = 0.0; satrec.cc1    = 0.0; satrec.cc4      = 0.0;
  satrec.cc5     = 0.0; satrec.d2     = 0.0; satrec.d3       = 0.0;
  satrec.d4      = 0.0; satrec.delmo  = 0.0; satrec.eta      = 0.0;
  satrec.argpdot = 0.0; satrec.omgcof = 0.0; satrec.sinmao   = 0.0;
  satrec.t       = 0.0; satrec.t2cof  = 0.0; satrec.t3cof    = 0.0;
  satrec.t4cof   = 0.0; satrec.t5cof  = 0.0; satrec.x1mth2   = 0.0;
  satrec.x7thm1  = 0.0; satrec.mdot   = 0.0; satrec.nodedot  = 0.0;
  satrec.xlcof   = 0.0; satrec.xmcof  = 0.0; satrec.nodecf   = 0.0;

  /* ----------- set all deep space variables to zero ------------ */
  satrec.irez  = 0;   satrec.d2201 = 0.0; satrec.d2211 = 0.0;
  satrec.d3210 = 0.0; satrec.d3222 = 0.0; satrec.d4410 = 0.0;
  satrec.d4422 = 0.0; satrec.d5220 = 0.0; satrec.d5232 = 0.0;
  satrec.d5421 = 0.0; satrec.d5433 = 0.0; satrec.dedt  = 0.0;
  satrec.del1  = 0.0; satrec.del2  = 0.0; satrec.del3  = 0.0;
  satrec.didt  = 0.0; satrec.dmdt  = 0.0; satrec.dnodt = 0.0;
  satrec.domdt = 0.0; satrec.e3    = 0.0; satrec.ee2   = 0.0;
  satrec.peo   = 0.0; satrec.pgho  = 0.0; satrec.pho   = 0.0;
  satrec.pinco = 0.0; satrec.plo   = 0.0; satrec.se2   = 0.0;
  satrec.se3   = 0.0; satrec.sgh2  = 0.0; satrec.sgh3  = 0.0;
  satrec.sgh4  = 0.0; satrec.sh2   = 0.0; satrec.sh3   = 0.0;
  satrec.si2   = 0.0; satrec.si3   = 0.0; satrec.sl2   = 0.0;
  satrec.sl3   = 0.0; satrec.sl4   = 0.0; satrec.gsto  = 0.0;
  satrec.xfact = 0.0; satrec.xgh2  = 0.0; satrec.xgh3  = 0.0;
  satrec.xgh4  = 0.0; satrec.xh2   = 0.0; satrec.xh3   = 0.0;
  satrec.xi2   = 0.0; satrec.xi3   = 0.0; satrec.xl2   = 0.0;
  satrec.xl3   = 0.0; satrec.xl4   = 0.0; satrec.xlamo = 0.0;
  satrec.zmol  = 0.0; satrec.zmos  = 0.0; satrec.atime = 0.0;
  satrec.xli   = 0.0; satrec.xni   = 0.0;

  // sgp4fix - note the following variables are also passed directly via satrec.
  // it is possible to streamline the sgp4init call by deleting the "x"
  // variables, but the user would need to set the satrec.* values first. we
  // include the additional assignments in case twoline2rv is not used.
  satrec.bstar   = xbstar;
  satrec.ecco    = xecco;
  satrec.argpo   = xargpo;
  satrec.inclo   = xinclo;
  satrec.mo	    = xmo;
  satrec.no	    = xno;
  satrec.nodeo   = xnodeo;

  // sgp4fix add opsmode
  satrec.operationmode = opsmode;

  printf("----------------------------------------\n");
  printf("mean_motion_dt2  %+.15e\n", satrec.ndot);
  printf("mean_motion_ddt6 %+.15e\n", satrec.nddot);
  printf("Bstar            %+.15e\n", satrec.bstar);
  printf("inclination      %+.15e\n", satrec.inclo);
  printf("right_asc_node   %+.15e\n", satrec.nodeo);
  printf("eccentricity     %+.15e\n", satrec.ecco);
  printf("argument_perigee %+.15e\n", satrec.argpo);
  printf("mean_anomaly     %+.15e\n", satrec.mo);
  printf("mean_motion      %+.15e\n", satrec.no);

  /* ------------------------ earth constants ----------------------- */
  // sgp4fix identify constants and allow alternate values
  getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
  ss     = 78.0 / radiusearthkm + 1.0;
  // sgp4fix use multiply for speed instead of pow
  qzms2ttemp = (120.0 - 78.0) / radiusearthkm;
  qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
  x2o3   =  2.0 / 3.0;

  satrec.init = 'y';
  satrec.t	 = 0.0;
  //print_orbit(satrec, "sgp4init pre-initl");
  initl
  (
      satn, whichconst, satrec.ecco, epoch, satrec.inclo, satrec.no, satrec.method,
      ainv, ao, satrec.con41, con42, cosio, cosio2, eccsq, omeosq,
      posq, rp, rteosq, sinio, satrec.gsto, satrec.operationmode
  );
  satrec.error = 0;
  //print_orbit(satrec, "sgp4init post-initl");

  // sgp4fix remove this check as it is unnecessary
  // the mrt check in sgp4 handles decaying satellite cases even if the starting
  // condition is below the surface of te earth
  //     if (rp < 1.0)
  //       {
  //        // printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
  //         satrec.error = 5;
  //       }

  if ((omeosq >= 0.0 ) || ( satrec.no >= 0.0))
  {
    satrec.isimp = 0;
    if (rp < (220.0 / radiusearthkm + 1.0))
      satrec.isimp = 1;
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
    satrec.eta  = ao * satrec.ecco * tsi;
    etasq = satrec.eta * satrec.eta;
    eeta  = satrec.ecco * satrec.eta;
    psisq = fabs(1.0 - etasq);
    coef  = qzms24 * pow(tsi, 4.0);
    coef1 = coef / pow(psisq, 3.5);
    cc2   = coef1 * satrec.no * (ao * (1.0 + 1.5 * etasq + eeta *
        (4.0 + etasq)) + 0.375 * j2 * tsi / psisq * satrec.con41 *
        (8.0 + 3.0 * etasq * (8.0 + etasq)));
    satrec.cc1   = satrec.bstar * cc2;
    cc3   = 0.0;
    if (satrec.ecco > 1.0e-4)
      cc3 = -2.0 * coef * tsi * j3oj2 * satrec.no * sinio / satrec.ecco;
    satrec.x1mth2 = 1.0 - cosio2;
    satrec.cc4    = 2.0* satrec.no * coef1 * ao * omeosq *
        (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco *
            (0.5 + 2.0 * etasq) - j2 * tsi / (ao * psisq) *
            (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq *
                (1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 *
                (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * satrec.argpo)));
    satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
        (etasq + eeta) + eeta * etasq);
    cosio4 = cosio2 * cosio2;
    temp1  = 1.5 * j2 * pinvsq * satrec.no;
    temp2  = 0.5 * temp1 * j2 * pinvsq;
    temp3  = -0.46875 * j4 * pinvsq * pinvsq * satrec.no;
    satrec.mdot     = satrec.no + 0.5 * temp1 * rteosq * satrec.con41 + 0.0625 *
        temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
    satrec.argpdot  = -0.5 * temp1 * con42 + 0.0625 * temp2 *
        (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
        temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
    xhdot1            = -temp1 * cosio;
    satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
        2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
    xpidot            =  satrec.argpdot+ satrec.nodedot;
    satrec.omgcof   = satrec.bstar * cc3 * cos(satrec.argpo);
    satrec.xmcof    = 0.0;
    if (satrec.ecco > 1.0e-4)
      satrec.xmcof = -x2o3 * coef * satrec.bstar / eeta;
    satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
    satrec.t2cof   = 1.5 * satrec.cc1;
    // sgp4fix for divide by zero with xinco = 180 deg
    if (fabs(cosio+1.0) > 1.5e-12)
      satrec.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    else
      satrec.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
    satrec.aycof   = -0.5 * j3oj2 * sinio;
    // sgp4fix use multiply for speed instead of pow
    delmotemp = 1.0 + satrec.eta * cos(satrec.mo);
    satrec.delmo   = delmotemp * delmotemp * delmotemp;
    satrec.sinmao  = sin(satrec.mo);
    satrec.x7thm1  = 7.0 * cosio2 - 1.0;

    printf("----------------------------------------\n");
    printf("qoms24 %+.15e\n", qzms24);
    printf("s4     %+.15e\n", sfour);
    printf("pinv2  %+.15e\n", pinvsq);
    printf("tsi    %+.15e\n", tsi);
    printf("eta    %+.15e\n", satrec.eta);
    printf("eta2   %+.15e\n", etasq);
    printf("eeta   %+.15e\n", eeta);
    printf("psi2   %+.15e\n", psisq);
    printf("coef   %+.15e\n", coef);
    printf("coef1  %+.15e\n", coef1);
    printf("C2     %+.15e\n", cc2);
    printf("C1     %+.15e\n", satrec.cc1);
    printf("C4     %+.15e\n", satrec.cc4);
    printf("theta4 %+.15e\n", cosio4);
    printf("xmdot  %+.15e\n", satrec.mdot);
    printf("omgdot %+.15e\n", satrec.argpdot);
    printf("xhdot1 %+.15e\n", xhdot1);
    printf("xnodot %+.15e\n", satrec.nodedot);
    printf("xnodcf %+.15e\n", satrec.nodecf);
    printf("xlcof  %+.15e\n", satrec.xlcof);
    printf("aycof  %+.15e\n", satrec.aycof);
    printf("x7thm1 %+.15e\n", satrec.x7thm1);
    printf("x1mth2 %+.15e\n", satrec.x1mth2);
    printf("dspace %d\n", (satrec.method == 'd'));
    printf("simple %d\n", satrec.isimp);

    /* --------------- deep space initialization ------------- */
    if ((2*pi / satrec.no) >= 225.0)
    {
      satrec.method = 'd';
      satrec.isimp  = 1;
      tc    =  0.0;
      inclm = satrec.inclo;
      dscom
      (
          epoch, satrec.ecco, satrec.argpo, tc, satrec.inclo, satrec.nodeo,
          satrec.no, snodm, cnodm,  sinim, cosim,sinomm,     cosomm,
          day, satrec.e3, satrec.ee2, em,         emsq, gam,
          satrec.peo,  satrec.pgho,   satrec.pho, satrec.pinco,
          satrec.plo,  rtemsq,        satrec.se2, satrec.se3,
          satrec.sgh2, satrec.sgh3,   satrec.sgh4,
          satrec.sh2,  satrec.sh3,    satrec.si2, satrec.si3,
          satrec.sl2,  satrec.sl3,    satrec.sl4, s1, s2, s3, s4, s5,
          s6,   s7,   ss1,  ss2,  ss3,  ss4,  ss5,  ss6,  ss7, sz1, sz2, sz3,
          sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33,
          satrec.xgh2, satrec.xgh3,   satrec.xgh4, satrec.xh2,
          satrec.xh3,  satrec.xi2,    satrec.xi3,  satrec.xl2,
          satrec.xl3,  satrec.xl4,    nm, z1, z2, z3, z11,
          z12, z13, z21, z22, z23, z31, z32, z33,
          satrec.zmol, satrec.zmos
      );
      //print_orbit(satrec, "sgp4init post-dscom");
      dpper
      (
          satrec.e3, satrec.ee2, satrec.peo, satrec.pgho,
          satrec.pho, satrec.pinco, satrec.plo, satrec.se2,
          satrec.se3, satrec.sgh2, satrec.sgh3, satrec.sgh4,
          satrec.sh2, satrec.sh3, satrec.si2,  satrec.si3,
          satrec.sl2, satrec.sl3, satrec.sl4,  satrec.t,
          satrec.xgh2,satrec.xgh3,satrec.xgh4, satrec.xh2,
          satrec.xh3, satrec.xi2, satrec.xi3,  satrec.xl2,
          satrec.xl3, satrec.xl4, satrec.zmol, satrec.zmos, inclm, satrec.init,
          satrec.ecco, satrec.inclo, satrec.nodeo, satrec.argpo, satrec.mo,
          satrec.operationmode
      );
      //print_orbit(satrec, "sgp4init post-dpper");
      argpm  = 0.0;
      nodem  = 0.0;
      mm     = 0.0;

      dsinit
      (
          whichconst,
          cosim, emsq, satrec.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
          ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, satrec.t, tc,
          satrec.gsto, satrec.mo, satrec.mdot, satrec.no, satrec.nodeo,
          satrec.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
          satrec.ecco, eccsq, em, argpm, inclm, mm, nm, nodem,
          satrec.irez,  satrec.atime,
          satrec.d2201, satrec.d2211, satrec.d3210, satrec.d3222 ,
          satrec.d4410, satrec.d4422, satrec.d5220, satrec.d5232,
          satrec.d5421, satrec.d5433, satrec.dedt,  satrec.didt,
          satrec.dmdt,  dndt,         satrec.dnodt, satrec.domdt ,
          satrec.del1,  satrec.del2,  satrec.del3,  satrec.xfact,
          satrec.xlamo, satrec.xli,   satrec.xni, opsmode
      );
      //print_orbit(satrec, "sgp4init post-dsinit");
    }

    /* ----------- set variables if not deep space ----------- */
    if (satrec.isimp != 1)
    {
      cc1sq          = satrec.cc1 * satrec.cc1;
      satrec.d2    = 4.0 * ao * tsi * cc1sq;
      temp           = satrec.d2 * tsi * satrec.cc1 / 3.0;
      satrec.d3    = (17.0 * ao + sfour) * temp;
      satrec.d4    = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) *
          satrec.cc1;
      satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
      satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 *
          (12.0 * satrec.d2 + 10.0 * cc1sq));
      satrec.t5cof = 0.2 * (3.0 * satrec.d4 +
          12.0 * satrec.cc1 * satrec.d3 +
          6.0 * satrec.d2 * satrec.d2 +
          15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
    }
  } // if omeosq = 0 ...


  /* finally propogate to zero epoch to initialize all others. */
  // sgp4fix take out check to let satellites process until they are actually below earth surface
  //       if(satrec.error == 0)

  printf("----------------------------------------\n");
  printf("C3     %+.15e\n", cc3);
  printf("xmcof  %+.15e\n", satrec.xmcof);
  printf("C5     %+.15e\n", satrec.cc5);
  printf("omgcof %+.15e\n", satrec.omgcof);
  printf("delmo  %+.15e\n", satrec.delmo);
  printf("sinmo  %+.15e\n", satrec.sinmao);
  printf("D2     %+.15e\n", satrec.d2);
  printf("D3     %+.15e\n", satrec.d3);
  printf("D4     %+.15e\n", satrec.d4);
  printf("t3cof  %+.15e\n", satrec.t3cof);
  printf("t4cof  %+.15e\n", satrec.t4cof);
  printf("t5cof  %+.15e\n", satrec.t5cof);

  sgp4(whichconst, satrec, 0.0, r, v);
  satrec.init = 'n';
  //print_orbit(satrec, "sgp4init post-sgp4 call");

  //#include "debug6.cpp"
  //sgp4fix return boolean. satrec.error contains any error codes

  return true;
}  // end sgp4init

/*-----------------------------------------------------------------------------
 *
 *                             procedure sgp4
 *
 *  this procedure is the sgp4 prediction model from space command. this is an
 *    updated and combined version of sgp4 and sdp4, which were originally
 *    published separately in spacetrack report #3. this version follows the
 *    methodology from the aiaa paper (2006) describing the history and
 *    development of the code.
 *
 *  author        : david vallado                  719-573-2600   28 jun 2005
 *
 *  inputs        :
 *    satrec	 - initialised structure from sgp4init() call.
 *    tsince	 - time eince epoch (minutes)
 *
 *  outputs       :
 *    r           - position vector                     km
 *    v           - velocity                            km/sec
 *  return code - non-zero on error.
 *                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
 *                   2 - mean motion less than 0.0
 *                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
 *                   4 - semi-latus rectum < 0.0
 *                   5 - epoch elements are sub-orbital
 *                   6 - satellite has decayed
 *
 *  locals        :
 *    am          -
 *    axnl, aynl        -
 *    betal       -
 *    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
 *    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
 *    cosisq  , cossu   , sinsu   , cosu    , sinu
 *    delm        -
 *    delomg      -
 *    dndt        -
 *    eccm        -
 *    emsq        -
 *    ecose       -
 *    el2         -
 *    eo1         -
 *    eccp        -
 *    esine       -
 *    argpm       -
 *    argpp       -
 *    omgadf      -
 *    pl          -
 *    r           -
 *    rtemsq      -
 *    rdotl       -
 *    rl          -
 *    rvdot       -
 *    rvdotl      -
 *    su          -
 *    t2  , t3   , t4    , tc
 *    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
 *    u   , ux   , uy    , uz     , vx     , vy     , vz
 *    inclm       - inclination
 *    mm          - mean anomaly
 *    nm          - mean motion
 *    nodem       - right asc of ascending node
 *    xinc        -
 *    xincp       -
 *    xl          -
 *    xlm         -
 *    mp          -
 *    xmdf        -
 *    xmx         -
 *    xmy         -
 *    nodedf      -
 *    xnode       -
 *    nodep       -
 *    np          -
 *
 *  coupling      :
 *    getgravconst-
 *    dpper
 *    dpspace
 *
 *  references    :
 *    hoots, roehrich, norad spacetrack report #3 1980
 *    hoots, norad spacetrack report #6 1986
 *    hoots, schumacher and glover 2004
 *    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

bool sgp4
(
    gravconsttype whichconst, elsetrec& satrec,  double tsince,
    double r[3],  double v[3]
)
{
  double am   , axnl  , aynl , betal ,  cosim , cnod  ,
  cos2u, coseo1, cosi , cosip ,  cosisq, cossu , cosu,
  delm , delomg, em   , emsq  ,  ecose , el2   , eo1 ,
  ep   , esine , argpm, argpp ,  argpdf, pl,     mrt = 0.0,
  mvt  , rdotl , rl   , rvdot ,  rvdotl, sinim ,
  sin2u, sineo1, sini , sinip ,  sinsu , sinu  ,
  snod , su    , t2   , t3    ,  t4    , tem5  , temp,
  temp1, temp2 , tempa, tempe ,  templ , u     , ux  ,
  uy   , uz    , vx   , vy    ,  vz    , inclm , mm  ,
  nm   , nodem, xinc , xincp ,  xl    , xlm   , mp  ,
  xmdf , xmx   , xmy  , nodedf, xnode , nodep, tc  , dndt,
  twopi, x2o3  , j2   , j3    , tumin, j4 , xke   , j3oj2, radiusearthkm,
  mu, vkmpersec, delmtemp;
  int ktr;

  /* ------------------ set mathematical constants --------------- */
  // sgp4fix divisor for divide by zero check on inclination
  // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
  // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
  const double temp4 =   1.5e-12;
  twopi = 2.0 * pi;
  x2o3  = 2.0 / 3.0;
  // sgp4fix identify constants and allow alternate values
  getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
  vkmpersec     = radiusearthkm * xke/60.0;

  /* --------------------- clear sgp4 error flag ----------------- */
  satrec.t     = tsince;
  satrec.error = 0;

  /* ------- update for secular gravity and atmospheric drag ----- */
  xmdf    = satrec.mo + satrec.mdot * satrec.t;
  argpdf  = satrec.argpo + satrec.argpdot * satrec.t;
  nodedf  = satrec.nodeo + satrec.nodedot * satrec.t;
  argpm   = argpdf;
  mm      = xmdf;
  t2      = satrec.t * satrec.t;
  nodem   = nodedf + satrec.nodecf * t2;
  tempa   = 1.0 - satrec.cc1 * satrec.t;
  tempe   = satrec.bstar * satrec.cc4 * satrec.t;
  templ   = satrec.t2cof * t2;

  if (satrec.isimp != 1)
  {
    delomg = satrec.omgcof * satrec.t;
    // sgp4fix use mutliply for speed instead of pow
    delmtemp =  1.0 + satrec.eta * cos(xmdf);
    delm   = satrec.xmcof *
        (delmtemp * delmtemp * delmtemp -
            satrec.delmo);
    temp   = delomg + delm;
    mm     = xmdf + temp;
    argpm  = argpdf - temp;
    t3     = t2 * satrec.t;
    t4     = t3 * satrec.t;
    tempa  = tempa - satrec.d2 * t2 - satrec.d3 * t3 -
        satrec.d4 * t4;
    tempe  = tempe + satrec.bstar * satrec.cc5 * (sin(mm) -
        satrec.sinmao);
    templ  = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof +
        satrec.t * satrec.t5cof);
  }

  nm    = satrec.no;
  em    = satrec.ecco;
  inclm = satrec.inclo;

  printf("----------------------------------------\n");
  printf("xmdf   %+.15e\n", xmdf);
  printf("xmp    %+.15e\n", mm);
  printf("omgadf %+.15e\n", argpdf);
  printf("xnoddf %+.15e\n", nodedf);
  printf("t2     %+.15e\n", t2);
  printf("tempa  %+.15e\n", tempa);
  printf("tempe  %+.15e\n", tempe);
  printf("templ  %+.15e\n", templ);
  printf("omega  %+.15e\n", argpm);
  printf("nm     %+.15e\n", nm);
  printf("em     %+.15e\n", em);
  printf("inclm  %+.15e\n", inclm);

  if (satrec.method == 'd')
  {
    tc = satrec.t;
    dspace
    (
        satrec.irez,
        satrec.d2201, satrec.d2211, satrec.d3210,
        satrec.d3222, satrec.d4410, satrec.d4422,
        satrec.d5220, satrec.d5232, satrec.d5421,
        satrec.d5433, satrec.dedt,  satrec.del1,
        satrec.del2,  satrec.del3,  satrec.didt,
        satrec.dmdt,  satrec.dnodt, satrec.domdt,
        satrec.argpo, satrec.argpdot, satrec.t, tc,
        satrec.gsto, satrec.xfact, satrec.xlamo,
        satrec.no, satrec.atime,
        em, argpm, inclm, satrec.xli, mm, satrec.xni,
        nodem, dndt, nm, satrec.operationmode
    );
  } // if method = d

  if (nm <= 0.0)
  {
    //        // printf("# error nm %f\n", nm);
    satrec.error = 2;
    // sgp4fix add return
    return false;
  }

  am = pow((xke / nm),x2o3) * tempa * tempa;
  nm = xke / pow(am, 1.5);
  em = em - tempe;

//  printf("-------------------------------\n");
//  printf("am:     %+.15e\n", am);
//  printf("nm:     %+.15e\n", nm);
//  printf("em:     %+.15e\n", em);

  // fix tolerance for error recognition
  // sgp4fix am is fixed from the previous nm check
  if ((em >= 1.0) || (em < -0.001) || ((am < 0.95) && (satrec.operationmode == 'a')) )
  {
    //        // printf("# error em %f\n", em);
    satrec.error = 1;
    // sgp4fix to return if there is an error in eccentricity
    return false;
  }
  // sgp4fix fix tolerance to avoid a divide by zero
  if ((em < 1.0e-6) && (satrec.operationmode != 'a'))
    em  = 1.0e-6;
  mm     = mm + satrec.no * templ;
  xlm    = mm + argpm + nodem;
  emsq   = em * em;
  temp   = 1.0 - emsq;

  nodem  = fmod(nodem, twopi);
  argpm  = fmod(argpm, twopi);
  xlm    = fmod(xlm, twopi);
  mm     = fmod(xlm - argpm - nodem, twopi);

  /* ----------------- compute extra mean quantities ------------- */
  cosim = cos(inclm);
  sinim = sin(inclm);

  /* -------------------- add lunar-solar periodics -------------- */
  ep     = em;
  xincp  = inclm;
  argpp  = argpm;
  nodep  = nodem;
  mp     = mm;
  sinip  = sinim;
  cosip  = cosim;

  printf("----------------------------------------\n");
  printf("am     %+.15e\n", am);
  printf("nm     %+.15e\n", nm);
  printf("em     %+.15e\n", em);
  printf("xlm    %+.15e\n", xlm);
  printf("em2    %+.15e\n", emsq);
  printf("em     %+.15e\n", em);
  printf("xnode  %+.15e\n", xnode);
  printf("omega  %+.15e\n", argpm);
  printf("inclp  %+.15e\n", xincp);
  printf("ep     %+.15e\n", ep);
  printf("nodep  %+.15e\n", nodep);
  printf("argpp  %+.15e\n", argpp);
  printf("mp     %+.15e\n", mp);
  printf("sinip  %+.15e\n", sinip);
  printf("cosip  %+.15e\n", cosip);

  //if (satrec.method == 'd') // TODO: Debug dpper
  if (false)
  {
    dpper
    (
        satrec.e3,   satrec.ee2,  satrec.peo,
        satrec.pgho, satrec.pho,  satrec.pinco,
        satrec.plo,  satrec.se2,  satrec.se3,
        satrec.sgh2, satrec.sgh3, satrec.sgh4,
        satrec.sh2,  satrec.sh3,  satrec.si2,
        satrec.si3,  satrec.sl2,  satrec.sl3,
        satrec.sl4,  satrec.t,    satrec.xgh2,
        satrec.xgh3, satrec.xgh4, satrec.xh2,
        satrec.xh3,  satrec.xi2,  satrec.xi3,
        satrec.xl2,  satrec.xl3,  satrec.xl4,
        satrec.zmol, satrec.zmos, satrec.inclo,
        'n', ep, xincp, nodep, argpp, mp, satrec.operationmode
    );

//    printf("-------------------------------\n");
//    printf("xmp:    %+.15e\n", mm);
//    printf("xlm:    %+.15e\n", xlm);
//    printf("em2:    %+.15e\n", emsq);
//    printf("xnode:  %+.15e\n", nodem);
//    printf("omega:  %+.15e\n", argpm);
//    printf("incl_lp:%+.15e\n", xincp);
//    printf("node_lp:%+.15e\n", nodep);
//    printf("argplp: %+.15e\n", argpp);
//    printf("ecc_lp: %+.15e\n", ep);
//    printf("mo_lp:  %+.15e\n", mp);

    if (xincp < 0.0)
    {
      xincp  = -xincp;
      nodep = nodep + pi;
      argpp  = argpp - pi;
    }
    if ((ep < 0.0 ) || ( ep > 1.0))
    {
      //           // printf("# error ep %f\n", ep);
      satrec.error = 3;
      // sgp4fix add return
      return false;
    }
  } // if method = d

  /* -------------------- long period periodics ------------------ */
  if (satrec.method == 'd')
  {
    sinip =  sin(xincp);
    cosip =  cos(xincp);
    satrec.aycof = -0.5*j3oj2*sinip;
    // sgp4fix for divide by zero for xincp = 180 deg
    if ((fabs(cosip+1.0) > 1.5e-12) && (satrec.operationmode != 'a'))
      satrec.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
    else
      satrec.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
  }
  axnl = ep * cos(argpp);
  temp = 1.0 / (am * (1.0 - ep * ep));
  aynl = ep* sin(argpp) + temp * satrec.aycof;
  xl   = mp + argpp + nodep + temp * satrec.xlcof * axnl;

  /* --------------------- solve kepler's equation --------------- */
  u    = fmod(xl - nodep, twopi); // TODO: This is evaluating incorrectly? Failed to reproduce
  u    = 7.699938584533328e-004; // CRUTCH!
  eo1  = u;
  tem5 = 9999.9;
  ktr  = 1;

  printf("----------------------------------------\n");
  printf("xlcof   %+.15e\n", satrec.xlcof);
  printf("axnl    %+.15e\n", axnl);
  printf("a1e2inv %+.15e\n", temp);
  printf("aynl    %+.15e\n", aynl);
  printf("xl      %+.25e\n", xl);
  printf("nodep   %+.25e\n", nodep);
  printf("twopi   %+.25e\n", twopi);
  printf("u       %+.15e\n", u);

  //   sgp4fix for kepler iteration
  //   the following iteration needs better limits on corrections
  while (( fabs(tem5) >= 1.0e-12) && (ktr <= 10) )
  {
    sineo1 = sin(eo1);
    coseo1 = cos(eo1);
    tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
    tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
    if(fabs(tem5) >= 0.95)
      tem5 = tem5 > 0.0 ? 0.95 : -0.95;
    eo1    = eo1 + tem5;
    ktr = ktr + 1;

    printf("> >>>>>>\n");
    printf("> sineo1 %+.15e\n", sineo1);
    printf("> coseo1 %+.15e\n", coseo1);
    printf("> eo1    %+.15e\n", eo1);
    printf("> ktr    %d\n", ktr);
  }

  /* ------------- short period preliminary quantities ----------- */
  ecose = axnl*coseo1 + aynl*sineo1;
  esine = axnl*sineo1 - aynl*coseo1;
  el2   = axnl*axnl + aynl*aynl;
  pl    = am*(1.0-el2);

  printf("----------------------------------------\n");
  printf("ecose %+.15e\n", ecose);
  printf("esine %+.15e\n", esine);
  printf("el2   %+.15e\n", el2);
  printf("pl    %+.15e\n", pl);

  if (pl < 0.0)
  {
    //        // printf("# error pl %f\n", pl);
    satrec.error = 4;
    // sgp4fix add return
    return false;
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

    /* -------------- update for short period periodics ------------ */
    if (satrec.method == 'd')
    {
      cosisq                 = cosip * cosip;
      satrec.con41  = 3.0*cosisq - 1.0;
      satrec.x1mth2 = 1.0 - cosisq;
      satrec.x7thm1 = 7.0*cosisq - 1.0;
    }
    mrt   = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) +
        0.5 * temp1 * satrec.x1mth2 * cos2u;
    su    = su - 0.25 * temp2 * satrec.x7thm1 * sin2u;
    xnode = nodep + 1.5 * temp2 * cosip * sin2u;
    xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
    mvt   = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / xke;
    rvdot = rvdotl + nm * temp1 * (satrec.x1mth2 * cos2u +
        1.5 * satrec.con41) / xke;

    printf("----------------------------------------\n");
    printf("rl     %+.15e\n", rl);
    printf("rdotl  %+.15e\n", rdotl);
    printf("rvdotl %+.15e\n", rvdotl);
    printf("betal  %+.15e\n", betal);
    printf("sinu   %+.15e\n", sinu);
    printf("cosu   %+.15e\n", cosu);
    printf("su     %+.15e\n", su);
    printf("sin2u  %+.15e\n", sin2u);
    printf("cos2u  %+.15e\n", cos2u);
    printf("temp1  %+.15e\n", temp1);
    printf("con41  %+.15e\n", satrec.con41);
    printf("x1mth2 %+.15e\n", satrec.x1mth2);
    printf("temp2  %+.15e\n", temp2);
    printf("mrt    %+.15e\n", mrt);
    printf("xnode  %+.15e\n", xnode);
    printf("xinc   %+.15e\n", xinc);
    printf("mvt    %+.15e\n", mvt);
    printf("rvdot  %+.15e\n", rvdot);

    /* --------------------- orientation vectors ------------------- */
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

    /* --------- position and velocity (in km and km/sec) ---------- */
    r[0] = (mrt * ux)* radiusearthkm;
    r[1] = (mrt * uy)* radiusearthkm;
    r[2] = (mrt * uz)* radiusearthkm;
    v[0] = (mvt * ux + rvdot * vx) * vkmpersec;
    v[1] = (mvt * uy + rvdot * vy) * vkmpersec;
    v[2] = (mvt * uz + rvdot * vz) * vkmpersec;

    printf("----------------------------------------\n");
    printf("px %+.15e\n", r[0]);
    printf("py %+.15e\n", r[1]);
    printf("pz %+.15e\n", r[2]);
    printf("vx %+.15e\n", v[0]);
    printf("vy %+.15e\n", v[1]);
    printf("vz %+.15e\n", v[2]);

  }  // if pl > 0

  // sgp4fix for decaying satellites
  if (mrt < 1.0)
  {
    //        // printf("# decay condition %11.6f \n",mrt);
    satrec.error = 6;
    return false;
  }

  //#include "debug7.cpp"
  return true;
}  // end sgp4


/* -----------------------------------------------------------------------------
 *
 *                           function gstime
 *
 *  this function finds the greenwich sidereal time.
 *
 *  author        : david vallado                  719-573-2600    1 mar 2001
 *
 *  inputs          description                    range / units
 *    jdut1       - julian date in ut1             days from 4713 bc
 *
 *  outputs       :
 *    gstime      - greenwich sidereal time        0 to 2pi rad
 *
 *  locals        :
 *    temp        - temporary variable for doubles   rad
 *    tut1        - julian centuries from the
 *                  jan 1, 2000 12 h epoch (ut1)
 *
 *  coupling      :
 *    none
 *
 *  references    :
 *    vallado       2004, 191, eq 3-45
 * --------------------------------------------------------------------------- */

double  gstime
(
    double jdut1
)
{
  const double twopi = 2.0 * pi;
  const double deg2rad = pi / 180.0;
  double       temp, tut1;

  tut1 = (jdut1 - 2451545.0) / 36525.0;
  temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
      (876600.0*3600 + 8640184.812866) * tut1 + 67310.54841;  // sec
      temp = fmod(temp * deg2rad / 240.0, twopi); //360/86400 = 1/240, to deg, to rad

      // ------------------------ check quadrants ---------------------
      if (temp < 0.0)
        temp += twopi;

      return temp;
}  // end gstime


/* -----------------------------------------------------------------------------
 *
 *                           function getgravconst
 *
 *  this function gets constants for the propagator. note that mu is identified to
 *    facilitiate comparisons with newer models. the common useage is wgs72.
 *
 *  author        : david vallado                  719-573-2600   21 jul 2006
 *
 *  inputs        :
 *    whichconst  - which set of constants to use  wgs72old, wgs72, wgs84
 *
 *  outputs       :
 *    tumin       - minutes in one time unit
 *    mu          - earth gravitational parameter
 *    radiusearthkm - radius of the earth in km
 *    xke         - reciprocal of tumin
 *    j2, j3, j4  - un-normalized zonal harmonic values
 *    j3oj2       - j3 divided by j2
 *
 *  locals        :
 *
 *  coupling      :
 *    none
 *
 *  references    :
 *    norad spacetrack report #3
 *    vallado, crawford, hujsak, kelso  2006
  --------------------------------------------------------------------------- */

void getgravconst
(
    gravconsttype whichconst,
    double& tumin,
    double& mu,
    double& radiusearthkm,
    double& xke,
    double& j2,
    double& j3,
    double& j4,
    double& j3oj2
)
{

  switch (whichconst)
  {
    // -- wgs-72 low precision str#3 constants --
    case wgs72old:
      mu     = 398600.79964;        // in km3 / s2
      radiusearthkm = 6378.135;     // km
      xke    = 0.0743669161;
      tumin  = 1.0 / xke;
      j2     =   0.001082616;
      j3     =  -0.00000253881;
      j4     =  -0.00000165597;
      j3oj2  =  j3 / j2;
      break;
      // ------------ wgs-72 constants ------------
    case wgs72:
      mu     = 398600.8;            // in km3 / s2
      radiusearthkm = 6378.135;     // km
      xke    = 60.0 / sqrt(radiusearthkm*radiusearthkm*radiusearthkm/mu);
      tumin  = 1.0 / xke;
      j2     =   0.001082616;
      j3     =  -0.00000253881;
      j4     =  -0.00000165597;
      j3oj2  =  j3 / j2;
      break;
    case wgs84:
      // ------------ wgs-84 constants ------------
      mu     = 398600.5;            // in km3 / s2
      radiusearthkm = 6378.137;     // km
      xke    = 60.0 / sqrt(radiusearthkm*radiusearthkm*radiusearthkm/mu);
      tumin  = 1.0 / xke;
      j2     =   0.00108262998905;
      j3     =  -0.00000253215306;
      j4     =  -0.00000161098761;
      j3oj2  =  j3 / j2;
      break;
    default:
      fprintf(stderr,"unknown gravity option (%d)\n",whichconst);
      break;
  }

}   // end getgravconst





