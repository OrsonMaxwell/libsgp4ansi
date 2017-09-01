#include <stdio.h>
#include <time.h>
#include "types.h"

void print_orbit(orbit* sat, char* caption)
{
  FILE *f = fopen("orbit.txt", "a");
  fprintf(f, "-------- %s --------\n", caption);
  fprintf(f, "-- TLE portion\n");
  fprintf(f, "number:\t\t%d\n", sat->number);
  fprintf(f, "epoch:\t\t%f\n", ctime(&sat->epoch));
  fprintf(f, "epoch_ms:\t%d\n", sat->epoch_ms);
  fprintf(f, "julepoch:\t%f\n", sat->julepoch);
  fprintf(f, "nprimediv2:\t%f\n", sat->nprimediv2);
  fprintf(f, "ndprimediv6:\t%f\n", sat->ndprimediv6);
  fprintf(f, "Bstar:\t\t%f\n", sat->Bstar);
  fprintf(f, "i:\t\t%f\n", sat->i);
  fprintf(f, "alpha:\t\t%f\n", sat->alpha);
  fprintf(f, "e:\t\t%f\n", sat->e);
  fprintf(f, "omega:\t\t%f\n", sat->omega);
  fprintf(f, "Mo:\t\t%f\n", sat->Mo);
  fprintf(f, "no:\t\t%f\n", sat->no);
  fprintf(f, "-- Standard elements\n");
  fprintf(f, "a:\t\t%f\n", sat->a);
  fprintf(f, "altapoR:\t\t%f\n", sat->altapoR);
  fprintf(f, "altperR:\t\t%f\n", sat->altperR);
  fprintf(f, "aycof:\t\t%f\n", sat->aycof);
  fprintf(f, "C1:\t\t%f\n", sat->C1);
  fprintf(f, "C4:\t\t%f\n", sat->C4);
  fprintf(f, "C5:\t\t%f\n", sat->C5);
  fprintf(f, "con41:\t\t%f\n", sat->con41);
  fprintf(f, "d2:\t\t%f\n", sat->d2);
  fprintf(f, "d3:\t\t%f\n", sat->d3);
  fprintf(f, "d4:\t\t%f\n", sat->d4);
  fprintf(f, "delMo:\t\t%f\n", sat->delMo);
  fprintf(f, "eta:\t\t%f\n", sat->eta);
  fprintf(f, "mdot:\t\t%f\n", sat->mdot);
  fprintf(f, "nodecf:\t\t%f\n", sat->nodecf);
  fprintf(f, "nodedot:\t%f\n", sat->nodedot);
  fprintf(f, "omegaprime:\t%f\n", sat->omegaprime);
  fprintf(f, "omgcof:\t\t%f\n", sat->omgcof);
  fprintf(f, "sinMo:\t\t%f\n", sat->sinMo);
  fprintf(f, "t2cof:\t\t%f\n", sat->t2cof);
  fprintf(f, "t3cof:\t\t%f\n", sat->t3cof);
  fprintf(f, "t4cof:\t\t%f\n", sat->t4cof);
  fprintf(f, "t5cof:\t\t%f\n", sat->t5cof);
  fprintf(f, "x1mth2:\t\t%f\n", sat->x1mth2);
  fprintf(f, "x7thm1:\t\t%f\n", sat->x7thm1);
  fprintf(f, "xlcof:\t\t%f\n", sat->xlcof);
  fprintf(f, "xmcof:\t\t%f\n", sat->xmcof);
  fprintf(f, "GSTo:\t\t%f\n", sat->GSTo);
  fprintf(f, "-- Deepspace elements\n");
  fprintf(f, "e3:\t\t%f\n", sat->e3);
  fprintf(f, "ee2:\t\t%f\n", sat->ee2);
  fprintf(f, "peo:\t\t%f\n", sat->peo);
  fprintf(f, "pgho:\t\t%f\n", sat->pgho);
  fprintf(f, "pho:\t\t%f\n", sat->pho);
  fprintf(f, "pinco:\t\t%f\n", sat->pinco);
  fprintf(f, "plo:\t\t%f\n", sat->plo);
  fprintf(f, "se2:\t\t%f\n", sat->se2);
  fprintf(f, "se3:\t\t%f\n", sat->se3);
  fprintf(f, "sgh2:\t\t%f\n", sat->sgh2);
  fprintf(f, "sgh3:\t\t%f\n", sat->sgh3);
  fprintf(f, "sgh4:\t\t%f\n", sat->sgh4);
  fprintf(f, "sh2:\t\t%f\n", sat->sh2);
  fprintf(f, "sh3:\t\t%f\n", sat->sh3);
  fprintf(f, "si2:\t\t%f\n", sat->si2);
  fprintf(f, "si3:\t\t%f\n", sat->si3);
  fprintf(f, "sl2:\t\t%f\n", sat->sl2);
  fprintf(f, "sl3:\t\t%f\n", sat->sl3);
  fprintf(f, "sl4:\t\t%f\n", sat->sl4);
  fprintf(f, "xgh2:\t\t%f\n", sat->xgh2);
  fprintf(f, "xgh3:\t\t%f\n", sat->xgh3);
  fprintf(f, "xgh4:\t\t%f\n", sat->xgh4);
  fprintf(f, "xh2:\t\t%f\n", sat->xh2);
  fprintf(f, "xh3:\t\t%f\n", sat->xh3);
  fprintf(f, "xi2:\t\t%f\n", sat->xi2);
  fprintf(f, "xi3:\t\t%f\n", sat->xi3);
  fprintf(f, "xl2:\t\t%f\n", sat->xl2);
  fprintf(f, "xl3:\t\t%f\n", sat->xl3);
  fprintf(f, "xl4:\t\t%f\n", sat->xl4);
  fprintf(f, "zmol:\t\t%f\n", sat->zmol);
  fprintf(f, "zmos:\t\t%f\n", sat->zmos);
  fprintf(f, "-- Flags\n");
  fprintf(f, "Deepsp:\t\t%d\n", sat->isdeepspace);
  fprintf(f, "Loworb:\t\t%d\n", sat->islowperigee);
  fclose(f);
}

void print_elsetrec(elsetrec* sat, char* caption)
{
  FILE *f = fopen("elsetrec.txt", "a");
  fprintf(f, "-------- %s --------\n", caption);
  fprintf(f, "-- TLE portion\n");
  fprintf(f, "number:\t\t%d\n", sat->satnum);
  fprintf(f, "epochyr:\t%d\n", sat->epochyr);
  fprintf(f, "epochdays:\t%f\n", sat->epochdays);
  fprintf(f, "jdsatepoch:\t%f\n", sat->jdsatepoch);
  fprintf(f, "nprimediv2:\t%f\n", sat->ndot);
  fprintf(f, "ndprimediv6:\t%f\n", sat->nddot);
  fprintf(f, "Bstar:\t\t%f\n", sat->bstar);
  fprintf(f, "i:\t\t%f\n", sat->inclo);
  fprintf(f, "alpha:\t\t%f\n", sat->nodeo);
  fprintf(f, "e:\t\t%f\n", sat->ecco);
  fprintf(f, "omega:\t\t%f\n", sat->argpo);
  fprintf(f, "Mo:\t\t%f\n", sat->mo);
  fprintf(f, "no:\t\t%f\n", sat->no);
  fprintf(f, "-- Standard elements\n");
  fprintf(f, "a:\t\t%f\n", sat->a);
  fprintf(f, "altapoR:\t\t%f\n", sat->alta);
  fprintf(f, "altperR:\t\t%f\n", sat->altp);
  fprintf(f, "aycof:\t\t%f\n", sat->aycof);
  fprintf(f, "C1:\t\t%f\n", sat->cc1);
  fprintf(f, "C4:\t\t%f\n", sat->cc4);
  fprintf(f, "C5:\t\t%f\n", sat->cc5);
  fprintf(f, "con41:\t\t%f\n", sat->con41);
  fprintf(f, "d2:\t\t%f\n", sat->d2);
  fprintf(f, "d3:\t\t%f\n", sat->d3);
  fprintf(f, "d4:\t\t%f\n", sat->d4);
  fprintf(f, "delMo:\t\t%f\n", sat->delmo);
  fprintf(f, "eta:\t\t%f\n", sat->eta);
  fprintf(f, "mdot:\t\t%f\n", sat->mdot);
  fprintf(f, "nodecf:\t\t%f\n", sat->nodecf);
  fprintf(f, "nodedot:\t%f\n", sat->nodedot);
  fprintf(f, "omegaprime:\t%f\n", sat->argpdot);
  fprintf(f, "omgcof:\t\t%f\n", sat->omgcof);
  fprintf(f, "sinMo:\t\t%f\n", sat->sinmao);
  fprintf(f, "t2cof:\t\t%f\n", sat->t2cof);
  fprintf(f, "t3cof:\t\t%f\n", sat->t3cof);
  fprintf(f, "t4cof:\t\t%f\n", sat->t4cof);
  fprintf(f, "t5cof:\t\t%f\n", sat->t5cof);
  fprintf(f, "x1mth2:\t\t%f\n", sat->x1mth2);
  fprintf(f, "x7thm1:\t\t%f\n", sat->x7thm1);
  fprintf(f, "xlcof:\t\t%f\n", sat->xlcof);
  fprintf(f, "xmcof:\t\t%f\n", sat->xmcof);
  fprintf(f, "GSTo:\t\t%f\n", sat->gsto);
  fprintf(f, "-- Deepspace elements\n");
  fprintf(f, "e3:\t\t%f\n", sat->e3);
  fprintf(f, "ee2:\t\t%f\n", sat->ee2);
  fprintf(f, "peo:\t\t%f\n", sat->peo);
  fprintf(f, "pgho:\t\t%f\n", sat->pgho);
  fprintf(f, "pho:\t\t%f\n", sat->pho);
  fprintf(f, "pinco:\t\t%f\n", sat->pinco);
  fprintf(f, "plo:\t\t%f\n", sat->plo);
  fprintf(f, "se2:\t\t%f\n", sat->se2);
  fprintf(f, "se3:\t\t%f\n", sat->se3);
  fprintf(f, "sgh2:\t\t%f\n", sat->sgh2);
  fprintf(f, "sgh3:\t\t%f\n", sat->sgh3);
  fprintf(f, "sgh4:\t\t%f\n", sat->sgh4);
  fprintf(f, "sh2:\t\t%f\n", sat->sh2);
  fprintf(f, "sh3:\t\t%f\n", sat->sh3);
  fprintf(f, "si2:\t\t%f\n", sat->si2);
  fprintf(f, "si3:\t\t%f\n", sat->si3);
  fprintf(f, "sl2:\t\t%f\n", sat->sl2);
  fprintf(f, "sl3:\t\t%f\n", sat->sl3);
  fprintf(f, "sl4:\t\t%f\n", sat->sl4);
  fprintf(f, "xgh2:\t\t%f\n", sat->xgh2);
  fprintf(f, "xgh3:\t\t%f\n", sat->xgh3);
  fprintf(f, "xgh4:\t\t%f\n", sat->xgh4);
  fprintf(f, "xh2:\t\t%f\n", sat->xh2);
  fprintf(f, "xh3:\t\t%f\n", sat->xh3);
  fprintf(f, "xi2:\t\t%f\n", sat->xi2);
  fprintf(f, "xi3:\t\t%f\n", sat->xi3);
  fprintf(f, "xl2:\t\t%f\n", sat->xl2);
  fprintf(f, "xl3:\t\t%f\n", sat->xl3);
  fprintf(f, "xl4:\t\t%f\n", sat->xl4);
  fprintf(f, "zmol:\t\t%f\n", sat->zmol);
  fprintf(f, "zmos:\t\t%f\n", sat->zmos);
  fprintf(f, "-- Flags\n");
  fprintf(f, "Deepsp:\t\t%c\n", sat->method);
  fprintf(f, "Loworb:\t\t%d\n", sat->isimp);
  fclose(f);
}
