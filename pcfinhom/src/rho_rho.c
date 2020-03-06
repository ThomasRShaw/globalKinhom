#include <Rmath.h>
#include <R_ext/Utils.h>

#define TWOPI 6.2831853071795

void rho_rho(nquery, xq, yq, ndata, xd, yd, xh, yh, rmaxi, sig, result) 
  /* inputs */
  int *nquery;            /* number of locations to be interrogated */
  double *xq, *yq;    /* (x,y) coordinates to be interrogated */
  int *ndata;            /* number of data points */
  double *xd, *yd;    /* (x,y) coordinates of data */
  double *xh, *yh;  /* (x,y) coordinates of h (one for now) */
  double *rmaxi;    /* maximum distance at which points contribute */
  double *sig;      /* Gaussian sd */
  /* output */
  double *result;   /* vector of computed density values */
{
  double coef, resulti; 
  double sigma, twosig2; 
  int i,j,k, jleft,kleft;
  double rmax, r2max;
  int nq,nd;
  double xpi, ypi, xqi, yqi, d2, d2p, dx, dy, dxp, dyp;
  double x1left, x2left;

  sigma = *sig;
  twosig2 = 2.0 * sigma * sigma;
  coef = 1.0/(TWOPI * sigma * sigma);
  coef = coef*coef;

  rmax = *rmaxi;
  r2max = rmax*rmax;

  nq = *nquery;
  nd = *ndata;

  if(nq == 0 || nd == 0)
    return;

  /* jleft[i] <= jleft[i+1], likewise kleft, so only initialize once */
  jleft=0;
  kleft=0;
  /* i indexes the query points xq,yq */
  for (i = 0; i < nq; i++) {
    xqi = xq[i];
    yqi = yq[i];

    xpi = xqi + *xh;
    ypi = yqi + *yh;

    /* for sorted xq, xp is also sorted */
    /* NOTE: this assumes window is unit square */
    if (xpi > 1) return;
    if (ypi > 1 || ypi < 0 || xpi < 0) continue;

    resulti = 0; /* re-zero the accumulator */

    x1left = xqi - rmax;
    while((xd[jleft] < x1left) && (jleft + 1 < nd)) ++jleft;

    for (j = jleft; j < nd; j++) {
      dx = xd[j] - xqi;
      if (dx > rmax) break;
      dy = yd[j] - yqi;

      d2 = dx*dx + dy*dy;
      if (d2 <= r2max) {

        x2left=xpi - rmax;
        while((xd[kleft] < x2left) && (kleft + 1 < nd)) ++kleft;

        for (k=kleft; k<nd; k++) {
          if (k == j) continue;

          dxp = xd[k] - xpi;
          if (dxp > rmax) break;
          dyp = yd[k] - ypi;
          d2p = dxp * dxp + dyp * dyp;
          if (d2p < r2max) {
            /* now contribute */
            resulti += exp(-(d2 + d2p)/twosig2);
          }
        }
      }

    }

    result[i] = resulti*coef;

  }

}
