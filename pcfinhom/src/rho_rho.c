#include <Rmath.h>
#include <R_ext/Utils.h>

#define TWOPI 6.2831853071795

void rho_rho(nquery, xq, yq, ndata, xd, yd, nsep, xh, yh, rmaxi, sig, result) 
  /* inputs */
  int *nquery;            /* number of locations to be interrogated */
  double *xq, *yq;    /* (x,y) coordinates to be interrogated */
  int *ndata;            /* number of data points */
  double *xd, *yd;    /* (x,y) coordinates of data */
  int *nsep;
  double *xh, *yh;  /* (x,y) coordinates of h (one for now) */
  double *rmaxi;    /* maximum distance at which points contribute */
  double *sig;      /* Gaussian sd */
  /* output */
  double *result;   /* vector of computed density values */
{
  double coef, resulti; 
  double sigma, twosig2; 
  int i,j,k, jleft,kleft, h, hleft;
  double rmax, r2max;
  int nq,nd,nh;
  double xpi, ypi, xqi, yqi, d2, d2p, dx, dy, dxp, dyp;
  double x1left, x2left, hxleft;
  double hminx, hmaxx;
  int hright, jright;

  sigma = *sig;
  twosig2 = 2.0 * sigma * sigma;
  coef = 1.0/(TWOPI * sigma * sigma);
  coef = coef*coef;

  rmax = *rmaxi;
  r2max = rmax*rmax;

  nq = *nquery;
  nd = *ndata;
  nh = *nsep;

  /* since they're sorted */
  hminx = xh[0];
  hmaxx = xh[nh-1];

  if(nq == 0 || nd == 0 || nh == 0)
    return;

  /* jleft[i] <= jleft[i+1], likewise kleft, so only initialize once */
  /* i indexes the query points xq,yq */
  for (i = 0; i < nq; i++) {
    xqi = xq[i];

    /* no hs that suit this xq */
    //if (hmaxx + xqi < 0) continue;
    /* no hs that suit any xq bigger than this one */
    //if (hminx + xqi > 1) break;

    yqi = yq[i];

    /* given xq, range of h and range of x1 are determined. */
    /* index for h is h, index for x1 is j */
    /* hxleft for finding left-most h in the window */
    hxleft = -xqi;
    hleft = 0;

    /* get relevant range of j for this xq */
    x1left = xqi - rmax;
    jleft=0;
    while((xd[jleft] < x1left) && (jleft + 1 < nd)) ++jleft;
    jright = nd;
    while(xd[jright - 1] > xqi + rmax) --jright;

    /* get leftmost h s.t. xq = xpi - xh[h] is in the window (> 0) */
    while ((xh[hleft] < hxleft) && (hleft + 1 < nh)) ++hleft;
    hright = nh;
    while (xh[hright - 1] > (1 - xqi)) --hright;

    for (h = hleft; h < hright; h++) {
      xpi = xqi + xh[h];
      ypi = yqi + yh[h];

      /* for sorted xq, xp is also sorted */
      /* NOTE: this assumes window is unit square */
      if (xpi > 1) break;
      if (ypi > 1 || ypi < 0) continue;

      resulti = 0; /* re-zero the accumulator */

      for (j = jleft; j < jright; j++) {
        dx = xd[j] - xqi;
        dy = yd[j] - yqi;

        d2 = dx*dx + dy*dy;
        if (d2 <= r2max) {
          kleft = 0;

          x2left = xpi - rmax;
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
      result[i*nh + h] = resulti*coef;
    }
  }
}
