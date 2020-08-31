#include <Rmath.h>
#include <R_ext/Utils.h>

#define TWOPI 6.2831853071795


void rho_rho_oneh(nquery, xq, yq, ndata, xd, yd, xh, yh, rmaxi, sig, result, result2) 
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
  double *result2;   /* vector of computed density values */
{
  double coef,resulti,resulti2,contrib;
  double sigma, twosig2; 
  int i,j,k, jleft,kleft;
  double rmax, r2max;
  int nq,nd;
  double xpi, ypi, xqi, yqi, d2, d2p, dx, dy, dxp, dyp;
  double x1left, x2left;
  int jright;

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

  resulti=0;
  resulti2=0;


  /* jleft[i] <= jleft[i+1], likewise kleft, so only initialize once */
  /* i indexes the query points xq,yq */
  for (i = 0; i < nq; i++) {
    xqi = xq[i];
    yqi = yq[i];

    /* given xq, range of h and range of x1 are determined. */
    /* index for h is h, index for x1 is j */
    /* get relevant range of j for this xq */
    x1left = xqi - rmax;
    jleft=0;
    jright = nd;
    while((xd[jleft] < x1left) && (jleft + 1 < nd)) ++jleft;
    while(xd[jright - 1] > xqi + rmax) --jright;

    xpi = xqi + *xh;
    ypi = yqi + *yh;

    for (j = jleft; j < jright; j++) { /* first data loop */
      dx = xd[j] - xqi;
      if (dx > rmax) break;
      dy = yd[j] - yqi;

      d2 = dx*dx + dy*dy;
      if (d2 <= r2max) {

        kleft = 0;
        x2left = xpi - rmax;
        while((xd[kleft] < x2left) && (kleft + 1 < nd)) ++kleft;

        for (k=kleft; k<nd; k++) { /* second data loop */
          if (k == j) continue;

          dxp = xd[k] - xpi;
          if (dxp > rmax) break;
          dyp = yd[k] - ypi;
          d2p = dxp * dxp + dyp * dyp;
          if (d2p < r2max) {
            /* now contribute */
            contrib =  exp(-(d2 + d2p)/twosig2);
            resulti += contrib;
            resulti2 += contrib*contrib;
          }
        } /* end second data loop */
      }
    } /* first data loop */
  }

  *result = resulti*coef;
  *result2 = resulti2*coef*coef;
}
