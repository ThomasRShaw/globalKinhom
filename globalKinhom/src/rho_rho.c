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
    if (hmaxx + xqi < 0) {
      continue; /* to next xq */
    }
    /* no hs that suit any xq bigger than this one */
    if (hminx + xqi > 1) {
      break; /* We're done */
    }
    yqi = yq[i];

    /* given xq, range of h and range of x1 are determined. */
    /* index for h is h, index for x1 is j */
    /* get relevant range of j for this xq */
    x1left = xqi - rmax;
    jleft=0;
    jright = nd;
    while((xd[jleft] < x1left) && (jleft + 1 < nd)) ++jleft;
    while(xd[jright - 1] > xqi + rmax) --jright;

    /* get leftmost h s.t. xpi = xqi + xh[h] is in the window (> 0) */
    hxleft = -xqi;
    hleft = 0;
    hright = nh;
    while ((xh[hleft] < hxleft) && (hleft + 1 < nh)) ++hleft;
    while (xh[hright - 1] > (1 - xqi)) --hright;

    /* for (h=0; h<nh; h++) result[i*nh + h] = 0.0;  everything should be 0 already */

    for (j = jleft; j < jright; j++) { /* first data loop */
      dx = xd[j] - xqi;
      dy = yd[j] - yqi;

      d2 = dx*dx + dy*dy;
      if (d2 <= r2max) {

        for (h = hleft; h < hright; h++) { /* h loop */
          xpi = xqi + xh[h];
          ypi = yqi + yh[h];

          /* NOTE: this assumes window is unit square */
          //if (xpi > 1) break; /* should never happen */
          if (ypi > 1 || ypi < 0) continue;

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
              result[i*nh + h] += exp(-(d2 + d2p)/twosig2);
            }
          } /* end second data loop */
        } /* end h loop */
      }
    } /* first data loop */
    
    for (h=hleft; h < hright; h++) result[i*nh + h] *= coef;
  }
}

void rho_rho_excess(nquery, xq, yq, ndata, xd, yd, nsep, xh, yh, rmaxi, sig, result) 
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
  int i, j, jleft, h, hleft;
  double rmax, r2max;
  int nq,nd,nh;
  double xpi, ypi, xqi, yqi, d2, d2p, dx, dy, dxp, dyp;
  double x1left, hxleft;

  sigma = *sig;
  twosig2 = 2.0 * sigma * sigma;
  coef = 1.0/(TWOPI * sigma * sigma);
  coef = coef*coef;

  rmax = *rmaxi;
  r2max = rmax*rmax;

  nq = *nquery;
  nd = *ndata;
  nh = *nsep;

  if(nq == 0 || nd == 0 || nh == 0)
    return;

  /* jleft[i] <= jleft[i+1], likewise kleft, so only initialize once */
  /* i indexes the query points xq,yq */
  jleft=0;
  hleft = nh;

  for (i = 0; i < nq; i++) {
    xqi = xq[i];
    yqi = yq[i];

    /* given xq, range of h and range of x1 are determined. */
    /* index for h is h, index for x1 is j */
    /* get relevant range of j for this xq */
    x1left = xqi - rmax;
    while((xd[jleft] < x1left) && (jleft + 1 < nd)) ++jleft;

    /* get leftmost h s.t. xpi = xqi + xh[h] is in the window (> 0) */
    hxleft = -xqi;
    while ((xh[hleft] > hxleft) && (hleft > 0)) --hleft;

    /* TODO: only check ||h|| < rmax? */
    for (h = hleft; h < nh; h++) { /* h loop */
      xpi = xqi + xh[h];
      if (xpi > 1) break;
      ypi = yqi + yh[h];
      /* NOTE: this assumes window is unit square */
      if (ypi > 1 || ypi < 0) continue;

      resulti = 0;

      for (j = jleft; j < nd; j++) { /* first data loop */
        dx = xd[j] - xqi;
        if (dx > rmax) break;
        dy = yd[j] - yqi;

        d2 = dx*dx + dy*dy;
        if (d2 <= r2max) {
          dxp = dx - xh[h];
          if (dxp > rmax) break;
          dyp = dy - yh[h];

          d2p = dxp * dxp + dyp * dyp;
          if (d2p < r2max) {
            /* now contribute */
            resulti += exp(-(d2 + d2p)/twosig2);
          }
        } /* first data loop */
      }
      result[nh*i + h] = resulti*coef;
    } /* end h loop */
  }
}
