#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <stdio.h>

SEXP expectedCrossPairs(SEXP hx, SEXP hy, SEXP tol) {
/* rho1 and rho2 are funxy, but missing for now
 * for now assume rho1(x,y) = rho2(x,y) = 100 + 100*x
 * hx and hy are numeric of the same length, and tol is how small
 * you want the standard deviation to be before terminating
 * Oh yeah and the Window is the unit square: [0,1] x [0,1]
 */
    int nh, *n_U;
    double *d_hx, *d_hy, *d_tol, *sd_est;
    double x, xp, y, yp, rho;
    double *eps, *ep2s;
    double *f;
    int reached_tol=FALSE;
    SEXP ans;
    PROTECT(hx = AS_NUMERIC(hx));
    PROTECT(hy = AS_NUMERIC(hy));
    PROTECT(tol = AS_NUMERIC(tol));

    nh = LENGTH(hx);
    if (nh != LENGTH(hy)) error("hx and hy must be of same length");

    d_hx = NUMERIC_POINTER(hx);
    d_hy = NUMERIC_POINTER(hy);
    d_tol = NUMERIC_POINTER(tol);

    PROTECT(ans = allocVector(REALSXP, nh));
    eps = malloc(nh*sizeof(double));
    ep2s = malloc(nh*sizeof(double));
    sd_est = malloc(nh*sizeof(double));
    n_U = malloc(nh*sizeof(int));

    for (int j=0; j<nh; j++) {
        eps[j] = 0;
    }

    while (!reached_tol) {
        for (int i=0; i<1000; i++) {
            x = runif(0, 1);
            y = runif(0, 1);
            rho = x * 100 + 100;

            for (int j=0; j<nh; j++) {
                xp = x + d_hx[j];
                yp = y + d_hy[j];
                if (xp > 0 && xp < 1 && yp > 0 && yp < 1) {
                    eps[j] += (rho * (xp * 100 + 100));
                    ep2s[j] += R_pow_di(rho * (xp * 100 + 100), 2);
                    n_U[j]++;
                }
                
                xp = x - d_hx[j];
                yp = y - d_hy[j];
                if (xp > 0 && xp < 1 && yp > 0 && yp < 1) {
                    eps[j] += (rho * (xp * 100 + 100));
                    ep2s[j] += R_pow_di(rho * (xp * 100 + 100), 2);
                    n_U[j]++;
                }
            }
        }
        reached_tol = TRUE;
        for (int j=0; j<nh; j++) {
            sd_est[j] = sqrt((ep2s[j]*n_U[j] - R_pow_di(eps[j], 2)) / ((double) n_U[j] - 1))
                                        / eps[j];
            reached_tol = reached_tol && (sd_est[j] < *d_tol);
        }
    }


    for (int j=0; j<nh; j++) {
        REAL(ans)[j] = eps[j]/((double) n_U[j]);
        Rprintf("%0.2e\t", sd_est[j]);
    }

    free(eps); free(ep2s); free(sd_est); free(n_U);
    UNPROTECT(4);

    return (ans);

}

