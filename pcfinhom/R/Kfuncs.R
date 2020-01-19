###############################################################################
# K_global from the paper
# argument list is based on spatstat Kinhom, though many of those arguments
# are not relevant here. ... are passed to density.ppp or densityfun.ppp
# if no lambda, use an analytical kernel version of expectedCrossPairs
Kglobal <-
function(X, lambda=NULL, ..., sigma=bw.CvL(X), r=NULL, rmax=NULL, breaks=NULL,
            normtol=.001, discrete.h=FALSE, isotropic=FALSE, leaveoneout=FALSE) {
    # Check inputs
    verifyclass(X, "ppp")
    W <- as.owin(X)
    npts <- npoints(X)
    AreaW <- area(W)

    # What should the rs be?
    rfixed <- !missing(r) || !missing(breaks)
    rmaxdefault <- if (!is.null(rmax)) rmax else rmax.rule("K", W, npts/AreaW)
    breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max
#     rmax <- .25
#     if (missing(r)) r <- seq(.001, rmax, by=.001)
#     else rmax <- max(r)

    if (missing(lambda)) {
        analytical <- TRUE
    } else {
        analytical <- FALSE

        Wl <- as.owin(lambda)
        stopifnot(W$xrange == Wl$xrange && W$yrange == Wl$yrange)
    }

    pairs <- closepairs(X, rmax, twice=FALSE)
    hx <- pairs$xi - pairs$xj
    hy <- pairs$yi - pairs$yj
    pairdist <- pairs$d

    if (isotropic) {
        rh <- sqrt(hx^2 + hy^2)
        if (discrete.h) {
            rchecks <- seq(0, max(rh), length.out=100)
            fcheck <- expectedPairs_iso(lambda, rchecks, tol=normtol) / (2 * pi * rchecks)
            f <- approx(rchecks, fcheck, xout=rh)$y
        } else {
            f <- expectedPairs_iso(lambda, rh, tol=normtol) / (2 * pi * rh)
        }
    } else {
        if (discrete.h) {
            dhx <- (r[2] - r[1])/2
            npt <- ceil(rmax/dhx)
            xs <- (-npt:npt)*dhx
            lathx <- outer(xs, xs, function(x,y) x)
            lathy <- outer(xs, xs, function(x,y) y)
            if (analytical) latf <- expectedPairs_kernel(X, lathx, lathy, sigma)
            else latf <- expectedPairs(lambda, lathx, lathy, tol=normtol)

            dim(latf) <- c(2*npt + 1, 2*npt + 1)
            latf.im <- as.im(latf, xrows=xs, ycols=ys)

            f <- interp.im(latf.im, hx, hy)
        } else {
            if (analytical) {
                f <- if (leaveoneout) expectedPairs_kernel(X, hx, hy, sigma)
                    else expectedCrossPairs_kernel(X,X,hx, hy, sigma)
            } else f <- expectedPairs(lambda, hx, hy, tol=normtol)
        }
    }

    bins <- .bincode(pairdist, r, include.lowest=TRUE)
    K <- numeric(length(r))

    for (i in 2:length(r)) {
        K[i] <- sum(1/f[bins == i - 1])
    }
    K <- cumsum(K)*2 # make up for twice=FALSE above
    Kf <- data.frame(r=r, theo=pi*r^2, global=K)

    fv(Kf, argu="r", ylab=quote(K(r)), valu="global", fmla= . ~ r,
        alim=c(0,rmax), labl=c("r", "%s[Pois](r)", "%s[global](r)"),
        desc=c("distance argument r", "theoretical poison %s",
        "global correction %s"), fname="K")
}

# If lambdaX and lambdaY are both NULL, use the analytical kernel
# otherwise, use Monte Carlo estimates for f
Kglobalcross <-
function(X, Y, lambdaX=NULL, lambdaY=NULL, ..., sigma=bw.CvL(X), r=NULL,
            rmax=NULL, breaks=NULL, normtol=.001, discrete.lambda=FALSE,
            discrete.h=FALSE, isotropic=FALSE) {
    # Check inputs
    verifyclass(X, "ppp")
    verifyclass(Y, "ppp")
    W <- as.owin(X)
    W2 <- as.owin(Y)
    stopifnot(W$xrange == W2$xrange && W$yrange == W2$yrange)

    npts <- npoints(X)
    AreaW <- area(W)

    # What should the rs be?
    rfixed <- !missing(r) || !missing(breaks)
    rmaxdefault <- if (!is.null(rmax)) rmax else rmax.rule("K", W, npts/AreaW)
    breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max

    # How to compute f?
    # if both lambdaX and lambdaY are missing, use an analytical kernel-based
    # estimator
    if (is.null(lambdaX) && is.null(lambdaY)) {
        analytical <- TRUE
    } else {
        analytical <- FALSE

        lambdaX <- fixLambda(lambdaX, X, discrete.lambda, sigma, ...)
        lambdaY <- fixLambda(lambdaY, Y, discrete.lambda, sigma, ...)
    }

    pairs <- crosspairs(X, Y, rmax)
    hx <- pairs$xi - pairs$xj
    hy <- pairs$yi - pairs$yj
    pairdist <- pairs$d

    if (isotropic) {
        rh <- sqrt(hx^2 + hy^2)
        if (discrete.h) {
            rchecks <- seq(0, max(rh), length.out=100)
            fcheck <- expectedCrossPairs_iso(lambdaX, lambdaY, rchecks,
                                            tol=normtol) / (2 * pi * rchecks)
            f <- approx(rchecks, fcheck, xout=rh)$y
        } else {
            f <- expectedCrossPairs_iso(lambdaX, lambdaY, rh, tol=normtol) /
                                                                (2 * pi * rh)
        }
    } else {
        if (discrete.h) {
            dhx <- (r[2] - r[1])/2
            npt <- ceil(rmax/dhx)
            xs <- (-npt:npt)*dhx
            lathx <- outer(xs, xs, function(x,y) x)
            lathy <- outer(xs, xs, function(x,y) y)
            latf <- expectedCrossPairs(lambdaX, lambdaY, lathx, lathy, tol=normtol)

            dim(latf) <- c(2*npt + 1, 2*npt + 1)
            latf.im <- as.im(latf, xrows=xs, ycols=ys)

            f <- interp.im(latf.im, hx, hy)
        } else {
            if (analytical) {
                f <- expectedCrossPairs_kernel(X,Y,hx, hy, sigma)
            } else {
                f <- expectedCrossPairs(lambdaX, lambdaY, hx, hy, tol=normtol)
            }
        }
    }

    bins <- .bincode(pairdist, r, include.lowest=TRUE)
    K <- numeric(length(r))

    for (i in 2:length(r)) {
        K[i] <- sum(1/f[bins == i - 1])
    }
    K <- cumsum(K)

    Kf <- data.frame(r=r, theo=pi*r^2, global=K)

    fv(Kf, argu="r", ylab=quote(K(r)), valu="global", fmla= . ~ r,
        alim=c(0,rmax), labl=c("r", "%s[Pois](r)", "%s[global](r)"),
        desc=c("distance argument r", "theoretical poison %s",
        "global correction %s"), fname="K")
}


fixLambda <- function(lambdaX, X, discrete.lambda, sigma, ...) {
    W <- as.owin(X)
    if (is.null(lambdaX)) {
        if (discrete.lambda) {
            lambdaX.im <- density(X, sigma, ...)
            lambdaX <- funxy(function(x,y) interp.im(lambdaX.im, x,y), W)
        } else {
            lambdaX <- densityfun(X, sigma, ...)
        }
    } else {
        Wl <- as.owin(lambdaX)
        stopifnot(Wl$xrange == W$xrange && Wl$yrange == W$yrange)
    }

    lambdaX
}
