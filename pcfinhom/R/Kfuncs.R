###############################################################################
# K_global from the paper
# argument list is based on spatstat Kinhom, though many of those arguments
# are not relevant here. ... are passed to density.ppp or densityfun.ppp
# if no lambda, use an analytical kernel version of expectedCrossPairs
Kglobal <-
function(X, lambda=NULL, ..., sigma=bw.CvL(X), r=NULL, rmax=NULL, breaks=NULL,
            analytical=NULL, normtol=.001, discrete.lambda=FALSE,
            interpolate=FALSE, interpolate.fac=10, isotropic=FALSE,
            leaveoneout=FALSE) {
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

    # do the analytical method if no lambda is provided
    if (is.null(analytical)) {
        analytical <- is.null(lambda)
    }

    # Set up lambda if applicable
    if (!analytical) {
        if (!is.null(sig_tmp <- attr(lambda, "sigma"))) {
            sigma <- sig_tmp
        }
        lambda <- fixLambda(lambda, X, discrete.lambda, sigma, ...)
    }

    pairs <- closepairs(X, rmax, twice=FALSE)
    hx <- pairs$xi - pairs$xj
    hy <- pairs$yi - pairs$yj
    pairdist <- pairs$d

    # Get fs, depending on isotropic and interpolate options
    if (isotropic) {
        rh <- sqrt(hx^2 + hy^2)
        if (interpolate) {
            dr <- sigma/interpolate.fac
            rcheck <- seq(0, max(rh) + dr, by=dr)
        } else {
            rcheck <- rh
        }

        fcheck <- if (analytical) {
                Y <- if (leaveoneout) NULL else X
                expectedCrossPairs_kernel_iso(X,Y,rcheck, sigma=sigma)
            } else expectedPairs_iso(lambda, rcheck, tol=normtol)
        f <- fcheck / (2 * pi * rcheck)
        if (interpolate) {
            f <- approx(rchecks, f, xout=rh)$y
        }
    } else { # !isotropic
        if (interpolate) {
            dhx <- sigma/interpolate.fac
            npt <- ceil(rmax/dhx)
            xs <- (-npt:npt)*dhx
            lathx <- outer(xs, xs, function(x,y) x)
            lathy <- outer(xs, xs, function(x,y) y)
        } else {
            lathx <- hx
            lathy <- hy
        }

        if (analytical) {
            Y <- if (leaveoneout) NULL else X
            f <- expectedCrossPairs_kernel(X, Y, lathx, lathy, sigma)
        } else {
            f <- expectedPairs(lambda, lathx, lathy, tol=normtol)
        }

        if (interpolate) {
            dim(f) <- c(2*npt + 1, 2*npt + 1)
            latf.im <- as.im(f, xrows=xs, ycols=ys)

            f <- interp.im(latf.im, hx, hy)
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
Kinhomcross <-
function(X, Y, lambdaX=NULL, lambdaY=NULL, ..., sigma=bw.CvL(X), r=NULL,
            rmax=NULL, breaks=NULL, normtol=.001, analytical=NULL,
            discrete.lambda=FALSE, interpolate=FALSE, isotropic=FALSE,
            interpolate.fac=10, leaveoneout=TRUE) {
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
    if (is.null(analytical)) {
        analytical <- is.null(lambdaX) && is.null(lambdaY)
    }

    if (!analytical) {
        if (!is.null(sig_tmp <- attr(lambdaX, "sigma"))) {
            sigma <- sig_tmp
        }

        lambdaX <- fixLambda(lambdaX, X, discrete.lambda, sigma, ...)
        lambdaY <- fixLambda(lambdaY, Y, discrete.lambda, sigma, ...)
    }

    pairs <- crosspairs(X, Y, rmax)
    hx <- pairs$xi - pairs$xj
    hy <- pairs$yi - pairs$yj
    pairdist <- pairs$d

    # Get fs, depending on isotropic and interpolate options
    if (isotropic) {
        rh <- sqrt(hx^2 + hy^2)
        if (interpolate) {
            dr <- sigma/interpolate.fac
            rcheck <- seq(0, max(rh) + dr, by=dr)
        } else {
            rcheck <- rh
        }

        fcheck <- if (analytical) {
                expectedCrossPairs_iso_kernel(X,Y,rcheck, sigma=sigma)
            } else expectedCrossPairs_iso(lambdaX, lambdaY, rcheck, tol=normtol)
        f <- fcheck
        if (interpolate) {
            f <- approx(rcheck, f, xout=rh)$y
        }
    } else { # !isotropic
        if (interpolate) {
            dhx <- sigma/interpolate.fac
            npt <- ceil(rmax/dhx)
            xs <- (-npt:npt)*dhx
            lathx <- outer(xs, xs, function(x,y) x)
            lathy <- outer(xs, xs, function(x,y) y)
        } else {
            lathx <- hx
            lathy <- hy
        }

        if (analytical) {
            f <- expectedCrossPairs_kernel(X, Y, lathx, lathy, sigma)
        } else {
            f <- expectedCrossPairs(lambdaX, lambdaY, lathx, lathy, tol=normtol)
        }

        if (interpolate) {
            dim(f) <- c(2*npt + 1, 2*npt + 1)
            latf.im <- as.im(f, xrows=xs, ycols=ys)

            f <- interp.im(latf.im, hx, hy)
        }
    }

    lambdaXs <- density.ppp(X, ..., sigma, at="points", leaveoneout=leaveoneout)
    lambdaYs <- density.ppp(Y, ..., sigma, at="points", leaveoneout=leaveoneout)
    lambda2s <- lambdaXs[pairs$i]*lambdaYs[pairs$j]

    bins <- .bincode(pairdist, r, include.lowest=TRUE)
    K <- numeric(length(r))
    Kl <- numeric(length(r))

    for (i in 2:length(r)) {
        K[i] <- sum(1/f[bins == i - 1])
        Kl[i] <- sum(1/lambda2s[bins== i-1])
    }

    K <- cumsum(K)
    Kl <- cumsum(Kl)

    Kf <- data.frame(r=r, theo=pi*r^2, global=K, local=Kl)

    out <- fv(Kf, argu="r", ylab=quote(K(r)), valu="global", fmla= . ~ r,
        alim=c(0,rmax), labl=c("r", "%s[Pois](r)", "%s[global](r)", "%s[local](r)"),
        desc=c("distance argument r", "theoretical poison %s",
        "global correction %s", "local correction %s"), fname="K")

    attr(out, "fs") <- f
    attr(out, "prs") <- pairs
    attr(out, "rh") <- rh
    out
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
