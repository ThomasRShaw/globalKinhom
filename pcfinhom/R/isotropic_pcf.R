# isotropic estimate for pcf. bw is the bandwidth for pcf estimation, _NOT_ for
# intensity estimation. ... are passed to densityfun.ppp()
global_pcf_iso <-
function(X, lambda=NULL, ..., sigma=bw.CvL(X), r=NULL, rmax=NULL,
            kernel="epanechnikov",
            bw=NULL, stoyan=0.15, normtol=.001, ratio=FALSE,
            discrete.lambda=FALSE, divisor=c("r", "d"), analytical=NULL,
            leaveoneout=TRUE, interpolate=TRUE, interpolate.fac=10) {
    verifyclass(X, "ppp")
    W <- as.owin(X)
    areaW <- area(W)
    npts <- npoints(X)

    rmaxdefault <- if (is.null(rmax)) rmax <- rmax.rule("K", W, npts/areaW)
    breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmax)
    r <- breaks$r
    rmax <- breaks$max

    alim <- c(0, min(rmax, rmaxdefault))

    kernel <- match.kernel(kernel)

    if (is.null(bw) && (kernel == "epanechnikov")) {
        h <- stoyan/sqrt(npts/areaW)
        hmax <- h
        bw <- h/sqrt(5)
    }
    else if (is.numeric(bw)) {
        hmax <- 3 * bw
    }
    else {
        hmax <- 2 * stoyan/sqrt(npts/areaW)
    }
    denargs <- list(kernel=kernel, bw=bw, n=length(r), from=0, to=rmax)

    if (is.null(analytical)) {
        analytical <- is.null(lambda) # do the analytical method if no lambda is provided
    }

    if (interpolate) {
        dr <- sigma / interpolate.fac
        if (dr < r[2]) {
            interpolate=FALSE
            r_test <- r
        } else {
            r_test <- seq(r[2], rmax - r[2] + dr, by=dr)
        }
    }

    if (analytical) {
        Y <- if (leaveoneout) NULL else X
        gammas <- expectedCrossPairs_iso_kernel(X, Y, r_test, sigma=sigma)
    } else {
        lambda <- fixLambda(lambda, X, discrete.lambda, sigma, ...)
        gammas <- expectedPairs_iso(lambda, r_test, tol=normtol)
    }

    prs <- closepairs(X, rmax + 2*bw, what='ijd')

    if (interpolate) {
        spl <- smooth.spline(r_test, gammas/r_test, df=length(r_test))
        gammas <- predict(spl, r)$y * r
    }

    df <- data.frame(r=r, theo=rep.int(1, length(r)))
    out <- ratfv(df, NULL, gammas, "r", quote(g(r)), "theo",
                NULL, alim, c("r", "%s[Pois](r)"), c("distance argument r",
                "theoretical Poisson %s"), fname="g", ratio=ratio)

    bw.used <- NULL

    kdenG <- sewpcf(prs$d, 1, denargs, 1, divisor=divisor)
    gG <- kdenG$g*2*pi*r/gammas
    bw.used <- attr(kdenG, "bw")
    if (!ratio) {
        out <- bind.fv(out, data.frame(global=gG), "hat(%s)[global](r)",
                        "Global intensity reweighted estimate of %s", "global")
    } else {
        out <- bind.ratfv(out, data.frame(global=gG * gammas), gammas,
                        "hat(%s)[global](r)",
                        "Global intensity reweighted estimate of %s", "global")
    }

    if (is.null(lambda)) {
        lambdas <- density.ppp(X, at="points", ..., sigma=sigma, leaveoneout=leaveoneout)
        #lambdas_noloo <- density.ppp(X, at="points", ..., sigma=sigma, leaveoneout=FALSE)
    } else {
        lambdas <- lambda(X)
    }
    lambda2 <- lambdas[prs$i]*lambdas[prs$j]

    edgewt <- edge.Trans(W=W, dx=X$x[prs$i] - X$x[prs$j], dy=X$y[prs$i]-X$y[prs$j],
                            paired=TRUE)
    kdenL <- sewpcf(prs$d, edgewt/lambda2, denargs, areaW, divisor="r")
    gL <- kdenL$g
    if (!ratio) {
        out <- bind.fv(out, data.frame(local=gL), "hat(%s)[local](r)",
                        "Local intensity reweighted estimate of %s", "local")
    } else {
        out <- bind.ratfv(out, data.frame(global=gL * lambda2/edgewt), lambda2/edgewt,
                        "hat(%s)[local](r)",
                        "Local intensity reweighted estimate of %s", "local")
    }

    formula(out) <- . ~ r
    fvnames(out, ".") <- setdiff(rev(colnames(out)), c("r", "v"))
    unitname(out) <- unitname(X)
    if (ratio)
        out <- conform.ratfv(out)
    attr(out, "bw") <- bw.used

    out
}

#TODO: this only works if bw equals the spacing of the rs!!
global_cross_pcf_iso <- function(X,Y, lambdaX=NULL, lambdaY=NULL, ...,
    r=NULL, rmax=NULL, kernel=c("box", "epanechnikov"), bw=NULL, stoyan=0.15,
    normtol=.001, discrete.lambda=FALSE) {
    verifyclass(X, "ppp")
    verifyclass(Y, "ppp")
    W <- as.owin(X)
    W2 <- as.owin(Y)
    stopifnot(W$xrange == W2$xrange && W$yrange == W2$yrange)
    areaW <- area(W)
    npts <- npoints(X)


    rmaxdefault <- if (is.null(rmax)) rmax <- rmax.rule("K", W, npts/areaW)
    breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmax)
    r <- breaks$r
    rmax <- breaks$max

    denargs <- resolve.defaults(list(kernel=kernel, bw=bw),
                    list(n=length(r), from=0, to=rmax))

    if (is.null(lambdaX)) {
        if (discrete.lambda) {
            lambdaX.im <- density(X, eps=.005, ...)
            lambdaX <- funxy(function(x,y) interp.im(lambdaX.im,x,y), W)
        } else {
            lambdaX <- densityfun(X, ...)
        }
    } else {
        Wl <- as.owin(lambdaX)
        stopifnot(W$xrange == Wl$xrange && W$yrange == Wl$yrange)
    }
    if (is.null(lambdaY)) {
        if (discrete.lambda) {
            lambdaY.im <- density(Y, eps=.005, ...)
            lambdaY <- funxy(function(x,y) interp.im(lambdaY.im,x,y), W)
        } else {
            lambdaY <- densityfun(Y, ...)
        }
    } else {
        Wl <- as.owin(lambdaY)
        stopifnot(W$xrange == Wl$xrange && W$yrange == Wl$yrange)
    }

    gammas <- expectedCrossPairs_iso(lambdaX, lambdaY, r, tol=normtol)


    if (is.null(bw)) {
        bw <- diff(r)
        bw <- c(bw, bw[length(bw)])
        breaks <- c(0, r + bw/2)
        breaksl <- 1:length(r)
        breaksr <- breaksl + 1
    } else {
        breaks <- sort(unique(c(r - bw/2, r + bw/2)))
        breaksl <- numeric(length(r))
        breaksr <- numeric(length(r))
        for (i in 1:length(r)) {
            breaksl[i] <- which(breaks == r[i] - bw/2)
            breaksr[i] <- which(breaks == r[i] + bw/2)
        }
    }

    prs <- crosspairs(X, Y, rmax + bw[length(bw)]/2, what='ijd')

    bins <- .bincode(prs$d, breaks, include.lowest=TRUE)

    c <- numeric(length(r))
    for (i in 1:length(r)) {
        if (length(bw) == 1) thisbw <- bw else thisbw <- bw[i]
        c[i] <- sum(bins < breaksr[i] & bins >= breaksl[i]) / thisbw / gammas[i]
    }

    df <- data.frame(r=r, theo=rep(1, length(r)), global=c)

    out <- fv(df, argu="r", ylab=quote(c(r)), "global", fmla= . ~ r,
                alim=c(0,rmax), labl=c("r", "%s[pois](r)", "%s[global](r)"),
                desc=c("distance argument r", "theoretical poisson %s",
                "global correction %s"),
                fname="c")

    out
}
