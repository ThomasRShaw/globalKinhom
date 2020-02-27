# isotropic estimate for pcf. bw is the bandwidth for pcf estimation, _NOT_ for
# intensity estimation. ... are passed to densityfun.ppp()
pcfinhom <-
function(X, lambda=NULL, ..., sigma=bw.CvL(X), r=NULL, rmax=NULL,
    kernel="epanechnikov", bw=NULL, stoyan=0.15, normtol=.001, ratio=FALSE,
    discrete.lambda=FALSE, divisor=c("r", "d"), analytical=NULL,
    leaveoneout=TRUE, interpolate=TRUE, interpolate.fac=10, exp_prs=NULL,
    interpolate.maxdx=diameter(as.owin(X))/100, dump=FALSE) {

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
        analytical <- is.null(lambda) && is.null(exp_prs) # do the analytical method if no lambda is provided
    }

    if (interpolate) {
        dr <- min(sigma/ interpolate.fac, interpolate.maxdx)
        if (dr < r[2])
            interpolate=FALSE
    }
    if (interpolate) {
        r_test <- seq(0, rmax - r[2] + dr, by=dr)
    } else {
        r_test <- r
    }

    if (analytical) {
        Y <- if (leaveoneout) NULL else X
        gammas <- expectedCrossPairs_iso_kernel(X, Y, r_test, sigma=sigma)
    } else if (is.null(exp_prs)) {
        ll <- fixLambda(lambda, X, discrete.lambda, sigma, ...)
        gammas <- expectedPairs_iso(ll, r_test, tol=normtol)
    } else {
        if (is.function(exp_prs)) {
            gammas <- exp_prs(r_test)
        } else {
            stop("exp_prs is unknown format")
        }
    }

    if (interpolate) {
        gammas <- approx(r_test, gammas, r, rule=2)$y
#         spl <- smooth.spline(r_test, gammas, df=length(r_test))
#         gammas <- predict(spl, r)$y
    }

    prs <- closepairs(X, rmax + hmax, what='ijd')

    df <- data.frame(r=r, theo=rep.int(1, length(r)))
    out <- ratfv(df, NULL, gammas, "r", quote(g(r)), "theo",
                NULL, alim, c("r", "%s[Pois](r)"), c("distance argument r",
                "theoretical Poisson %s"), fname="g", ratio=ratio)

    bw.used <- NULL

    kdenG <- sewpcf(prs$d, 1, denargs, 1, divisor=divisor)
    gG <- kdenG$g/gammas
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
    } else {
        lambdas <- lambda(X)
    }
    lambda2 <- lambdas[prs$i]*lambdas[prs$j]
    edgewt <- edge.Trans(W=W, dx=X$x[prs$i] - X$x[prs$j], dy=X$y[prs$i]-X$y[prs$j],
                            paired=TRUE)
    wIJ <- edgewt/lambda2
    kdenL <- sewpcf(prs$d, wIJ, denargs, areaW, divisor="r")
    gL <- kdenL$g
    if (!ratio) {
        out <- bind.fv(out, data.frame(local=gL), "hat(%s)[local](r)",
                        "Local intensity reweighted estimate of %s", "local")
    } else {
        out <- bind.ratfv(out, data.frame(local=gL * lambda2/edgewt), lambda2/edgewt,
                        "hat(%s)[local](r)",
                        "Local intensity reweighted estimate of %s", "local")
    }

    formula(out) <- . ~ r
    fvnames(out, ".") <- setdiff(rev(colnames(out)), c("r", "v"))
    unitname(out) <- unitname(X)
    if (ratio)
        out <- conform.ratfv(out)
    attr(out, "bw") <- bw.used

    if (dump) {
        attr(out, "prs") <- prs
        attr(out, "wIJ") <- prs
        attr(out, "fs") <- gammas
    }

    out
}

pcfcrossinhom <- function(X,Y, lambdaX=NULL, lambdaY=NULL, ...,
    sigma=bw.CvL(X), r=NULL, rmax=NULL, kernel="epanechnikov", bw=NULL,
    stoyan=0.15, normtol=.001, ratio=FALSE, discrete.lambda=FALSE,
    divisor=c("r", "d"), analytical=NULL, interpolate=TRUE,
    interpolate.fac=10, leaveoneout=TRUE, exp_prs=NULL,
    interpolate.maxdx=diameter(as.owin(X))/100, dump=FALSE) {

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

    alim <- c(0, min(rmax, rmaxdefault))

    kernel <- match.kernel(kernel)

    # Bandwidth for the pcf estimate
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
        analytical <- is.null(lambdaX) && is.null(lambdaY) && is.null(exp_prs)
    }

    if (interpolate) {
        dr <- min(sigma/ interpolate.fac, interpolate.maxdx)
        if (dr < r[2])
            interpolate=FALSE
    }
    if (interpolate) {
        r_test <- seq(0, rmax - r[2] + dr, by=dr)
    } else {
        r_test <- r
    }

    if (analytical) {
        gammas <- expectedCrossPairs_iso_kernel(X,Y, r_test, sigma=sigma)
    } else if (is.null(exp_prs)){
        lX <- fixLambda(lambdaX, X, discrete.lambda, sigma, ...)
        lY <- fixLambda(lambdaY, Y, discrete.lambda, sigma, ...)

        gammas <- expectedCrossPairs_iso(lX, lY, r_test, tol=normtol)
    } else {
        if (is.function(exp_prs)) {
            gammas <- exp_prs(r_test)
        } else {
            stop("exp_prs given in unknown form")
        }
    }

    if (interpolate) {
        gammas <- approx(r_test, gammas, r, rule=2)$y
#         spl <- smooth.spline(r_test, gammas, df=length(r_test))
#         gammas <- predict(spl, r)$y
    }

    prs <- crosspairs(X, Y, rmax + hmax, what='all')

    df <- data.frame(r=r, theo=rep(1, length(r)))
    out <- ratfv(df, NULL, gammas, "r", quote(c(r)), "theo", NULL, alim,
                c("r", "%s[Pois](r)"),
                c("distance argument r", "theoretical poisson %s"),
                fname="c", ratio=ratio)

    bw.used <- NULL

    kdenG <- sewpcf(prs$d, 1, denargs, 1, divisor=divisor)
    cG <- kdenG$g/gammas
    bw.used <- attr(kdenG, "bw")
    if (!ratio) {
        out <- bind.fv(out, data.frame(global=cG), "hat(%s)[global](r)",
                "Global intensity reweighted estimate of %s", "global")
    } else {
        out <- bind.ratfv(out, data.frame(global=cG*gammas), gammas,
                "hat(%s)[global](r)", "Global intensity reweighted estimate of %s",
                "global")
    }

    # Do the local version
    if (is.null(lambdaX)) {
        lX <- density.ppp(X, at="points", ..., sigma=sigma, leaveoneout=leaveoneout)
    } else {
        lX <- lambdaX(X)
    }
    if (is.null(lambdaY)) {
        lY <- density.ppp(Y, at="points", ..., sigma=sigma, leaveoneout=leaveoneout)
    } else {
        lY <- lambdaY(Y)
    }
    lambda2 <- lX[prs$i]*lY[prs$j]
    edgewt <- edge.Trans(W=W, dx=prs$dx, dy=prs$dy, paired=TRUE)
    wIJ <- edgewt/lambda2

    kdenL <- sewpcf(prs$d, wIJ, denargs, areaW, divisor=divisor)
    cL <- kdenL$g
    if (!ratio) {
        out <- bind.fv(out, data.frame(local=cL), "hat(%s)[local](r)",
                "Local intensity reweighted estimate of %s", "local")
    } else {
        out <- bind.ratfv(out, data.frame(local=cL*lambda2/edgewt), lambda2/edgewt,
                "hat(%s)[local](r)",
                "Local intensity reweighted estimate of %s", "local")
    }

    formula(out) <- . ~ r
    fvnames(out, ".") <- setdiff(rev(colnames(out)), c("r", "v"))
    unitname(out) <- unitname(X)
    if (ratio)
        out <- conform.ratfv(out)

    attr(out, "bw") <- bw.used
    if (dump) {
        attr(out, "wIJ") <- wIJ
        attr(out, "prs") <- prs
        attr(out, "fs") <- gammas
    }

    out
}
