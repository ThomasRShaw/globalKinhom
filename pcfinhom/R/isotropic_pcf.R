# isotropic estimate for pcf. bw is the bandwidth for pcf estimation, _NOT_ for
# intensity estimation. ... are passed to densityfun.ppp()
global_pcf_iso <-
function(X, lambda=NULL, ..., r=NULL, rmax=NULL,
            kernel=c("box", "epanechnikov"), bw=NULL, stoyan=0.15, normtol=.001,
            discrete.lambda=FALSE) {
    verifyclass(X, "ppp")
    W <- as.owin(X)
    areaW <- area(W)
    npts <- npoints(X)

    rmaxdefault <- if (is.null(rmax)) rmax <- rmax.rule("K", W, npts/areaW)
    breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmax)
    r <- breaks$r
    rmax <- breaks$max

    if (is.null(lambda)) {
        if (discrete.lambda) {
            lambda.im <- density(X, eps=.005, ...)
            lambda <- funxy(function(x,y) interp.im(lambda.im,x,y), W)
        } else {
            lambda <- densityfun(X, ...)
        }
    } else {
        W2 <- as.owin(lambda)
        stopifnot(W$xrange == W2$xrange && W$yrange == W2$yrange)
    }

    gammas <- expectedPairs_iso(lambda, r, tol=normtol)

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

    #maxr <- breaks[length(breaks)]
    prs <- closepairs(X, rmax, what='ijd')

    #TODO: for consistency, this could be done with density.default, kernel, bw
    bins <- .bincode(prs$d, breaks, include.lowest=TRUE)

    g <- numeric(length(r))
    for (i in 1:length(r)) {
        if (length(bw) == 1) thisbw <- bw else thisbw <- bw[i]
        g[i] <- sum(bins < breaksr[i] & bins >= breaksl[i]) / thisbw / gammas[i]
    }
    df <- data.frame(r=r, theo=rep(1, length(r)), global=g)

    out <- fv(df, argu="r", ylab=quote(g(r)), "global", fmla= . ~ r,
                alim=c(0,rmax), labl=c("r", "%s[Pois](r)", "%s[global](r)"),
                desc=c("distance argument r", "theoretical poisson %s",
                "global correction %s"),
                fname="g")

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

    prs <- crosspairs(X, Y, rmax, what='ijd')

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
