###############################################################################
# These are numerical estimates of the quatities f and f_12 in the paper.
# Perhaps a better name would be "expected number of pairs for poisson process(es)".
# Hence, naming these functions expectedPairs and the like

# expectedPairs is just a wrapper for expectedCrossPairs that takes the right args
expectedPairs <- function(rho, hx, hy=NULL, method=c("mc", "lattice"), tol=.005,
                            dx=diff(as.owin(rho1)$xrange)/200, maxeval=1e6,
                            maxsamp=5e3) {
    expectedCrossPairs(rho1=rho, rho2=NULL, hx, hy, method, tol, dx,
                        maxsamp=maxsamp, maxeval=maxeval)
}

expectedCrossPairs <- function(rho1, rho2=NULL, hx, hy=NULL, method=c("mc", "lattice"),
                            tol=.005, dx=diff(as.owin(rho1)$xrange)/200,
                            maxeval=1e6, maxsamp=5e3) {

    # Validate arguments
    stopifnot(is.im(rho1) || inherits(rho1, "funxy"))
    if (is.im(rho1)) rho1 <- funxy(function(x,y) interp.im(rho1, x, y), as.owin(rho1))
    cross <- !is.null(rho2)
    if (cross) {
        stopifnot(is.im(rho2) || inherits(rho2, "funxy"))
        if (is.im(rho2)) rho2 <- funxy(function(x,y) interp.im(rho2, x, y), as.owin(rho2))
    } else {
        rho2 <- rho1 # for convenience later
    }
    method <- match.arg(method, c("mc", "lattice"))

    use.mc <- FALSE
    use.lattice <- FALSE
    if (method == "mc") {
        use.mc <- TRUE
        stopifnot(is.numeric(tol) && tol > 0)
    } else if (method == "lattice") {
        use.lattice <- TRUE
        stopifnot(is.numeric(dx) && dx > 0 && dx < as.owin(rho1))
    }

    # Check windows...
    W <- as.owin(rho1)
    if (cross) {
        W2 <- as.owin(rho2)
        stopifnot(W$xrange == W2$xrange && W$yrange == W2$yrange)
    }

    # Get coordinates of h
    xy <- xy.coords(hx, hy)
    allhx <- xy$x
    allhy <- xy$y
    nh <- length(hx)

    # allocate results
    epr <- numeric(length(hx))
    epr2 <- numeric(length(hx))
    sampwts <- numeric(length(hx))

    if (is.rectangle(W)) {
        winwts <- (abs(diff(W$xrange)) - abs(allhx))*(abs(diff(W$yrange)) - abs(allhy))
    } else {
        winwts <- numeric(length(allhx))
        for (i in 1:length(allhx)) {
            W2 <- shift(W, - c(allhx[i], allhy[i]))

            winwts[i] <- overlap.owin(W, W2)
        }
    }

    # state of looping
    valid_h <- which(winwts > 0) # if winwts<=0 then f(h) = 0
    nloop <- 0
    inds <- valid_h # hs that haven't converged
    nh_active <- length(valid_h)
    if (nh_active == 0) return(eps) # none of the hs are valid

    # Main loop
    while (nh_active > 0) {
        weights <- sampwts[inds]
        hx <- allhx[inds]
        hy <- allhy[inds]

        # how many samples to take this iteration?
        n_samp <- max(1, min(floor(maxeval/nh_active), maxsamp))

        # Monte carlo sample of points in W
        U <- runifpoint(n_samp, W)

        # Corresponding intensity(s)
        rho1U <- rho1(U)
        rho2U <- if (cross) rho2(U) else rho1U

        # Locations of second point, for each h
        Uplushx <- outer(U$x, hx, `+`)
        Uplushy <- outer(U$y, hy, `+`)
        Uminushx <- outer(U$x, hx, `-`)
        Uminushy <- outer(U$y, hy, `-`)

        # Are these in the window?
        Uplus_inside <- inside.owin(Uplushx, Uplushy, W)
        Uminus_inside <- inside.owin(Uminushx, Uminushy, W)

        # Where yes, corresponding intensity(s)
        rhoUplus <- matrix(0, nrow=n_samp, ncol=nh_active)
        rhoUminus <- matrix(0, nrow=n_samp, ncol=nh_active)
        rhoUplus[Uplus_inside] <- rho2(Uplushx[Uplus_inside], Uplushy[Uplus_inside])
        rhoUminus[Uminus_inside] <- rho1(Uminushx[Uminus_inside], Uminushy[Uminus_inside])

        # Number of new points for each h?
        weights <- weights + colSums(Uplus_inside + Uminus_inside)
        sampwts[inds] <- weights

        # Update estimates
        epr[inds] <- (epr[inds]
                        + (rho1U %*% rhoUplus)
                        + (rho2U %*% rhoUminus))

        epr2[inds] <- (epr2[inds]
                        + ((rho1U^2) %*% (rhoUplus^2))
                        + ((rho2U^2) %*% (rhoUminus^2)))

        # estimate the standard error of the estimates
        sd_est <- ( sqrt( (weights*epr2[inds] - epr[inds]^2) / (weights - 1))
                            / epr[inds])

        passed <- sd_est < tol
        inds <- inds[!passed]
        nh_active <- length(inds)
        nloop <- nloop + 1


        print(c(nh_active, epr[inds[1]], epr2[inds[1]], min(weights), max(weights), min(sd_est), max(sd_est)),
                digits=2)
    }
    # Final output. sampwts is the number of applicable monte carlo samples
    # winwts is the area of W \cap W_{-h}
    epr[valid_h] <- epr[valid_h] * winwts[valid_h] / sampwts[valid_h]
    epr
}

cross_f_lattice <- function(rho1, rho2, hx, hy=NULL, dx=.01, ...) {

    # Get coordinates of h
    xy <- xy.coords(hx, hy)
    hx <- xy$x
    hy <- xy$y

    val <- numeric(length(hx))

    # get number of points on an appropriate lattice
    stopifnot(is.im(rho2))
    W <- as.owin(rho1)
    W2 <- as.owin(rho2)
    stopifnot(W$xrange == W2$xrange && W$yrange == W2$yrange)

    xrange <- W$xrange
    yrange <- W$yrange
    nx <- ceiling(diff(xrange)/dx)
    ny <- ceiling(diff(yrange)/dx)
    dx <- diff(xrange)/nx
    dy <- diff(yrange)/ny

    #compute points of lattice
    xs <- seq(xrange[1] + dx/2, xrange[2], by=dx);
    ys <- seq(yrange[1] + dy/2, yrange[2], by=dy);
    
    latticex <- as.vector(outer(xs, ys, function(x,y) x))
    latticey <- as.vector(outer(xs, ys, function(x,y) y))

    # and associated rho
    latticez <- interp.im(rho1, latticex, latticey)

    # figure out how many hs per call to rho2:
    ncall <- min(ceiling(length(latticez)* length(hx) / 1e6), length(hx))
    npercall <- ceiling(length(hx)/ncall)

    for (i in 1:ncall) {
        # Which inds to get \rho for this time:
        hinds <- (1 + (i-1)*npercall):min(i*npercall, length(hx))

        allxs <- outer(latticex, hx[hinds], `+`)
        allys <- outer(latticey, hy[hinds], `+`)
    
        allrhos <- interp.im(rho2, as.vector(allxs), as.vector(allys))

        dim(allrhos) <- dim(allxs)

        # foreach xs
        for (j in 1:length(hinds)) {
            increment <- sum(latticez*dx*dy*allrhos[,j], na.rm=TRUE)
            val[hinds[j]] <- val[hinds[j]] + increment
        }
    }

    return(val)
}
