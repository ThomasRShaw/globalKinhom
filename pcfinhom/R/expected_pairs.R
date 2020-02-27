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
    if (is.im(rho1)) {
        rho1.im <- rho1
        rho1 <- funxy(function(x,y) interp.im(rho1.im, x, y), as.owin(rho1.im))
    }
    cross <- !is.null(rho2)
    if (cross) {
        stopifnot(is.im(rho2) || inherits(rho2, "funxy"))
        if (is.im(rho2)) {
            rho2.im <- rho2
            rho2 <- funxy(function(x,y) interp.im(rho2.im, x, y), as.owin(rho2.im))
        }
    } else {
        rho2 <- rho1 # for convenience later
    }
    method <- match.arg(method)

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


        #print(c(nh_active, min(weights), max(weights),
        #            min(sd_est), max(sd_est)), digits=2)
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

expectedPairs_iso <- function(rho, r, tol=.001, maxeval=1e6, maxsamp=5e3) {
    expectedCrossPairs_iso(rho, NULL, r, tol, maxeval, maxsamp)
}

expectedCrossPairs_iso <- function(rho1, rho2=NULL, r,
                    tol=.001, maxeval=1e6, maxsamp=5e3) {

    # Validate arguments
    stopifnot(is.im(rho1) || inherits(rho1, "funxy"))
    if (is.im(rho1)) {
        rho1.im <- rho1
        rho1 <- funxy(function(x,y) interp.im(rho1.im, x,y), as.owin(rho1.im))
    }
    cross <- !is.null(rho2)
    if (cross) {
        stopifnot(is.im(rho2) || inherits(rho2, "funxy"))
        if (is.im(rho2)) {
            rho2.im <- rho2
            rho2 <- funxy(function(x,y) interp.im(rho2.im, x,y), as.owin(rho2.im))
        }
    } else {
        rho2 <- rho1
    }

    stopifnot(is.numeric(tol) && tol > 0)

    # Check windows
    W <- as.owin(rho1)
    if (cross) {
        W2 <- as.owin(rho2)
        stopifnot(W$xrange == W2$xrange && W$yrange== W2$yrange)
    }

    # allocate results
    epr <- numeric(length(r))
    epr2 <- numeric(length(r))
    sampwts <- numeric(length(r))

    # winwts account for the W \cap W_{-h} weighting of the window.
    # They are equal to \gamma_r if \rho_1 = \rho_2 = 1 everywhere
    if (is.rectangle(W)) {
        l <- diff(W$xrange)
        h <- diff(W$yrange)
        winwts <- (2*pi - 4*(r/l + r/h) + 2*r^2/(l*h)) * l*h/(2*pi)
    } else {
        #TODO: make them supported
        stop("non-rectangular windows not yet supported")
    }

    #looping state
    allr <- r
    valid_r <- which(winwts > 0)
    nloop <- 0
    inds <- valid_r
    nr_active <- length(valid_r)
    if (nr_active == 0) return(epr)

    se_from_sums <- function(s1, s2, n) sqrt( (n*s2/s1^2 - 1) / (n-1))

    # MC loop
    ndir <- 1
    while (nr_active > 0) {
        weights <- sampwts[inds]
        r <- allr[inds]

        # how many samples this iter?
        n_U <- max(1, min(floor(maxeval/nr_active), maxsamp))

        U <- runifpoint(n_U, W)
        rho1U <- rho1(U)

        # sample 10 random directions
        for (j in 1:ndir) {
            # get a random direction
            hx <- rnorm(1)
            hy <- rnorm(1)
            rh <- sqrt(hx^2 + hy^2)
            hx <- hx/rh
            hy <- hy/rh

            Uphx <- outer(U$x, hx*r, `+`)
            Uphy <- outer(U$y, hy*r, `+`)

            inside <- inside.owin(Uphx, Uphy, W)

            rhoUp <- matrix(0, nrow=n_U, ncol=nr_active)
            rhoUp[inside] <- rho2(Uphx[inside], Uphy[inside])

            weights <- weights + colSums(inside)

            epr[inds] <- (epr[inds] + (rho1U %*% rhoUp))
            epr2[inds] <- (epr2[inds] + (rho1U)^2 %*% (rhoUp)^2)
        }

        sampwts[inds] <- weights
        sd_est <- sqrt(ndir)*(sqrt( (weights*epr2[inds] - epr[inds]^2) / (weights - 1)) / epr[inds])

        passed <- sd_est < tol
        inds <- inds[!is.na(sd_est) & !passed]
        nr_active <- length(inds)
        nloop <- nloop + 1

        #print(c(nr_active, min(weights), max(weights), min(sd_est), max(sd_est)), digits=2)
    }

    epr[valid_r] <- epr[valid_r] * winwts[valid_r] / sampwts[valid_r]
    epr
}

ep_kernel_nodiggle <- function(X,hx,hy,sigma,tol=.005,excess.only=FALSE) {

#TODO: validate inputs. 
CoV <- rep(Inf,length(hx))

num_samples <- numeric(length(hx))
result <- numeric(length(hx))
result2 <- numeric(length(hx))

winwts <- (1-abs(hx))*(1-abs(hy))
validh <- which(winwts > 0)
inds <- validh

x <- X$x; y <- X$y;
n <- length(x)

allhx <- hx
allhy <- hy

fac <- (2*pi*sigma^2)^(-2)

while (any(CoV > tol)) {
    # MC sample u,v
    u <- runif(1)
    v <- runif(1)

    theseh <- (allhx > -u) & (allhx < 1-u) & (allhy > -v) & (allhy < 1-v) & (CoV > tol)

    if (!any(theseh)) next

    hx <- allhx[theseh]
    hy <- allhy[theseh]

    up <- u + hx
    vp <- v + hy

    all_weights <- diggle_weights2(c(u,up),c(v,vp),sigma)
    w <- all_weights[1]*all_weights[1 + (1:length(hx))]

    iterm <- exp(-((u - x)^2 + (v - y)^2)/(2*sigma^2))
    jterms <- outer(x, up, function(x, up) exp(-(up - x)^2/(2*sigma^2))) *
                outer(y, vp, function(x, up) exp(-(up - x)^2/(2*sigma^2)))

    if (!excess.only) {
        newterms <- numeric(length(hx))
        if (length(hx) == 1) {
            for (i in 1:n) {

                js <- (1:n) != i
                newterms <- newterms + iterm[i]*sum(jterms[js,])
            }
        } else {
            for (i in 1:n) {

                js <- (1:n) != i
                newterms <- newterms + iterm[i]*colSums(jterms[js,])
            }
        }
    } else {
        newterms <- colSums(iterm*jterms)
    }

    result[theseh] <- result[theseh] + newterms * w
    result2[theseh] <- result2[theseh] + (newterms * w)^2

    num_samples[theseh] <- num_samples[theseh] + 1

    sd_est <- ( sqrt( (num_samples[theseh]*result2[theseh] - result[theseh]^2) / (num_samples[theseh] - 1))
                        / result[theseh])

    sd_est[num_samples[theseh] < 10] <- Inf


    #if (min(num_samples) > 1000) break
    #CoV[theseh] <- 1000*tol/num_samples[theseh]
    CoV[theseh] <- sd_est
}

print(c(min(num_samples), mean(num_samples), max(num_samples), sd(num_samples)))

result/num_samples*winwts*fac

 
}

diggle_weights2 <- function(x, y, sigma) { #TODO: support windows that aren't unit square?
#    x <- coords(X)$x
#    y <- coords(X)$y
#    W <- as.owin(X)
#    stopifnot("rectangle" == W$type)
#    xrange <- W$xrange
#    yrange <- W$yrange
    xrange <- c(0,1)
    yrange <- c(0,1)
    1/((pnorm((xrange[2]-x)/sigma) - pnorm((xrange[1] - x)/sigma)) *
                (pnorm((yrange[2] - y)/sigma) - pnorm((yrange[1] - y)/sigma)))
}

expectedCrossPairs_tess <- function(X, rho1, rho2=NULL, hx, hy=NULL,
        method=c("mc", "lattice"), tol=.005, dx=diff(as.owin(rho1)$xrange)/200,
        maxeval=1e6, maxsamp=5e3, minus.excess=FALSE, nwin=10) {

    # Validate arguments
    stopifnot(is.im(rho1) || inherits(rho1, "funxy"))
    if (is.im(rho1)) {
        rho1.im <- rho1
        rho1 <- funxy(function(x,y) interp.im(rho1.im, x, y), as.owin(rho1.im))
    }
    cross <- !is.null(rho2)
    if (cross) {
        stopifnot(is.im(rho2) || inherits(rho2, "funxy"))
        if (is.im(rho2)) {
            rho2.im <- rho2
            rho2 <- funxy(function(x,y) interp.im(rho2.im, x, y), as.owin(rho2.im))
        }
    } else {
        rho2 <- rho1 # for convenience later
    }
    method <- match.arg(method)

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


    # Generate quadrats
    if (is.rectangle(W)) {
        nSW <- nwin^2
        xrange <- W$xrange
        yrange <- W$yrange
        dx <- diff(xrange)/nwin
        dy <- diff(yrange)/nwin

        xl <- seq(xrange[1], xrange[2] - dx, by=dx)
        xr <- xl + dx
        xlr <- rbind(xl, xr)
        yl <- seq(yrange[1], yrange[2] - dy, by=dy)
        yr <- yl + dy
        ylr <- rbind(yl, yr)

        SW <- vector("list", nSW)
        for (i in 1:nwin) for (j in 1:nwin)
            SW[[i + (j-1)*nwin]] <- owin(xlr[,i],ylr[,j])
    } else {
        #DT <- dirichlet(X)
        DT <- quadrats(W, nwin,nwin)
        #SW <- DT$tiles # sampling windows
        nSW <- DT$n
# THis is the slow part!
        SW <- lapply(as.list(1:nSW), function(i) as.owin(DT[i]))
    }

    #winwts[i,j] is the overlap of sample window i with W_{-h}
    winwts <- matrix(0, nrow=nSW, ncol=nh)
    for (j in 1:nh) {
        Wshift <- shift(W, c(-allhx[j], -allhy[j]))
        winwts[,j] <- sapply(SW, overlap.owin, Wshift)
    }

    # Generate MC sample points in each sample window
    nper <- 10
    winx <- sapply(SW, function(W) W$xrange[1])
    winy <- sapply(SW, function(W) W$yrange[1])
    samplex <- outer(winx, runif(nper, 0, dx), `+`)
    sampley <- outer(winy, runif(nper, 0, dy), `+`)

    # allocate results
    epr <- matrix(0, nrow=nSW, ncol=nh)
    epr2 <- matrix(0, nrow=nSW, ncol=nh)
    sampwts <- matrix(0, nrow=nSW, ncol=nh)

    valid_h <- (winwts > 0) # if winwts<=0 then f(h) = 0

    # state of looping
    k <- 0
    thous <- 0
    while (TRUE) {
        k <- k + 1

        # Generate more sample points if necessary
        if (k > ncol(samplex)) {
            samplePP <- Map(runifpoint, SW, n=nper)
            samplex <- outer(winx, runif(nper, 0, dx), `+`)
            sampley <- outer(winy, runif(nper, 0, dy), `+`)
            k <- 1
            thous <- thous + 1
        }

        ux <- samplex[,k]
        uy <- sampley[,k]

        # TODO: stop sampling hs that have converged?
        hx <- allhx
        hy <- allhy

        # Corresponding intensity(s)
        rho1U <- rho1(ux, uy)

        # Locations of second point, for each h
        Uplushx <- outer(ux, hx, `+`)
        Uplushy <- outer(uy, hy, `+`)

        # Are these in the window?
        Uplus_inside <- matrix(FALSE, nrow=nSW, ncol=nh) # so that it's logical
        Uplus_inside[valid_h] <- inside.owin(Uplushx[valid_h], Uplushy[valid_h], W)

        inds <- Uplus_inside

        # Where yes, corresponding intensity(s)
        rhoUplus <- matrix(0, nrow=nSW, ncol=nh)
        rhoUplus[inds] <- rho2(as.vector(Uplushx[inds]), Uplushy[inds])

        # Number of new points for each h/SW?
        sampwts[inds] <- sampwts[inds] + 1

        # Update estimates
        newterm <- (rho1U * rhoUplus)

        if (minus.excess) {
            ww <- diggle_weights2(ux, uy, sigma)
            wwp <- diggle_weights2(Uplushx, Uplushy, sigma)

            s <- matrix(0, nrow=nrow(Uplushx), ncol=ncol(Uplushx))
            for (l in 1:npoints(X)) {
                jterms <- matrix(0, nrow=nrow(Uplushx), ncol=ncol(Uplushx))
                x <- X$x[l]
                y <- X$y[l]
                
                # one for each sample u
                iterm <- exp(-((ux - x)^2 + (uy - y)^2)/(2*sigma^2))
                # each u by each h
                jterms[inds] <- exp(-((Uplushx[inds] - x)^2 +
                                        (Uplushy[inds] - y)^2)/(2*sigma^2))

                s <- s + iterm*jterms
                #browser()
            }
            newterm <- newterm - s*(ww*wwp)/(2*pi*sigma^2)^2
        }
        epr <- (epr + newterm)
        epr2 <- (epr2 + newterm^2)

        # Figure out if we're done
        if (any(apply(sampwts, 2, max) < 2)) next


        missingwinwts <- winwts
        missingwinwts[sampwts > 1] <- 0
        #print(colSums(missingwinwts)/colSums(winwts))
        if (any(colSums(missingwinwts)/colSums(winwts) > tol)) next

        # estimate the standard error of the estimates
        #sd_est <- winwts^2*( (sampwts*epr2 - epr^2) / (sampwts - 1)) / epr^2
        #sd_est[is.nan(sd_est)] <- 0
        var_est <- (1/(sampwts*(sampwts-1)))*(epr2 - epr/sampwts)
        var_est[is.nan(var_est) | is.infinite(var_est)] <- 0
        var_est[(sampwts < 2) & valid_h] <- max(var_est)
        var_est <- var_est*(winwts)^2

        tot_sd_est <- sqrt(colSums(var_est))

        mu_est <- winwts/sampwts * epr
        mu_est[is.nan(mu_est)] <- 0
        tot_mu_est <- colSums(mu_est)

        #browser()
        if (any(tot_mu_est == 0)) next
        if (any(tot_sd_est/tot_mu_est > tol)) next

        print(c(k+nper*thous, min(colSums(sampwts)), min(rowSums(sampwts)), min(sampwts[valid_h])))
        #print(tot_sd_est/tot_mu_est)
        break

    }
    # Final output. sampwts is the number of applicable monte carlo samples
    # winwts is the area of W \cap W_{-h}

    epr[valid_h] <- epr[valid_h] * winwts[valid_h] / sampwts[valid_h]
    epr[sampwts==0] <- 0
    #list(colSums(epr),winwts,sampwts)
    colSums(epr)
}
