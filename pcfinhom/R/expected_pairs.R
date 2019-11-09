###############################################################################
# These are numerical estimates of the quatities f and f_12 in the paper.
# Perhaps a better name would be "expected number of pairs for poisson process(es)".
# Hence, naming these functions expectedPairs and the like
cross_f <- function(rho1, rho2, hx, hy=NULL, dx=.01, ...) UseMethod("cross_f")

f_est <- function(rho, hx, hy=NULL, dx=.01) cross_f(rho, rho, hx, hy, dx)

cross_f.im <- function(rho1, rho2, hx, hy=NULL, dx=.01, ...) {

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

cross_f.funxy <- function(rho1, rho2, hx, hy=NULL, dx, ...) {

    # Get coordinates of h
    xy <- xy.coords(hx, hy)
    hx <- xy$x
    hy <- xy$y

    val <- numeric(length(hx))

    # Check windows
    stopifnot(inherits(rho2, "funxy"))
    W <- as.owin(rho1)
    W2 <- as.owin(rho2)
    stopifnot(W$xrange == W2$xrange && W$yrange == W2$yrange)

    # get number of points on an appropriate lattice
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

    # remove points not in window
    inside <- inside.owin(latticex, latticey, W)
    latticex <- latticex[inside]
    latticey <- latticey[inside]

    # and associated rho
    latticez <- rho1(latticex, latticey)

    # figure out how many hs per call to rho2:
    ncall <- min(ceiling(length(latticez)* length(hx) / 1e6), length(hx))
    npercall <- ceiling(length(hx)/ncall)

    for (i in 1:ncall) {
        # Which inds to get \rho for this time:
        hinds <- (1 + (i-1)*npercall):min(i*npercall, length(hx))

        # and rho(x + h)
        allxs <- outer(latticex, hx[hinds], `+`)
        allys <- outer(latticey, hy[hinds], `+`)
        sz <- dim(allxs)
        allxs <- as.vector(allxs)
        allys <- as.vector(allys)

        inside <- inside.owin(allxs, allys, W)
        allrhos <- numeric(length(allxs))
        allrhos[!inside] <- NA
        
        # densityfun warns when any points are outside of W
        tmprhos <- suppressWarnings( rho2(as.vector(allxs), as.vector(allys)) )

        if (length(tmprhos) == sum(inside)) {
            # densityfun only returns values if they're inside W
            allrhos[inside] <- tmprhos
        } else if (length(tmprhos) == length(allxs)) {
            # manually constructed funxy may return garbage outside of W
            allrhos[inside] <- tmprhos[inside]
        } else {
            stop('rho2 returned unexpected number of values')
        }

        dim(allrhos) <- sz

        # foreach xs
        for (j in 1:length(hinds)) {
            increment <- dx*dy*sum(latticez*allrhos[,j],na.rm=TRUE)
            val[hinds[j]] <- increment
        }
    }
    return(val)
}

# Another version of f, using monte-carlo
cross_f_mc.funxy <- function(rho1, rho2, hx, hy=NULL, tol=.005) {

    # get the coordinates
    xy <- xy.coords(hx, hy)
    hx <- xy$x
    hy <- xy$y
    o.out <- order(hx^2 + hy^2)
    hx <- hx[o.out]
    hy <- hy[o.out]

    val <- numeric(length(hx))
    val2 <- numeric(length(hx))
    allwts <- numeric(length(hx))

    # Check windows
    W <- as.owin(rho1)
    W2 <- as.owin(rho2)
    stopifnot(W$xrange == W2$xrange && W$yrange == W2$yrange)

    # Monte-carlo scheme
    nh <- length(hx)
    nloop <- 0
    discrep <- Inf
    firsth <- 1
    inds <- firsth:nh

    while (TRUE) {
        weights <- allwts[inds]
        nh_active <- length(inds)

        # so you never have too many fevals per iteration
        n_per_iter <- min(5000, ceiling(1e6/nh_active))
        
        U <- runifpoint(n_per_iter, W)

        rho1U <- rho1(U)
        rho2U <- rho2(U)

        thesehx <- hx[inds]
        thesehy <- hy[inds]

        # add and subtract h from the random samples
        Uphx <- outer(U$x, thesehx, `+`)
        Uphy <- outer(U$y, thesehy, `+`)
        Umhx <- outer(U$x, thesehx, `-`)
        Umhy <- outer(U$y, thesehy, `-`)

        insidep <- inside.owin(Uphx, Uphy, W)
        insidem <- inside.owin(Umhx, Umhy, W)
        rho2Uph <- matrix(0, nrow=n_per_iter, ncol=nh_active)
        rho1Umh <- matrix(0, nrow=n_per_iter, ncol=nh_active)

        #sequential version (one line)
        rho2Uph[insidep] <- rho2(Uphx[insidep], Uphy[insidep])
        rho1Umh[insidem] <- rho1(Umhx[insidem], Umhy[insidem])
        
        nperhp <- colSums(insidep)
        nperhm <- colSums(insidem)
        nperh <- nperhp + nperhm

        nperup <- rowSums(insidep)
        nperum <- rowSums(insidem)

        oldval <- val[inds]/weights
        val[inds] <-  (val[inds]
                        + (rho1U %*% rho2Uph)
                        + (rho2U %*% rho1Umh))

        val2[inds] <- (val2[inds]
                        + (rho1U)^2 %*% (rho2Uph)^2
                        + (rho2U)^2 %*% (rho1Umh)^2)

        weights <- weights + nperh

        allwts[inds] <- weights

        var_est <- sqrt(val2[inds] - val[inds]^2 /(weights-1)) / val[inds]

        passed <- var_est < tol
        npassed <- match(FALSE, passed) - 1
        nh_active <- sum(!passed)

        nloop <- nloop + 1
        print(c(n_per_iter, nh_active,
                min(nperum), max(nperum),
                min(var_est), max(var_est)), digits=2)

        if (is.na(npassed)) { #all elements passed
            break
        } else {
            firsth <- firsth + npassed
            inds <- inds[!passed] #firsth:nh
        }
    }

    winwts <- numeric(length(hx))
    if (is.rectangle(W)) {
        winwts <- (abs(diff(W$xrange)) - abs(hx))*(abs(diff(W$yrange)) - abs(hy))
    } else {
        for (i in 1:length(hx)) {
            W2 <- shift(W, - c(hx[i], hy[i]))

            winwts[i] <- overlap.owin(W, W2)
        }
    }

    val <- val * winwts / allwts
    val[o.out] <- val #Sort back to original order
    val
}
