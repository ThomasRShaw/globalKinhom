# Here are some tools for doing fourier-transform-based expectedCrossPairs

expectedCrossPairs_fourier <- function(X, Y, hx, hy, bx, by=bx, tol=.001) {
    ftX <- my_pp_ft(X, bx)
    ftY <- my_pp_ft(Y, by)

    const <- 1/(2*pi*(bx^2 + by^2))

    edgewts <- (1 - abs(hx))*(1- abs(hy))

    weights <- 0
    f <- numeric(length(hx))
    varf <- 0

    n <- 1000
    while(TRUE) {
        weights <- weights + n

        ux <- rnorm(n, sd=1/sqrt(bx^2 + by^2))
        uy <- rnorm(n, sd=1/sqrt(by^2 + by^2))
        ux <- ux
        uy <- uy

        u_part <- ftX(ux, uy) * Conj(ftY(ux, uy))

        ft_part <-  exp(0+1i * (outer(ux, hx) + outer(uy,hy)))

        all_vals <- Re(u_part * ft_part) # amounts counting u and -u; all_vals(-u) == all_vals(u)*

        f <- f + colSums(all_vals)

        varf <- varf + apply(all_vals, 2, var) * n

        s <- sqrt(varf)/abs(f)
        #print(c(min(abs(s)), max(abs(s))))
        if (max(abs(s)) < tol) break;

    }

    f*const/weights #*edgewts
}

sinc <- function(x, y) {
    v <- sin(pi * x) / x * sin(pi * y) /y / pi^2
    v[x*y == 0] <- 1
    v
}

my_pp_ft <- function(X, b=bw.CvL(X)) {
    xy <- coords(X); myx <- xy$x; myy <- xy$y;
    W <- as.owin(X)
    xrange <- W$xrange
    yrange <- W$yrange
    # diggle
    w <- 1/((pnorm((xrange[2]-myx)/b) - pnorm((xrange[1] - myx)/b)) *
                (pnorm((yrange[2] - myy)/b) - pnorm((yrange[1] - myy)/b)))
    # non-diggle
#     w <- ??
    # TODO: The reason this is wrong is that it's not edge corrected. i.e. it is
    # \int_{-\infty}^\infty \rho(x) \exp(-ikx) instead of \int_W ...
    # My research on this indicates that the edge corrections amounts to the
    # difference of two erfz s (with argument depending on the point), but I
    # haven't figured out what the arguments should be yet. all in all seems like
    # this method is unlikely to end up being faster.
    function(kx,ky=NULL) {
        kxy <- xy.coords(kx, ky)
        kx <- as.vector(kxy$x)
        ky <- as.vector(kxy$y)
        drop(w %*% exp(0-1i * (outer(myx, kx) + outer(myy, ky))))
    }
}

diggle_weights <- function(X, b) {
    x <- coords(X)$x
    y <- coords(X)$y
    xrange <- c(0,1)
    yrange <- c(0,1)
    1/((pnorm((xrange[2]-x)/b) - pnorm((xrange[1] - x)/b)) *
                (pnorm((yrange[2] - y)/b) - pnorm((yrange[1] - y)/b)))
}

auto_ep_analytical <- function(X, hx, hy, b = bw.CvL(X)) 
    cross_ep_analytical(X, NULL, hx, hy, b)

cross_ep_analytical <- function(X, Y=NULL, hx, hy, b = bw.CvL(X)) {
    autof <- is.null(Y)
    wX <- diggle_weights(X, b)
    wY <- if (autof) wX else diggle_weights(Y, b)
    if (autof) Y <- X

    res <- numeric(length(hx))

    doublecount <- if (autof) 2 else 1
    cfac <- 1/(4*pi*b^2) * doublecount

    lasti <- if (autof) length(wX) - 1 else length(wX)
    for (i in 1:lasti) {
        xX <- X$x[i]
        yX <- X$y[i]
        firstj <- if (autof) i+1 else 1

        for (j in firstj:length(wY)) {
            xY <- Y$x[j]
            yY <- Y$y[j]

            res <- res + (wX[i] * wY[j]) *
                exp(- (xX - xY + hx)^2/4/b^2) *
                exp(- (yX - yY + hy)^2/4/b^2) *
                (pnorm( sqrt(2)*(1 - (xX + xY + abs(hx))/2)/b) - pnorm(-sqrt(2)/b*(xX + xY - abs(hx))/2) )*
                (pnorm( sqrt(2)*(1 - (yX + yY + abs(hy))/2)/b) - pnorm(-sqrt(2)/b*(yX + yY - abs(hy))/2) )
        }
    }

    cfac * res
}
