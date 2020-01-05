# These are edge-correction weights for Gaussian kernel intensity estimators
# Assumes owin is a rectangle
# result is (\int_W \kappa(x - x') dx')^{-1} for each x in X
diggle_weights <- function(X, sigma) {
    x <- coords(X)$x
    y <- coords(X)$y
    W <- as.owin(X)
    stopifnot("rectangle" == W$type)
    xrange <- W$xrange
    yrange <- W$yrange
    1/((pnorm((xrange[2]-x)/sigma) - pnorm((xrange[1] - x)/sigma)) *
                (pnorm((yrange[2] - y)/sigma) - pnorm((yrange[1] - y)/sigma)))
}

expectedPairs_kernel <- function(X, hx, hy, sigma = bw.CvL(X)) 
    expectedCrossPairs_kernel(X, NULL, hx, hy, sigma)

expectedCrossPairs_kernel <- function(X, Y=NULL, hx, hy, sigma = bw.CvL(X)) {
    autof <- is.null(Y)
    wX <- diggle_weights(X, sigma)
    wY <- if (autof) wX else diggle_weights(Y, sigma)
    if (autof) Y <- X

    res <- numeric(length(hx))

    cfac <- 1/(4*pi*sigma^2)

    # determine integration bounds, from W \cap \W_{-h}
    xhi <- rep(1, length(hx))
    xhi[hx > 0] <- 1 - hx[hx > 0]
    xlo <- numeric(length(hx))
    xlo[hx < 0] <- -hx[hx < 0]

    yhi <- rep(1, length(hy))
    yhi[hy > 0] <- 1 - hy[hy > 0]
    ylo <- numeric(length(hy))
    ylo[hy < 0] <- -hy[hy < 0]

    # Double sum over points of X, points of Y
    # If X == Y (i.e. if (autof)), skip the i == j term
    for (i in 1:length(wX)) {
        xX <- X$x[i]
        yX <- X$y[i]

        js <- if (!autof) 1:length(wY) else which((1:length(wY)) != i)
        for (j in js) {
            xY <- Y$x[j]
            yY <- Y$y[j]

            res <- res + (wX[i] * wY[j]) *
                exp(- (xX - xY + hx)^2/4/sigma^2) *
                exp(- (yX - yY + hy)^2/4/sigma^2) *
                (pnorm( sqrt(2)/sigma*(xhi - (xX + xY - hx)/2)) - pnorm(sqrt(2)/sigma*(xlo - (xX + xY - hx)/2)) )*
                (pnorm( sqrt(2)/sigma*(yhi - (yX + yY - hy)/2)) - pnorm(sqrt(2)/sigma*(ylo - (yX + yY - hy)/2)) )
        }
    }

    cfac * res
}

