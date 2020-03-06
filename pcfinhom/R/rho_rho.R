rhorho <- function(ux, uy, hx, hy, X, sigma, cutoff) {

    inds <- order(X$x)
    coords(X) <- coords(X)[inds,]
    x <- X$x
    y <- X$y

    indsq <- order(ux)
    xq <- ux[indsq]
    yq <- uy[indsq]

    wq <- diggle_weights2(ux, uy, sigma)
    wp <- diggle_weights2(ux + hx, uy + hy, sigma)

    zz <- .C("rho_rho",
        nquery = as.integer(length(ux)),
        xq = as.double(xq),
        yq = as.double(yq),
        ndata = as.integer(length(x)),
        xd = as.double(x),
        yd = as.double(y),
        xh = as.double(hx),
        yh = as.double(hy),
        rmaxi = as.double(cutoff),
        sig= as.double(sigma),
        result = as.double(double(length(ux))),
        package = "pcfinhom")

    res <- zz$result*wq * wp

    res[indsq] <- res
    res

}
