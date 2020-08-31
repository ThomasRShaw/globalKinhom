ep_excess <- function(X,hx,hy,sigma,tol=.005, overcount=NULL) {

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

    newterms <- colSums(iterm*jterms)

    result[theseh] <- result[theseh] + newterms * w
    result2[theseh] <- result2[theseh] + (newterms * w)^2

    num_samples[theseh] <- num_samples[theseh] + 1

    if (is.null(overcount)) {
        sd_est <- ( sqrt( (num_samples[theseh]*result2[theseh] - result[theseh]^2) / (num_samples[theseh] - 1))
                        / result[theseh])
    } else {
        sd_est <- ( sqrt( (num_samples[theseh]*result2[theseh] - result[theseh]^2) / (num_samples[theseh] - 1))
                        / (num_samples[theseh]*overcount[theseh]/fac/winwts[theseh] - result[theseh]))
    }

    sd_est[num_samples[theseh] < 10] <- Inf


    #if (min(num_samples) > 1000) break
    #CoV[theseh] <- 1000*tol/num_samples[theseh]
    CoV[theseh] <- sd_est
}

print(c(min(num_samples), mean(num_samples), max(num_samples), sd(num_samples)))

result/num_samples*winwts*fac

 
}

ep_excess_iso <- function(X, r, sigma, tol=.005, overcount=NULL) {

    if (is.function(sigma)) sigma <- sigma(X)
    stopifnot(is.numeric(sigma))

    res <- numeric(length(r))
    res2 <- numeric(length(r))
    resn <- numeric(length(r))

    nr <- 2*ceiling(2*pi*r/sigma) + 1 # sample, nr equally spaced angles with spacing < sigma

    nr[is.infinite(nr)] <- 1

    if (!is.null(overcount)) {
        oc <- overcount
        overcount <- do.call(c, Map(function(n,oc) rep(oc,n), as.list(nr), as.list(overcount)))
    }


    dthetas <- vector('list', length(r))
    for (i in 1:length(r)) dthetas[[i]] <- 2*pi*(0:(nr[i]-1))/nr[i]

    last <- cumsum(nr)
    first <- last - nr + 1

    th0 <- runif(1,0,2*pi)
    th1 <- th0 + pi

    hx0 <- do.call(c, Map(function(dth, ri) ri*cos(th0 + dth), dthetas, r))
    hy0 <- do.call(c, Map(function(dth, ri) ri*sin(th0 + dth), dthetas, r))

    hx1 <- do.call(c, Map(function(dth, ri) ri*cos(th1 + dth), dthetas, r))
    hy1 <- do.call(c, Map(function(dth, ri) ri*sin(th1 + dth), dthetas, r))

    tmp0 <- ep_excess(X,hx0,hy0,sigma,tol, overcount)
    tmp1 <- ep_excess(X,hx1,hy1,sigma,tol, overcount)

    res0 <- numeric(length(r))
    res1 <- numeric(length(r))
    for (i in 1:length(r)) {
        res0[i] <- sum(tmp0[first[i]:last[i]])/nr[i]
        res1[i] <- sum(tmp1[first[i]:last[i]])/nr[i]
    }

    res <- res + (res0 + res1)/2

    if (is.null(overcount)) 
        nbad <- sum(abs(res0 - res1)/res > tol*sqrt(2))
    else
        nbad <- sum(abs(res0 - res1)/(oc - res) > tol*sqrt(2))

    if (nbad > 0) {
        warning("nbad ", nbad, ".")
    }

    resn <- resn + 1

    res / resn # note, no 2pi r here, to avoid problems at 0


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
