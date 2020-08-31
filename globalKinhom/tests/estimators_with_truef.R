require(globalKinhom)
require(parallel)

options(mc.cores=4)

set.seed(49)

type="poisson"
prof="gausspeak"
npoint=300

lambda <- funxy(function(x,y) 300*(.5*exp(-((x-.5)^2 + (y-.5)^2)*2) + .5), owin())

hxs <- seq(0, .3, length.out=151)
hx <- outer(hxs, hxs, function(x,y) x)
hy <- outer(hxs, hxs, function(x,y) y)
expected_pairs <- expectedPairs(lambda, hx, hy, tol=.001)
dim(expected_pairs) <- rep(length(hxs),2)
ep.im <- as.im(t(expected_pairs), xcol=hxs, yrow=hxs)
ep <- function(dx,dy) interp.im(ep.im, abs(dx),abs(dy))

r_check <- hxs
expected_pairs_iso <- expectedPairs_iso(lambda, r_check, tol=.001)
ep_iso.spl <- smooth.spline(r_check, expected_pairs_iso, df=length(r_check))
ep_iso <- function(r) predict(ep_iso.spl, r)$y
print("done with expected pairs")

nsim <- 1000
Xs <- rpoispp(lambda, nsim=nsim)
Ys <- rpoispp(lambda, nsim=nsim)
 
out <- data.frame(type=rep("poisson", nsim))
out$X <- Xs
out$Y <- Ys

out$g.truef <- Map(pcfinhom, Xs, MoreArgs=list(lambda=lambda, exp_prs=ep_iso))
print("done with g.truef")
out$c.truef <- Map(pcfcrossinhom, Xs, Ys,
        MoreArgs=list(lambdaX=lambda, lambdaY=lambda, exp_prs=ep_iso))
print("done with c.truef")
out$K.truef <- Map(Kglobal, Xs, MoreArgs=list(lambda=lambda, exp_prs=ep))
print("done with K.truef")
out$Kcross.truef <- Map(Kinhomcross, Xs, Ys,
        MoreArgs=list(lambdaX=lambda, lambdaY=lambda, exp_prs=ep))
print("done with Kcross.truef")

get.fv <- function(out, fn, nm)
    do.call(rbind, lapply(out[[fn]], function(s) s[[nm]]))

r <- out$g.truef[[1]]$r

# These are ~1se above 1
gg <- get.fv(out, "g.truef", "global")
gl <- get.fv(out, "g.truef", "local")
plot(r, colMeans(gg), type='l', ylim=c(.99, 1.01))
lines(r, colMeans(gl), lty=2)
lines(r, colMeans(gg) - apply(gg, 2, sd)/sqrt(nsim))
lines(r, colMeans(gg) + apply(gg, 2, sd)/sqrt(nsim))
lines(r, colMeans(gl) - apply(gl, 2, sd)/sqrt(nsim))
lines(r, colMeans(gl) + apply(gl, 2, sd)/sqrt(nsim))
lines(r, r/r, lty=3)

cg <- get.fv(out, "c.truef", "global")
cl <- get.fv(out, "c.truef", "local") 
plot(r, colMeans(cg), type='l', ylim=c(.9, 1.1))
lines(r, colMeans(cl), lty=2)
lines(r, colMeans(cg) - apply(cg, 2, sd)/sqrt(nsim))
lines(r, colMeans(cg) + apply(cg, 2, sd)/sqrt(nsim))
lines(r, colMeans(cl) - apply(cl, 2, sd)/sqrt(nsim))
lines(r, colMeans(cl) + apply(cl, 2, sd)/sqrt(nsim))
lines(r, r/r, lty=3)

#pi r^2 is well within the +/- sd of these
Kg <- get.fv(out, "K.truef", "global")
Kcg <- get.fv(out, "Kcross.truef", "global")
Kcl <- get.fv(out, "Kcross.truef", "local")


