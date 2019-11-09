require(spatstat)
require(pcfinhom)

rho.unif <- funxy(function(x,y) numeric(length(x)) + 1, owin());

f.true <- function(hx,hy) (1 - abs(hx))*(1 - abs(hy));

set.seed(41)
hx <- runif(100, -.2, .2);
hy <- runif(100, -.2, .2);
# f.comp.funxy <- f_est(rho.unif, hx, hy, dx=.005)
# 
# f.comp.funxy - f.true(hx,hy)
# 
# rho.unif.im <- as.im(matrix(1, nrow=1, ncol=1), owin())
# f.comp.im <- f_est(rho.unif.im, hx, hy, dx=.005)
# 
# f.comp.im - f.true(hx, hy)

f.comp.mc <- expectedPairs(rho.unif, hx, hy, tol=.004)

f.comp.mc - f.true(hx, hy)

sum((abs(f.comp.mc - f.true(hx, hy))) > .006) / length(hx)
