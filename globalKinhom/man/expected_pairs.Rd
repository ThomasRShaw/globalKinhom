\name{expectedPairs}
\alias{expectedPairs}
\alias{expectedCrossPairs}
\alias{cross_f_lattice}
\title{
Expected pairs in an inhomogeneous poisson process
}
\description{
Compute the expected number of pairs at a given displacement h in a poisson
process with a given intensity function. This corresponds to the integrals f and
f_{12} of our paper.
}
\usage{
expectedPairs(rho, hx, hy=NULL, method=c("mc", "lattice"),
                tol=.005, dx=diff(as.owin(rho1)$xrange)/200, maxeval=1e6,
                maxsamp=5e3)

expectedCrossPairs(rho1, rho2=NULL, hx, hy=NULL, method=c("mc", "lattice"),
                tol=.005, dx=diff(as.owin(rho1)$xrange)/200, maxeval=1e6,
                maxsamp=5e3)

cross_f_lattice(rho1, rho2, hx, hy=NULL, dx=.01, ...)
}
\arguments{
  \item{rho1, rho2, rho}{
    Intensity functions, either of class im or funxy.
    This may be produced by density.ppp or densityfun.ppp.
  }
  \item{hx, hy}{
    The displacements h to evaluate the integral at. These can be in any
    format supported by xy.coords.
  }
  \item{method}{
    Either mc (the default) or lattice. Compute integral using monte-carlo or
    on a lattice.
  }
  \item{tol}{
    A tolerance for how precise the integral should be. This is compared to a
    standard error for the mc estimate.
  }
  \item{maxeval}{
    Maximum number of evaluations of rho per iteration. Prevents memory-related
    crashes that can occur.
  }
  \item{maxsamp}{
    Maximum number of monte carlo samples per iteration. If this is too large,
    you may do more work than required to achieve tol.
  }
  \item{dx}{
    if method=="lattice", a lattice spacing for the computation. defaults to .01.
  }
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The return value is a numeric vector with length equal to the number of
  displacements h passed 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Thomas Shaw <shawtr@umich.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
hx <- runif(100, -.2, .2)
hy <- runif(100, -.2, .2)
rho <- funxy(function(x,y) x + 1, owin())
f <- expectedPairs(rho, hx, hy, tol=.004)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
