\name{Kglobal}
\alias{Kglobal}
\alias{Kcross.global}
\title{ (cross) K functions with a global intensity reweighting }
\description{
Compute \eqn{K_\textrm{global}}{K[global]}
}
\usage{
Kglobal(X, lambdaX=NULL, ... r=NULL, breaks=NULL, normtol=.001,
    discrete.lambda=FALSE, discrete.h=FALSE)

Kcross.global(X, Y, lambdaX=NULL, lambdaY=NULL, ..., r=NULL, breaks=NULL,
    normtol=.001, discrete.lambda=FALSE, discrete.h=FALSE)

}
\arguments{
  \item{X, Y}{
    point process of type ppp, on which to evaluate the (cross) K function
}
  \item{lambdaX, lambdaY}{
    intensity function estimates corresponding to X and Y. If omitted, intensity
    functions will be computed using density.ppp or densityfun.ppp (see discrete.lambda
    below)
}
  \item{...}{
    extra args passed to density.ppp or densityfun.ppp, if applicable.
}
  \item{r}{
    values of r to evaluate K(r) at. If omitted, a sensible default is chosen.
}
  \item{breaks}{
    For internal use only
}
  \item{normtol}{
    A tolerance to use for expectedPairs or expectedCrossPairs when computing
    the normalizing factor.
}
  \item{discrete.lambda}{
    if TRUE, and lambdaX or lambdaY are not supplied, evaluate lambda by
    interpolating the values on a discrete lattice (using interp.im and
    density.ppp), instead of exactly (using densityfun.ppp)
}
  \item{discrete.h}{
    if TRUE, evaluate the expectedCrossPairs on a lattice and interpolate,
    rather than at the exact displacements observed in the pattern
}
}
\value{
    The return value is a numeric vector corresponding to K(r)
}
\references{
}
\author{
Thomas Shaw <shawtr@umich.edu>
}
\seealso{
\code{\link{expectedPairs}}
}
\examples{
rho <- funxy(function(x,y) 100*(1+x), owin())
X <- rpoispp(rho)
K <- Kglobal(X, normtol=.006)
}
