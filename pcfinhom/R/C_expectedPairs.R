expCrossPairs_c <- function(hx, hy, tol) {
    .Call("expectedCrossPairs", hx, hy, tol)
}
