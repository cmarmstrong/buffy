#' Construct a rectangular buffer
#'
#' Construct a buffer that is a rectangle with a major and minor length equal to
#' twice the values of \code{d1} and \code{d2} respectively.
#'
#' If parameter \code{d2} is missing, then constructs a square with length
#' \code{d1}.
#'
#' @param d1 half length of top and
#' @param d2
#' @param b bearing of \code{d1}
rectBuff <- function(pnt, d, b) {
    pntE <- pnt + c(0, 1)
    pntN <- pnt + c(1, 0)
    spPnt <- as(pnt, 'Spatial')
    spPntE <- as(pntE, 'Spatial')
    spPntN <- as(pntN, 'Spatial')
    dE <- dist2Degrees(d, spPnt, spPntE)
    dN <- dist2Degrees(d, spPnt, spPntN)
    pntNE <- pnt + c(dN, dE)
    pntNW <- pnt + c(dN, -dE)
    pntSE <- pnt + c(-dN, dE)
    pntSW <- pnt + c(-dN, -dE)
    st_make_grid(st_sfc(c(pntNW, pntNE, pntSE, pntSW)), n=1)
}
