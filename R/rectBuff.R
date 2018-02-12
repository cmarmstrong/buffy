#' Construct a rectangular buffer
#'
#' Construct a buffer that is a rectangle with a major and minor length equal to
#' twice the values of \code{d1} and \code{d2} respectively.  The rectangle's
#' major axis will be aligned parallel with bearing \code{b}.
#'
#' If parameter \code{d2} is missing, then this function constructs a square
#' with length \code{d1}.
#'
#' @param x an object containing center points of buffers of class sfg, sfc, or
#'   sf
#' @param dMajor half length of top and
#' @param dMinor
#' @param b bearing of \code{dMajor}, see \link{\code{geosphere::bearing}}
#' @param d distance
#' @examples
#' p1 <- sf::st_sfc(sf::st_point(as.numeric(c(-92.44744828628, 34.566107548536))))
#' p2 <- p1 + rnorm(2, sd=0.1)
#' b <- geosphere::bearing(as(p1, 'Spatial'), as(p2, 'Spatial'))
#' rectBuff(p1, units::ud_units $mi, 2*units::ud_units $mi, b)
#' @export
rectBuff <- function(x, dMajor, dMinor, b) {
    browser()
    xMajor <- geosphere::destPoint(as(x, 'Spatial'), b, dMajor)
    xMinor <- geosphere::destPoint(as(x, 'Spatial'), b+90, dMinor)
    
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
