#' Construct rectangular buffers
#'
#' Construct rectangular buffers around points of \code{x} with side-lengths
#' twice the values of \code{d1} and \code{d2}.  The axis defined by \code{d1}
#' will be aligned parallel with bearing \code{b}.
#'
#' If parameter \code{d2} is missing, then this function constructs a square
#' with length \code{d1}.
#'
#' @param x an object containing longitude latitude points of class sfg, sfc,
#'   sf, or SpatialPoints
#' @param b bearing of \code{d1} in degrees
#' @param d1 half length of axis aligned direction
#' @param d2 half length of axis orthogonal direction
#' @return A SpatialPolygons object
#' @seealso \link{\code{geosphere::bearing}}
#' @examples
#' p11 <- sf::st_sfc(sf::st_point(as.numeric(c(-92.44744828628, 34.566107548536))))
#' p12 <- p11 + rnorm(2, sd=0.1)
#' p21 <- p12 + rnorm(2, sd=0.5)
#' p22 <- p21 + rnorm(2, sd=0.1)
#' b1 <- geosphere::bearing(as(p11, 'Spatial'), as(p12, 'Spatial'))
#' b2 <- geosphere::bearing(as(p21, 'Spatial'), as(p22, 'Spatial'))
#' sp::plot(rectBuff(p11, b1, units::ud_units $mi, 2*units::ud_units $mi))
#' rectBuff(c(p11, p21), c(b1, b2), units::ud_units $mi, 2*units::ud_units $mi)
#' @export
rectBuff <- function(x, b, d1, d2) {
    if(is.null(d2)) d2 <- d1
    x1 <- geosphere::destPoint(x <- as(x, 'Spatial'), b, d1)
    x2 <- geosphere::destPoint(x, b+90, d2)
    p1 <- geosphere::destPoint(x1, b+90, d2)
    p2 <- geosphere::destPoint(p1, b+180, 2*d1)
    p3 <- geosphere::destPoint(p2, b+270, 2*d2)
    p4 <- geosphere::destPoint(p3, b, 2*d1)
    buffers <- mapply(function(p1, p2, p3, p4) {
        sp::Polygon(rbind(p1, p2, p3, p4))
    }, split(p1, row(p1)), split(p2, row(p2)), split(p3, row(p3)), split(p4, row(p4)),
    SIMPLIFY=FALSE)
    sp::SpatialPolygons(list(sp::Polygons(buffers, 'rectangular buffers')))
    ## sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(rbind(p1, p2, p3, p4))), 'rectBuff')))
}
