#' Construct rectangular buffers
#'
#' Construct rectangular buffers around points of \code{x} with side-lengths
#' twice the values of \code{d1} and \code{d2}.  The axis defined by \code{d1}
#' will be aligned parallel with bearing \code{b}.  If \code{b} is missing,
#' then it is assumed to be 0 degrees.  If parameter \code{d2} is missing, then
#' this function constructs a square with length \code{d1}.
#'
#' @param x point(s) of longitude-latitude in a row vector or matrix of two
#'   columns, or an object of class SpatialPoints*
#' @param b bearing(s) of \code{d1} in degrees
#' @param d1 half length of bearing aligned distance(s) of class units
#' @param d2 half length of bearing orthogonal distances(s) of class units
#' @return A SpatialPolygons object
#' @seealso \link{\code{geosphere::bearing}}
#' @examples
#' p11 <- c(-92.44744828628, 34.566107548536)
#' p12 <- p11 + rnorm(2, sd=0.01)
#' p21 <- p12 + rnorm(2, sd=0.01)
#' p22 <- p21 + rnorm(2, sd=0.01)
#' b1 <- geosphere::bearing(p11, p12)
#' b2 <- geosphere::bearing(p21, p22)
#' sp::plot(rectBuff(p11, units::ud_units $mi, 2*units::ud_units $mi, b1))
#' sp::plot(with(units::ud_units, rectBuff(rbind(p11, p21), mi, 2*mi,
#'                                         c(b1, b2))))
#' @export
rectBuff <- function(x, d1, d2, b=0) {
    if(missing(d2)) d2 <- d1
    x1 <- geosphere::destPoint(x, b, d1)
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
}