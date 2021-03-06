% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/surfBuff.R
\name{surfBuff}
\alias{surfBuff}
\title{Surface buffer}
\usage{
surfBuff(x, p, d, nQuadSegs = 30)
}
\arguments{
\item{x}{points of class sfc}

\item{p}{polygons of class sf, optionally with a single column indicating
surface effects as distance multipliers; eg an effect of 6 increases
distance across a polygon by 6}

\item{d}{an integer or object of class units for radius distance; if units,
then the units must be convertible to units of p}

\item{nQuadSegs}{integer; number of segments per quadrant for buffers}
}
\value{
An sfc object of surface buffers for each of the points in x
}
\description{
Construct a buffer that is attenuated by effects of the surface.
}
\details{
If p & x have a coordinate reference system, then surfBuff will execute
the buffering geometry on the surface of the WSG84 ellipsoid, else the
buffering will use euclidean geometry.
}
\examples{
\dontrun{
outer = matrix(c(0,0,10,0,10,10,0,10,0,8,8,8,8,2,0,2,0,0),ncol=2, byrow=TRUE)
hole = matrix(c(4,4,4,6,6,6,6,4,4,4),ncol=2, byrow=TRUE)
p1 <- c(5, 11)
p2 <- c(2, 5)
p <- sf::st_sf(geom=sf::st_sfc(sf::st_polygon(list(outer)), sf::st_polygon(list(hole))), s=c(2, 3))
x <- sf::st_sfc(sf::st_point(p1), sf::st_point(p2))
surfBuff(x, p, 8)
}
}
