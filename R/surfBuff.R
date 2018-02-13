#' Construct a surface buffer
#'
#' Construct a buffer that is attenuated by effects of the surface.
#'
#' Parameters \code{x} and \code{p} must be \code{\link{sf}} objects.
#' 
#' Names of \code{s} must match factor IDs in \code{p}, values of \code{s}
#' indicate resistance through polygons with respective factor ID.
#'
#' Parameter \code{d} must be a \code{\link{units}} object
#'
#' @param x points of class sf
#' @param p polygons of class sf with factor IDs
#' @param s named vector of surface effects for each factor ID in \code{p}
#' @param d radius distance of class units
#' @examples
#' data(urbanUS)
#' x <- c(-92.44744828628, 34.566107548536)
#' mi <- units::ud_units $mi
#' d <- 10 # miles
#' dDegrees <- 5e3/geosphere::distGeo(x, x+c(0, 1)) # meters
#' aabb <- rectBuff(x, dDegrees)
#' sp::proj4string(aabb) <- sp::CRS('+init=espg:4326')
#' query <- osmdata::opq(sp::bbox(aabb))
#' query <- osmdata::add_osm_feature(query, 'shop', 'supermarket', value_exact=FALSE)
#' osmSupermarkets <- osmdata::osmdata_sf(query) # NOTICE switch to sf
#' forCrop <- sf::st_buffer(osmSupermarkets, (d+1)*mi)
#' rCropped <- crop(urbanUS, as(forCrop, 'Spatial'))
#' p <- rasterToPolygons(rCropped, digits=7, dissolve=TRUE)
#' p <- st_as_sf(p)
#' surfBuff(x, p, 0.1, d)
#' @export
surfBuff <- function(x, p, s, d) {
    browser()
    if(sf::st_crs(x)!=sf::st_crs(p)) stop('CRS of x and p do not match')
    else crsBuf <- sf::st_crs(x)
    buf <- sf::st_buffer(x, dist)
    coordsBuf <- sf::st_coordinates(buf)
    ## make linestrings from x to points on respective buf
    lCoords <- by(coordsBuf, coordsBuf[, 'L2'], apply, 1, function(M) {
        rbind(M[1:2], sf::st_coordinates(x[M[4], ]))
    })
    lM <- lapply(lCoords, function(coords) {
        lapply(split(t(coords), seq(NCOL(coords))), matrix, nrow=2)
    })
    lSf <- lapply(1:length(lM), function(i) {
        m <- lM[[i]]
        sfcLs <- do.call(st_sfc, lapply(m, sf::st_linestring))
        sf::st_sf(geometry=sfcLs, idLs=1:length(sfcLs), idFeature=i)
    })
    sfLs <- do.call(rbind, lSf)
    ## intersect linestrings with "terrain" polygons
    sfIn <- sf::st_intersection(sfLs, p)
    ## attenuate linestring lengths by "terrain"
    newBuffer <- mapply(function(lstring, mls) {
        coordsLs <- sf::st_coordinates(lstring)
        coordsLs[, 'L1'] <- c(0, 0)
        coordsLs <- cbind(coordsLs, L2=c(1, 1))
        coordsMls <- rbind(coordsLs[1, ], sf::st_coordinates(mls), coordsLs[2, ])
        coordsMls <- coordsMls[, c('X', 'Y')]
        sfLs <- st_sf(urban=c(rep(c(0, 1), (nrow(coordsMls)-1)%/%2), 0),
                      geometry=st_sfc(lapply(nrow(coordsMls):2, function(i) {
                          sf::st_linestring(rbind(coordsMls[i, ], coordsMls[i-1, ]))
                      })))
        sf::st_crs(sfLs) <- crsBuf
        lenLs <- sf::st_length(sfLs)
        csumLs <- cumsum(ifelse(sfLs $urban==1, lenLs, lenLs/10))
        xsLs <- with(ud_units, csumLs*m > mi)
        lsXs <- sfLs[xsLs, ][1, ]
        ## isUrban <- sfLs[xsLs, 'urban', drop=TRUE][1]
        ## isUrban <- lsXs $urban
        xsLen <- with(ud_units, csumLs[xsLs][1]*m - mi)
        startLen <- rev(csumLs[!xsLs])[1]
        okLen <- set_units(with(ud_units, mi - startLen*m), 'm')
        okLen <- ifelse(lsXs $urban==1, okLen, okLen*10)
        ## xsCoords <- st_coordinates(st_transform(sfLs[xsLs, ], 4326))[, c('X', 'Y')]
        xsCoords <- sf::st_coordinates(sf::st_transform(lsXs, 4326))[, c('X', 'Y')]
        b <- geosphere::bearing(xsCoords)
        ## NOTE: newPnt may not be correct
        newPnt <- geosphere::destPoint(xsCoords, b, okLen)[1, ]
        sfcPnt <- sf::st_sfc(st_point(newPnt))
        sf::st_crs(sfcPnt) <- 4326
        sf::st_coordinates(sf::st_transform(sfcPnt, crsBuf))
    }, split(sfLs, 1:nrow(sfLs)), split(sfIn, 1:nrow(sfIn)))
    newBuffer <- t(newBuffer)
    lNewBuffer <- lapply(split(newBuffer, rep(1:29, each=121)), matrix, ncol=2)
    sfcNewBuffers <- sf::st_sfc(lapply(lNewBuffer, function(newBuffer) {
        sf::st_polygon(list(rbind(newBuffer, newBuffer[1, ])))}))
    foodDeserts <- sf::st_union(sfcNewBuffers)
    st_crs(foodDeserts) <- crsBuf
}