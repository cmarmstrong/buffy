#' Construct a surface buffer
#'
#' Construct a buffer that is attenuated by effects of the surface.
#'
#' @param x points of class sf
#' @param p polygons of class sf with a factor indicating surface type
#' @param s a named list of named vectors; list names indicate factor names in
#'   \code{p}, vector names indicate factor values, and vector values indicate
#'   surface effects
#' @param d radius distance of class units
#' @examples
#' data(urbanUS)
#' x <- c(-92.44744828628, 34.566107548536)
#' mi <- units::ud_units $mi
#' m <- units::ud_units $m
#' aabb <- rectBuff(x, 7e3*m, 7e3*m)
#' sp::proj4string(aabb) <- sp::CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
#' query <- osmdata::opq(sp::bbox(aabb))
#' query <- osmdata::add_osm_feature(query, 'shop', 'supermarket', value_exact=FALSE)
#' osmSupermarkets <- osmdata::osmdata_sf(query) # NOTICE switch to sf
#' sfSupermarkets <- osmSupermarkets $osm_points
#' sfSupermarkets <- sf::st_transform(sfSupermarkets, 3083) # Albers EAC for buffering
#' forCrop <- sf::st_buffer(sfSupermarkets, 11*mi)
#' rCropped <- raster::crop(urbanUS, as(forCrop, 'Spatial'))
#' p <- raster::rasterToPolygons(rCropped, digits=7, dissolve=TRUE)
#' p <- sf::st_as_sf(p)
#' x <- sf::st_sf(geometry=sf::st_sfc(sf::st_point(x)))
#' sf::st_crs(x) <- 4326
#' x <- sf::st_transform(x, 3083)
#' s <- list(urbanUS=c('1'=0.1))
#' surfBuff(x, p, s, 10*mi)
#' @export
surfBuff <- function(x, p, s, d) {
    if(sf::st_crs(x)!=sf::st_crs(p)) stop('CRS of x and p do not match')
    else crsBuf <- sf::st_crs(x)
    buf <- sf::st_buffer(x, d)
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
        sfcLs <- do.call(sf::st_sfc, lapply(m, sf::st_linestring))
        sf::st_sf(geometry=sfcLs, idLs=1:length(sfcLs), idFeature=i)
    })
    sfLs <- do.call(rbind, lSf)
    sf::st_crs(sfLs) <- crsBuf
    ## intersect linestrings with "terrain" polygons
    sfIn <- sf::st_intersection(sfLs, p)
    ## attenuate linestring lengths by "terrain"
    newBuffer <- mapply(function(lstring, mls) {
        browser()
        coordsLs <- sf::st_coordinates(lstring)
        coordsLs[, 'L1'] <- c(0, 0)
        coordsLs <- cbind(coordsLs, L2=c(1, 1))
        coordsMls <- rbind(coordsLs[1, ], sf::st_coordinates(mls), coordsLs[2, ])
        coordsMls <- coordsMls[, c('X', 'Y')]
        ## 
        sfLs <- sf::st_sf(urban=c(rep(c(0, 1), (nrow(coordsMls)-1)%/%2), 0),
                          geometry=sf::st_sfc(lapply(nrow(coordsMls):2, function(i) {
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