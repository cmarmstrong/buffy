#' Construct a surface buffer
#'
#' Construct a buffer that is attenuated by effects of the surface.
#'
#' @param x points of class sf
#' @param p polygons of class sf with a factor indicating surface type
#' @param s a named list of named vectors; list names indicate factor names in
#'   \code{p}, vector names indicate factor values, and vector values indicate
#'   surface effects as distance multipliers
#' @param d radius distance of class units; the units must match the units of p
#' @return An sfc object holding surface buffers for each of the points in x
#' @examples
#' x <- c(-92.44744828628, 34.566107548536)
#' mi <- units::ud_units $mi
#' m <- units::ud_units $m
#' aabb <- sp::bbox(rectBuff(x, 1e4*m))
#' rx <- runif(20, aabb[1, 1], aabb[1, 2])
#' ry <- runif(20, aabb[2, 1], aabb[2, 2])
#' p <- rectBuff(cbind(rx, ry), 1e2*m)
#' ## TODO generate two random 's' categories for the polygons in p
#' rx <- runif(5, aabb[1, 1], aabb[1, 2])
#' ry <- runif(5, aabb[2, 1], aabb[2, 2])
#' x <- cbind(rx, ry)
#' surfBuff(x, p, s, 1e3*m)
#' 
#' data(urbanUS)
#' x <- c(-92.44744828628, 34.566107548536)
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
#' s <- list(urbanUS=c('1'=10))
#' d <- 10*mi
#' units(d) <- m
#' foodDeserts <- surfBuff(sfSupermarkets, p, s, d)
#' plot(rCropped, col='orange', legend=FALSE)
#' plot(sf::st_geometry(p), add=TRUE)
#' plot(sf::st_geometry(foodDeserts), add=TRUE)
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
    lM <- lapply(lCoords, function(coords) { # convert coords to 2x2 matrices
        lapply(split(t(coords), seq(NCOL(coords))), matrix, nrow=2)
    })
    lSf <- lapply(1:length(lM), function(i) { # make linestrings from 2x2 matrices
        m <- lM[[i]]
        sfcLs <- do.call(sf::st_sfc, lapply(m, sf::st_linestring))
        sf::st_sf(geometry=sfcLs, idLs=1:length(sfcLs), idFeature=i)
    })
    sfLs <- do.call(rbind, lSf)
    sf::st_crs(sfLs) <- crsBuf
    sfIn <- sf::st_intersection(sfLs, p)
    newBuffer <- mapply(function(lstring, mls) { # attenuate linestrings by surface effects
        browser()
        mls <- sf::st_cast(mls, 'MULTILINESTRING')
        coordsLs <- sf::st_coordinates(lstring) # start and end point
        coordsLs[, 'L1'] <- c(0, 0)
        coordsLs <- cbind(coordsLs, L2=c(1, 1)) # NOTE mls with > 1 L2 ?
        coordsMls <- rbind(coordsLs[1, ], sf::st_coordinates(mls), coordsLs[2, ])
        sfLs <- lapply(unique(coordsMls[, 'L2']), function(L2) {
            coordsSub <- coordsMls[coordsMls[, 'L2']==L2, c('X', 'Y')]
            s <- sapply(names(s), function(name) { # surface types
                ## NOTE insert distance factors, instead of identifiers?  see not in 'scaling'
                c(rep(c(0, mls[L2, name, drop=TRUE]), (nrow(coordsSub)-1)%/%2), 0)
            })
            geom <- lapply(nrow(coordsSub):2, function(i) {
                sf::st_linestring(rbind(coordsSub[i, ], coordsSub[i-1, ]))
            })
            sf::st_sf(geometry=sf::st_sfc(geom), s)
        })
        sfLs <- do.call(rbind, sfLs)
        sf::st_crs(sfLs) <- crsBuf
        lenLs <- sf::st_length(sfLs)
        lenLs <- sapply(names(s), function(name) { # scaling by surface types
            ## NOTE see note in 'surface types'
            isSurf <- sfLs[, name, drop=TRUE]==1
            lenLs[isSurf] <- lenLs[isSurf] * s[[name]]
            lenLs
        })
        csumLs <- cumsum(lenLs)
        units(csumLs) <- units(d)
        xsLs <- csumLs > d
        lsXs <- sfLs[xsLs, ][1, ]
        lenXs <- csumLs[xsLs][1] - d
        ## here
        lenStart <- rev(csumLs[!xsLs])[1]
        lenOk <- d - lenStart
        xsCoords <- sf::st_coordinates(sf::st_transform(lsXs, 4326))[, c('X', 'Y')]
        b <- geosphere::bearing(xsCoords)
        newPnt <- geosphere::destPoint(xsCoords, b, lenOk)[1, ]
        sfcPnt <- sf::st_sfc(sf::st_point(newPnt))
        sf::st_crs(sfcPnt) <- 4326
        sf::st_coordinates(sf::st_transform(sfcPnt, crsBuf))
    }, split(sfLs, 1:nrow(sfLs)), split(sfIn, 1:nrow(sfIn)))
    newBuffer <- t(newBuffer)
    lNewBuffer <- lapply(split(newBuffer, rep(1:29, each=121)), matrix, ncol=2)
    sfcNewBuffers <- sf::st_sfc(lapply(lNewBuffer, function(newBuffer) {
        sf::st_polygon(list(rbind(newBuffer, newBuffer[1, ])))}))
    sf::st_crs(sfcNewBuffers) <- crsBuf
    sfcNewBuffers
}
