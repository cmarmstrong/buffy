#' Construct a surface buffer
#'
#' Construct a buffer that is attenuated by effects of the surface.
#'
#' @param x points of class sf
#' @param p polygons of class sf, optionally with a single column indicating
#'   surface effects as distance multipliers; eg an effect of 6 increases
#'   distance across a polygon by 6
#' @param d radius distance of class units; the units must be convertible to
#'   units of p
#' @return An sfc object holding surface buffers for each of the points in x
#' @examples
#' mi <- units::ud_units $mi
#' m <- units::ud_units $m
#' x <- c(runif(1, -114, -82), runif(1, 31, 40))
#' aabb <- sp::bbox(rectBuff(x, 1e4*m))
#' rx <- runif(20, aabb[1, 1], aabb[1, 2])
#' ry <- runif(20, aabb[2, 1], aabb[2, 2])
#' p <- rectBuff(cbind(rx, ry), 2e3*m)
#' p <- sf::st_cast(sf::st_as_sf(p), 'POLYGON')
#' isect <- snap(p, 0)
#' p <- do.call(sf::st_sfc, sapply(1:max(isect), function(j) sf::st_union(p[isect==j, ])))
#' s <- data.frame(s=rpois(length(p), 3)+1)
#' p <- sf::st_sf(p, s)
#' ## will a single factor column handle overlapping polygons?
#' ## p <- sp::SpatialPolygonsDataFrame(sp::disaggregate(p), s)
#' ## p <- sf::st_as_sf(p)
#' rx <- runif(5, aabb[1, 1], aabb[1, 2])
#' ry <- runif(5, aabb[2, 1], aabb[2, 2])
#' x <- sf::st_sfc(sf::st_multipoint(cbind(rx, ry)))
#' sf::st_crs(x) <- 4326
#' sf::st_crs(p) <- 4326
#' x <- sf::st_transform(x, 3083)
#' p <- sf::st_transform(p, 3083)
#' x <- sf::st_sf(sf::st_cast(x, 'POINT'))
#' ## TODO union p and average s values; overlapping polygons may cause a problem
#' surfBuff(x, p, 3e3*m)
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
surfBuff <- function(x, p, d) {
    ## browser()
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
        coordsMls <- sf::st_coordinates(mls)
        coordsLs <- cbind(coordsLs, L2=c(coordsMls[1, 'L2'], coordsMls[nrow(coordsMls), 'L2']))
        coordsMls <- rbind(coordsLs[1, ], coordsMls, coordsLs[2, ])
        sfLs <- lapply(unique(coordsMls[, 'L2']), function(L2) {
            ## browser()
            coordsSub <- coordsMls[coordsMls[, 'L2']==L2, c('X', 'Y')]
            s <- c(rep(c(1, mls[L2, ] $s), (nrow(coordsSub)-1)%/%2), 1)
            ## s <- sapply(names(s), function(name) { # surface types
            ##     ## NOTE insert distance factors, instead of identifiers?  see note in 'scaling'
            ##     c(rep(c(0, mls[L2, name, drop=TRUE]), (nrow(coordsSub)-1)%/%2), 0)
            ## })
            geom <- lapply(nrow(coordsSub):2, function(i) {
                sf::st_linestring(rbind(coordsSub[i, ], coordsSub[i-1, ]))
            })
            sf::st_sf(geometry=sf::st_sfc(geom), s=s)
        })
        sfLs <- do.call(rbind, sfLs)
        sf::st_crs(sfLs) <- crsBuf
        lenLs <- sf::st_length(sfLs)
        lenLs <- lenLs * sfLs $s
        ## lenLs <- sapply(names(s), function(name) { # scaling by surface types
        ##     ## NOTE see note in 'surface types'
        ##     isSurf <- sfLs[, name, drop=TRUE]==1
        ##     lenLs[isSurf] <- lenLs[isSurf] * s[[name]]
        ##     lenLs
        ## })
        csumLs <- cumsum(lenLs)
        units(csumLs) <- units(d)
        xsLs <- csumLs > d
        lsXs <- sfLs[xsLs, ][1, ]
        lenXs <- csumLs[xsLs][1] - d
        ## here
        iStart <- rev(which(!xsLs))[1]
        ## lenStart <- rev(csumLs[!xsLs])[1]
        lenStart <- csumLs[iStart]
        lenOk <- d - lenStart
        xsCoords <- sf::st_coordinates(sf::st_transform(lsXs, 4326))[, c('X', 'Y')]
        b <- geosphere::bearing(xsCoords)
        ## descale lenOk: sfLs $s[iStart] ?
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
