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
#' \dontrun{
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
#' surfBuff(x, p, 3e3*m)
#' }
#' @export
surfBuff <- function(x, p, d, nQuadSegs=30, ...) { ## TODO: handles overlapping polygons?
    ## buffer and make linestrings from center to buffer points
    if(sf::st_crs(x)!=sf::st_crs(p)) stop('CRS of x and p do not match')
    else crsBuf <- sf::st_crs(x)
    sfBuf <- sf::st_buffer(x, d, nQuadSegs)
    coordsBuf <- sf::st_coordinates(sfBuf)
    ## make linestrings from x to points on respective buf
    lCoords <- by(coordsBuf, coordsBuf[, 'L2'], apply, 1, function(M) {
        rbind(sf::st_coordinates(x[M[4], ]), M[1:2]) # bind center and buffer point coords
    })
    lM <- lapply(lCoords, function(coords) { # convert coords to 2x2 matrices
        lapply(split(t(coords), seq(NCOL(coords))), matrix, nrow=2)
    })
    lSf <- lapply(1:length(lM), function(i) { # make linestrings from 2x2 matrices
        m <- lM[[i]]
        sfcLs <- do.call(sf::st_sfc, lapply(m, sf::st_linestring))
        sf::st_sf(geometry=sfcLs, idP=1:length(sfcLs), L2=i)
    })
    sfLs <- do.call(rbind, lSf)
    sf::st_crs(sfLs) <- crsBuf
    ## intersect linestrings with p and ensure order matches
    sfIn <- sf::st_intersection(sfLs, p)
    sfLsIn <- merge(sfLs, sfIn[, c('idP', 'L2') ,drop=TRUE], by=c('idP', 'L2'))
    sfLsIn <- sfLsIn[with(sfLsIn, order(L2, idP)), ]
    sfIn <- sfIn[with(sfIn, order(L2, idP)), ]
    ## trace linestring intersections and calculate total length given surface effects
    mNew <- mapply(function(lstring, mls) { # attenuate linestrings by surface effects
        mls <- sf::st_cast(mls, 'MULTILINESTRING')
        coordsLs <- sf::st_coordinates(lstring) # start and end points
        coordsLs[, 'L1'] <- c(0, 0)
        coordsMls <- sf::st_coordinates(mls)
        coordsLs <- cbind(coordsLs, L2=c(coordsMls[1, 'L2'], coordsMls[nrow(coordsMls), 'L2']))
        coordsMls <- rbind(coordsLs[1, ], coordsMls, coordsLs[2, ])
        lLs <- lapply(unique(coordsMls[, 'L2']), function(L2) { # build linestrings and effects
            coordsSub <- coordsMls[coordsMls[, 'L2']==L2, c('X', 'Y')]
            s <- c(rep(c(1, mls[L2, ] $s), (nrow(coordsSub)-1)%/%2), 1)
            geom <- lapply(2:nrow(coordsSub), function(i) { # switch iteration!?!
                sf::st_linestring(rbind(coordsSub[i-1, ], coordsSub[i, ])) # switched i's!?!
            })
            sf::st_sf(geometry=sf::st_sfc(geom), s=s)
        })
        sfLs <- do.call(rbind, lLs)
        sf::st_crs(sfLs) <- crsBuf
        lenLs <- sf::st_length(sfLs)
        lenLs <- lenLs * sfLs $s # scale lengths by surface effects
        csumLs <- cumsum(lenLs)
        units(csumLs) <- units(d)
        xsLs <- csumLs > d # the linestrings exceeding d
        if(!any(xsLs)) { # if linestrings < d: return unattenuated buffer point
            return(sf::st_coordinates(sfLs)[nrow(sfLs), c('X', 'Y')])
        }
        lsXs <- sfLs[xsLs, ][1, ] # the xs linestring
        lenXs <- csumLs[xsLs][1] - d # the xs length
        iStart <- rev(which(!xsLs))[1] # index of last ok linestring
        lenStart <- csumLs[iStart] # total length of ok linestrings
        lenOk <- d - lenStart # ok length of xs linestring
        lenOk <- lenOk / sfLs $s[iStart+1] # descale length by surface effect
        xsCoords <- sf::st_coordinates(sf::st_transform(lsXs, 4326))[, c('X', 'Y')]
        b <- geosphere::bearing(xsCoords)
        newPnt <- geosphere::destPoint(xsCoords, b, lenOk)[1, ]
        sfcPnt <- sf::st_sfc(sf::st_point(newPnt))
        sf::st_crs(sfcPnt) <- 4326
        browser() ## TODO: intersections with multiple surface effects
        sf::st_coordinates(sf::st_transform(sfcPnt, crsBuf))
    }, split(sfLsIn, 1:nrow(sfLsIn)), split(sfIn, 1:nrow(sfIn)))
    ## replace original buffer points with attenuated points and make new buffer geometries
    mNew <- t(mNew)
    colnames(mNew) <- c('X', 'Y')
    dfrNew <- cbind(mNew, L1=1, sfIn[, c('L2', 'idP'), drop=TRUE])
    dfrBuf <- data.frame(coordsBuf, idP=rep(1:(1+4*nQuadSegs), nrow(x)))
    dfrAntiNew <- dplyr::anti_join(dfrBuf, dfrNew, by=c('L2', 'idP'))
    dfrNewBuf <- rbind(dfrAntiNew, dfrNew)
    dfrNewBuf <- with(dfrNewBuf, dfrNewBuf[order(idP, L2), ])

    lNewBuf <- by(dfrNewBuf[, c('X', 'Y')], dfrNewBuf $L2, function(x) {
        sf::st_polygon(list(as.matrix(x)))
    }, simplify=FALSE)
    sfcNewBuf <- sf::st_sfc(lNewBuf)
    sf::st_crs(sfcNewBuf) <- crsBuf
    sfcNewBuf
}
