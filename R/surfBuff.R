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
#' outer = matrix(c(0,0,10,0,10,10,0,10,0,8,8,8,8,2,0,2,0,0),ncol=2, byrow=TRUE)
#' hole = matrix(c(4,4,4,6,6,6,6,4,4,4),ncol=2, byrow=TRUE)
#' ## lstr <- matrix(c(5,11,5,-1),ncol=2,byrow=TRUE)
#' p1 <- c(5, 11)
#' p2 <- c(2, 5)
#' p <- st_sf(geom=st_sfc(st_polygon(list(outer)), st_polygon(list(hole))), s=c(2, 3))
#' x <- st_sfc(st_point(p1), st_point(p2))
#' surfBuff(x, p, 4)
#' }
#' @export
surfBuff <- function(x, p, d, nQuadSegs=30, ...) { ## TODO: handles overlapping polygons?
    browser() # debug with unprojected data
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
    ## trace linestring intersections and attenuate linestrings by surface effects
    mNew <- mapply(function(lstring, mls) {
        ## browser() # debug switch
        mls <- sf::st_cast(mls, 'MULTILINESTRING')
        coords4326 <- sf::st_coordinates(sf::st_transform(lstring[1, ], 4326)) # start and end points; 4326
        b <- geosphere::bearing(coords4326) # bearing requires 4326
        quad <- ceil(b * 0.011111111111111) # degrees to quadrant
        coordsLs <- sf::st_coordinates(lstring[1, ]) # start and end points
        coordsLs[, 'L1'] <- c(0, 0)
        coordsMls <- sf::st_coordinates(mls)
        coordsLs <- cbind(coordsLs, L2=c(coordsMls[1, 'L2'], coordsMls[nrow(coordsMls), 'L2']))
        coordsMls <- rbind(coordsLs[1, ], coordsMls, coordsLs[2, ])
        coordsMls <- switch(quad, # sort mls by distance along bearing
               coordsMls[order(coordsMls[, 'X']), ],
               coordsMls[-order(coordsMls[, 'X']), ],
               coordsMls[order(coordsMls[, 'Y']), ],
               coordsMls[-order(coordsMls[, 'Y']), ])
        lLs <- lapply(2:nrow(coordsMls), function(i) { # build linestrings with surface effect
            sfgLs <- sf::st_linestring(rbind(coordsMls[i-1, ], coordsMls[i, ]))
            sf::st_sf(sf::st_sfc(sfgLs), s=mls[coordsMls[, 'L2'], ] $s)
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
        coordsXs <- sf::st_coordinates(sf::st_transform(lsXs, 4326))[, c('X', 'Y')]
        b <- geosphere::bearing(coordsXs)
        newPnt <- geosphere::destPoint(coordsXs, b, lenOk)[1, ]
        sfcPnt <- sf::st_sfc(sf::st_point(newPnt))
        sf::st_crs(sfcPnt) <- 4326
        ## browser() ## TODO: intersections with multiple surface effects
        sf::st_coordinates(sf::st_transform(sfcPnt, crsBuf))
    }, split(sfLsIn, with(sfLsIn, list(idP, L2))), split(sfIn, with(sfIn, list(idP, L2))))
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
