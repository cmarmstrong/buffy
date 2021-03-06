#' Surface buffer
#'
#' Construct a buffer that is attenuated by effects of the surface.
#' 
#' If p & x have a coordinate reference system, then surfBuff will execute
#' the buffering geometry on the surface of the WSG84 ellipsoid, else the
#' buffering will use euclidean geometry.
#'
#' @param x points of class sfc
#' @param p polygons of class sf, optionally with a single column indicating
#'   surface effects as distance multipliers; eg an effect of 6 increases
#'   distance across a polygon by 6
#' @param d an integer or object of class units for radius distance; if units,
#'   then the units must be convertible to units of p
#' @param nQuadSegs integer; number of segments per quadrant for buffers
#' @return An sfc object of surface buffers for each of the points in x
#' @examples
#' \dontrun{
#' outer = matrix(c(0,0,10,0,10,10,0,10,0,8,8,8,8,2,0,2,0,0),ncol=2, byrow=TRUE)
#' hole = matrix(c(4,4,4,6,6,6,6,4,4,4),ncol=2, byrow=TRUE)
#' p1 <- c(5, 11)
#' p2 <- c(2, 5)
#' p <- sf::st_sf(geom=sf::st_sfc(sf::st_polygon(list(outer)), sf::st_polygon(list(hole))), s=c(2, 3))
#' x <- sf::st_sfc(sf::st_point(p1), sf::st_point(p2))
#' surfBuff(x, p, 8)
#' }
#' @export
surfBuff <- function(x, p, d, nQuadSegs=30) { ## TODO: handles overlapping polygons?
    ## buffer and make linestrings from center to buffer points
    if(sf::st_crs(x)!=sf::st_crs(p)) stop('CRS of x and p do not match')
    else crsBuf <- sf::st_crs(x)
    sfBuf <- sf::st_buffer(x, d, nQuadSegs)
    coordsBuf <- sf::st_coordinates(sfBuf)
    ## make linestrings from x to points on respective buf
    lCoords <- by(coordsBuf, coordsBuf[, 'L2'], apply, 1, function(M) {
        rbind(sf::st_coordinates(x[M[4]]), M[1:2]) # bind center and buffer point coords
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
    lNew <- mapply(function(lstring, mls) {
        mls <- sf::st_cast(mls, 'MULTILINESTRING') # one mls per s
        if(is.na(crsBuf)) { # use trig to get bearing
            p1 <- sf::st_geometry(sf::st_cast(lstring, 'POINT')[1, ])[[1]] # start
            p2 <- sf::st_geometry(sf::st_cast(lstring, 'POINT')[2, ])[[1]] # end
            lP <- as.list(p2 - p1)
            b <- rad2deg(do.call(atan2, rev(lP))) # arctan(y/x)
            b <- (450 - b) %% 360 # convert unit circle to compass
        } else { # get bearing on ellipsoid
            coords4326 <- sf::st_coordinates(sf::st_transform(lstring, 4326))
            coords4326 <- coords4326[, c('X', 'Y')]
            b <- geosphere::bearing(coords4326)[1]
            b <- b%%360
        }
        quad <- floor(b * 0.011111111111111 + 1) # degrees to quadrant
        coordsLs <- sf::st_coordinates(lstring[1, ]) # start and end points
        coordsLs[, 'L1'] <- c(1, 2) # first and second elements of the 'start and end' feature
        coordsLs <- cbind(coordsLs, L2=c(0, 0)) # 'start and end' is feature 0
        coordsMls <- sf::st_coordinates(mls)
        coordsMls <- rbind(coordsLs[1, ], coordsMls, coordsLs[2, ])
        ## order lines close to 180/360 by Y to avoid errors of quadrant difference b/w projections
        ax <- ifelse(isTRUE(all.equal(0, abs(b-180), tolerance=2)), 'Y',
              ifelse(isTRUE(all.equal(0, (b-360)%%360, tolerance=2)), 'Y', 'X'))
        decreasing <- FALSE # set ascending/descending order by quadrant
        if(ax=='Y' & (quad==2|quad==3)) decreasing <- TRUE
        if(ax=='X' & quad > 2) decreasing <- TRUE
        coordsMls <- switch(quad, # order mls by distance along bearing (quads advance along compass)
                            coordsMls[order(coordsMls[, ax], decreasing=decreasing), ],
                            coordsMls[order(coordsMls[, ax], decreasing=decreasing), ],
                            coordsMls[order(coordsMls[, ax], decreasing=decreasing), ],
                            coordsMls[order(coordsMls[, ax], decreasing=decreasing), ])
        lLs <- lapply(2:nrow(coordsMls), function(i) { # build linestrings with surface effect
            coordsL <- rbind(coordsMls[i-1, ], coordsMls[i, ])
            sfgLs <- sf::st_linestring(coordsL[, c('X', 'Y')])
            L2 <- unique(coordsL[, 'L2']) # length always == 1 ?
            ## if line within intersected feature: set s, else s=1 (no surface effect)
            identical_L2 <- do.call(identical, as.list(coordsL[, 'L2'])) # feature
            identical_L1 <- do.call(identical, as.list(coordsL[, 'L1'])) # intersections of feature
            s <- ifelse(!identical_L2 | !identical_L1, 1, mls[L2, ] $s)
            sf::st_sf(geom=sf::st_sfc(sfgLs), s=s)
        })
        sfLs <- do.call(rbind, lLs)
        sf::st_crs(sfLs) <- crsBuf
        lenLs <- sf::st_length(sfLs) # length returns units object, if appropriate
        lenLs <- lenLs * sfLs $s # scale lengths by surface effects
        csumLs <- cumsum(lenLs)
        xsLs <- csumLs > d # the linestrings exceeding d
        if(!any(xsLs)) { # if linestrings < d: return unattenuated buffer point
            sfcP <- sf::st_geometry(sf::st_cast(lstring, 'POINT')[2, ])
            return(sf::st_sf(geom=sfcP, mls[1, c('idP', 'L2'), drop=TRUE]))
        }
        lsXs <- sfLs[xsLs, ][1, ] # the xs linestring
        lenXs <- csumLs[xsLs][1] - d # the xs length
        iStart <- rev(which(!xsLs))[1] # index of last ok linestring
        lenStart <- csumLs[iStart] # total length of ok linestrings
        lenOk <- d - lenStart # ok length of xs linestring
        lenOk <- lenOk / sfLs $s[iStart+1] # descale length by surface effect
        if(is.na(crsBuf)) { # use trig to get new point
            b <- deg2rad((90 - b) %% 360) # convert compass to unit circle
            sfcOk <- sf::st_geometry(sf::st_cast(lsXs, 'POINT')[1, ])
            sfcNew <- sfcOk + c(lenOk * cos(b), lenOk * sin(b)) # class(sfcNew) == class(sfcOk)
        } else { # get new point on ellipsoid
            coordsXs <- sf::st_coordinates(sf::st_transform(lsXs, 4326))[, c('X', 'Y')]
            b <- geosphere::bearing(coordsXs)
            coordsNew <- geosphere::destPoint(coordsXs, b, units::set_units(lenOk, 'm'))[1, ]
            sfcNew <- sf::st_sfc(sf::st_point(coordsNew))
            sf::st_crs(sfcNew) <- 4326
            sfcNew <- sf::st_transform(sfcNew, crsBuf)
        }
        sf::st_sf(geom=sfcNew, mls[1, c('idP', 'L2'), drop=TRUE])
    }, split(sfLsIn, with(sfLsIn, list(idP, L2)), drop=TRUE), # args
    split(sfIn, with(sfIn, list(idP, L2)), drop=TRUE),
    SIMPLIFY=FALSE)
    ## replace original buffer points with attenuated points and make new buffer geometries
    sfNew <- do.call(rbind, lNew)
    dfrBuf <- data.frame(coordsBuf, idP=rep(1:(1+4*nQuadSegs), length(sf::st_geometry(x))))
    dfrAntiNew <- dplyr::anti_join(dfrBuf, sfNew, by=c('L2', 'idP'))
    lAntiNew <- lapply(split(as.matrix(dfrAntiNew[, c('X', 'Y')]), 1:nrow(dfrAntiNew)), sf::st_point)
    sfcAntiNew <- do.call(sf::st_sfc, lAntiNew)
    sf::st_crs(sfcAntiNew) <- crsBuf
    sfAntiNew <- sf::st_sf(geom=sfcAntiNew, dfrAntiNew[, c('idP', 'L2')])
    sfNewBuf <- rbind(sfAntiNew, sfNew)
    sfNewBuf <- sfNewBuf[with(sfNewBuf, order(L2, idP)), ]
    lNewBuf <- tapply(sf::st_geometry(sfNewBuf), sfNewBuf $L2, function(x) {
        sf::st_polygon(list(sf::st_coordinates(x)))
    }, simplify=FALSE)
    sfcNewBuf <- sf::st_sfc(lNewBuf)
    sf::st_crs(sfcNewBuf) <- crsBuf
    sfcNewBuf
}
