#' Snap polygons
#'
#' Group polygons together that are within some distance of each other
#'
#' @param x an sf, sfc, or sfg
#' @param d an object of class units
#' @return a vector of group identifiers
snap <- function(x, d) {
    dist <- sf::st_distance(x)
    hc <- hclust(as.dist(dist>d), method='single')
    cutree(hc, h=0.5)
}
