#' Radians to degrees
#'
#' Convert radians to degrees
#'
#' @param x numeric radians or degrees
#' @return numeric degrees or radians
rad2deg <- function(rad) (rad * 180) / (pi)

#' @describeIn rad2deg Convert degrees to radians
deg2rad <- function(deg) {(deg * pi) / (180)}
