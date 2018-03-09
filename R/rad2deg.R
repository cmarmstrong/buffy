#' Radians to degrees
#'
#' Convert radians to degrees or vice-versa
#'
#' @param x numeric radians or degrees
#' @return numeric degrees or radians
rad2deg <- function(x) (x * 180) / (pi)

#' @describeIn rad2deg Convert degrees to radians
deg2rad <- function(x) (x * pi) / (180)
