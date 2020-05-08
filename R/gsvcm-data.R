#' Crash data in Florida
#'
#' Data includes census block groups in Florida with variables as follows: \cr
#' \itemize{
#' \item Offcrsh (off-roadway crash frequencies) \cr
#' \item log.VMT (log of vehicle miles traveled) \cr
#' \item log.Pop (log of total population) \cr
#' \item Rmale (proportion of males) \cr
#' \item Rold (proportion of people age 65 and older) \cr
#' \item Rhisp (proportion of Hispanics) \cr
#' \item Runemploy (unemployed rate) \cr
#' \item Lon (longitude of a location) \cr
#' \item Lat (latitude of a location)
#' }
#'
#' @docType data
#'
#' @usage data(Crash_Florida)
#'
#' @format A data frame with 11,249 rows and 9 variables.
#'
#' @keywords datasets
#'
#' @references Myungjin Kim and Lily Wang (2020+), Generalized Spatially Varying Coefficient Models. Journal of Computational and Graphical Statistics.
#'
#' @source
#' \url{https://www.census.gov/geo/maps-data/data/tiger-data.html} \cr
#' \url{https://www.arcgis.com/home/}\cr
#' \url{http://www.fdot.gov/statistics/gis/}\cr
#'
#' @examples
#' data(Crash_Florida)
#' crash <- Crash_Florida$Offcrsh
#' hist(crash)
#' summary(crash)

"Crash_Florida"
