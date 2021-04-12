#' Spartaco co-clustering
#'
#' This function will return the estimated parameters and the co-cluster labels.
#'
#' @import SpatialExperiment
#' @export
#'
#' @examples
#' library(SpatialExperiment)
#' example(SpatialExperiment)
#' spartaco(se)
spartaco <- function(spe,
                     K = 2,
                     R = 4,
                     method = c("default", "marginal")
                     ) {
    x <- assay(spe)
    coordinates <- spatialCoords(spe)

    Dist <- as.matrix(stats::dist(coordinates))

    if(method == "default") {
        main(x, K, R, coordinates, Dist)
    } else if(method == "marginal") {
        stop("Not implemented yet.")
    }

}
