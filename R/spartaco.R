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
#' # spartaco(se)
spartaco <- function(spe,
                     K = 2,
                     R = 4,
                     method = c("default", "marginal"),
                     traceRatio = 10,
                     max.iter = 1000,
                     metropolis.iterations = 150,
                     estimate.iterations = 10,
                     verbose = FALSE
                     ) {
    x <- assay(spe)
    coordinates <- spatialCoords(spe)

    Dist <- as.matrix(stats::dist(coordinates))

    if(method == "default") {
        main(x, K, R, coordinates, Dist,
             traceRatio = traceRatio,
             max.iter = max.iter,
             metropolis.iterations = metropolis.iterations,
             estimate.iterations = estimate.iterations,
             verbose = verbose)
    } else if(method == "marginal") {
        stop("Not implemented yet.")
    }

}
