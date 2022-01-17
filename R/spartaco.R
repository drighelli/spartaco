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
spartaco <- function(x,
                     coordinates,
                     K = 2,
                     R = 4,
                     traceRatio = 10,
                     max.iter = 1000,
                     metropolis.iterations = 150,
                     estimate.iterations = 10,
                     input.values = NULL,
                     verbose = FALSE
                     ) {

    Dist <- as.matrix(stats::dist(coordinates))

    main(x = x, 
         Dist = Dist,
         K = K, R = K, 
         traceRatio = traceRatio,
         max.iter = max.iter,
         metropolis.iterations = metropolis.iterations,
         estimate.iterations = estimate.iterations,
         input.values = input.values,
         verbose = verbose)
}
