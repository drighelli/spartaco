library(devtools)
load_all()
load("dev/Simulation_combtau1_combscales1_difftauscales.Rdata")

Dist <- as.matrix(stats::dist(Simulation$coordinates))

# version without save.image: 1891
# save.image version: 1913
set.seed(221)
system.time(res2 <- main(Simulation$x, K=2, R=4, Simulation$coordinates, Dist,
            traceRatio = 10,
            max.iter = 10,
            metropolis.iterations = 15,
            estimate.iterations = 10,
            verbose = TRUE))

# saveRDS(res, file = "dev/res_v1.rds")
