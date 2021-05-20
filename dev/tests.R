library(devtools)
load_all()
load("data/400x400.Rdata")

Dist <- as.matrix(stats::dist(Simulation$coordinates))

set.seed(221)
system.time(res <- main(Simulation$x, K=2, R=2, Simulation$coordinates, Dist,
            traceRatio = 10,
            max.iter = 10,
            metropolis.iterations = 150,
            estimate.iterations = 10,
            verbose = TRUE))

table(Simulation$original.Cs, res$Cs)
table(Simulation$original.Ds, res$Ds)

plot(Simulation$coordinates, pch=19, col=Simulation$original.Ds)
plot(Simulation$coordinates, pch=19, col=res$Ds)

#saveRDS(res, file = "dev/old_res.rds")

old_res <- readRDS("dev/old_res.rds")
table(res$Cs, old_res$Cs)
table(res$Ds, old_res$Ds)
