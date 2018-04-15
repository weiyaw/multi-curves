rm(list = ls())
setwd("~/Dropbox/master/algo/")
sitka <- read.table("data/sitka.txt", header = T)
sitka <- with(sitka, data.frame(x = days / 674,
                                  y = log.size,
                                  grp.sub = id.num,
                                  grp.pop = ozone))
source("main-tpf.R")
set.seed(100)
system.time(fm1 <- SubjectsTpfMul(sitka, 5, deg = 2, shape = "increasing", size = 10000, burn = 0, verbose = TRUE))
saveRDS(fm1, "simulations/multi/multi-sitka-uncon.rds")








