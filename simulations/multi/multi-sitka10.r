rm(list = ls())
setwd("~/Dropbox/master/algo/")
sitka10 <- read.table("data/sitka10.txt", header = T)
sitka10 <- with(sitka10, data.frame(x = days / 674,
                                    y = log.size,
                                    grp.sub = id.num,
                                    grp.pop = ozone))
source("main-tpf.R")
set.seed(100)
system.time(fm1 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 10000, burn = 0, verbose = TRUE))
saveRDS(fm1, "simulations/multi/multi-sitka10-uncon.rds")








