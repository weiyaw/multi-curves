rm(list = ls())
setwd("~/Dropbox/master/algo/")
growth <- reshape2::melt(fda::growth[-3])
growth <- with(growth, data.frame(x = Var1,
                                  y = value,
                                  grp.sub = Var2,
                                  grp.pop = L1))
source("main-tpf.R")
set.seed(100)
system.time(fm2 <- SubjectsTpfMul(growth, 5, deg = 2, shape = "increasing", size = 10000, burn = 0, verbose = TRUE))
saveRDS(fm2, "simulations/multi/multi-growth-uncon.rds")

