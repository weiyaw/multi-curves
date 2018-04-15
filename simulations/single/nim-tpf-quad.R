source("main-tpf.R")
sitka <- read.table("data/sitka5.txt", header = T)
sitka <- data.frame(x = sitka$days / 674,
                    y = sitka$log.size,
                    grps = factor(sitka$id.num))

system.time(fm1 <- SubjectsTpf(sitka, 5, deg = 2, size = 10000, burn = 1))
saveRDS(fm1, "tpf-quad.rds")
