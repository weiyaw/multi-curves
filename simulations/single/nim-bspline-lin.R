source("main-bspline.R")
sitka <- read.table("data/sitka5.txt", header = T)
sitka <- data.frame(x = sitka$days / 674,
                    y = sitka$log.size,
                    grps = factor(sitka$id.num))

system.time(fm1 <- SubjectsBs(sitka, 5, deg = 1, size = 10000, burn = 1))
saveRDS(fm1, "bspline-lin.rds")
