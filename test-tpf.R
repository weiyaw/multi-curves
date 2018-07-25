library(nlme)
rm(list = ls())
setwd("~/Dropbox/master/algo/")
source("main-tpf.R")
source("subor.R")
sitka10 <- read.table("data/sitka10.txt", header = T)
##sitka <- read.table("data/sitka.txt", header = T)



sitka10 <- read.table("data/sitka10.txt", header = T)
sitka10 <- with(sitka10, data.frame(x = days / 674,
                                    y = log.size / 6,
                                    grp.sub = factor(id.num),
                                    grp.pop = ozone))
levels(sitka10$grp.sub) <- c(levels(sitka10$grp.sub)[3:10], levels(sitka10$grp.sub)[1:2])
## linear spline
lm1 <- lmm_tpf(sitka10, K = 5, deg = 1)

## quadratic spline (algorithm diverges)
## only works if the random effect of the quadratic polynomial term is removed
lmm_tpf(sitka10, K = 5, deg = 2)



## Fit a lme model (quadratic)
## sub.pd.q <- pdBlocked(list(pdSymm(~ time.q + I(time.q^2)), pdIdent(~ Z.q - 1)))
## sub.pd.q <- pdSymm(~ time.q + I(time.q^2))
## sub.pd.q <- pdBlocked(list(pdSymm(~ time.q), pdIdent(~ Z.q - 1)))
## quad.fm <- lme(fixed = y ~ time.q + I(time.q^2),
##                random = list(pop.level = pop.pd.q, sub.level = sub.pd.q))



## TEST MIXED MODEL LINEAR SPLINE (SubjectsTpf)
rm(list = ls())
sitka <- read.table("data/sitka5.txt", header = T)
sitka <- with(sitka, data.frame(x = days / 674,
                                y = log.size,
                                grps = id.num))
source("main-tpf.R")
## system.time(fm2 <- get_tpf_old(sitka, 5, 2, size = 100, burn = 0))
saveRDS(fm2, "tpf-long.rds")
source("graphs.R")
fm1 <- readRDS("simulations/single/tpf-lin.rds")
fm2 <- readRDS("simulations/single/tpf-quad.rds")
PlotSpline(fm1, range(sitka$x), sitka)
PlotSpline(fm2, range(sitka$x), sitka)

fm3 <- readRDS("simulations/bspline-lin.rds")
fm4 <- readRDS("simulations/bspline-quad.rds")
PlotSpline(fm3, range(sitka$x), sitka)
PlotSpline(fm4, range(sitka$x), sitka)

## TEST MIXED MODEL LINEAR SPLINE WITH MULTIPLE POPULATION (SubjectsTpfMul)
## LOAD DATA sitka10 and sitka
rm(list = ls())
sitka10 <- read.table("data/sitka10.txt", header = T)
sitka10 <- with(sitka10, data.frame(x = days / 674,
                                    y = log.size,
                                    grp.sub = id.num,
                                    grp.pop = ozone))
sitka10$grp.sub <- factor(sitka10$grp.sub,
                          levels = c("1", "2", "3", "4", "5",
                                     "60", "59", "56", "57", "58"))

sitka <- read.table("data/sitka.txt", header = T)
sitka <- with(sitka, data.frame(x = days / 674,
                                y = log.size,
                                grp.sub = id.num,
                                grp.pop = ozone))

source("main-tpf.R")
source("graphs.R")
## single population (sitka10)
system.time(fm1 <- SubjectsTpf(sitka10, 5, deg = 2, shape = "increasing", size = 1000))
plot_spline(fm1)

## single population (sitka)
system.time(fm2 <- SubjectsTpf(sitka, 5, deg = 2, shape = "increasing", size = 10000, burn = 0))
plot_spline(fm2)
plot(fm2$samples$population[1, ])
plot(fm2$samples$population[2, ])
plot(fm2$samples$population[3, ])
plot(fm2$samples$population[4, ])
plot(fm2$samples$population[5, ])
plot(fm2$samples$population[6, ])
plot(fm2$samples$population[7, ])
plot(fm2$samples$population[8, ])

plot(fm2$samples$subjects[1, 1, ])
plot(fm2$samples$subjects[2, 1, ])
plot(fm2$samples$subjects[3, 1, ])
plot(fm2$samples$subjects[4, 1, ])
plot(fm2$samples$subjects[5, 1, ])
plot(fm2$samples$subjects[6, 1, ])
plot(fm2$samples$subjects[7, 1, ])
plot(fm2$samples$subjects[8, 1, ])


## multiple population (sitka)
system.time(fm3 <- SubjectsTpfMul(sitka, 5, deg = 2, shape = "increasing", size = 2000, burn = 0))

plot(fm1$samples$population$`1`[1, ])
plot(fm1$samples$population$`1`[2, ])
plot(fm1$samples$population$`1`[3, ])
plot(fm1$samples$population$`1`[4, ])
plot(fm1$samples$population$`1`[5, ])
plot(fm1$samples$population$`1`[6, ])
plot(fm1$samples$population$`1`[7, ])
plot(fm1$samples$population$`1`[8, ])

plot(fm1$samples$subjects$`1`[1, 1, ])
plot(fm1$samples$subjects$`1`[2, 1, ])
plot(fm1$samples$subjects$`1`[3, 1, ])
plot(fm1$samples$subjects$`1`[4, 1, ])
plot(fm1$samples$subjects$`1`[5, 1, ])
plot(fm1$samples$subjects$`1`[6, 1, ])
plot(fm1$samples$subjects$`1`[7, 1, ])
plot(fm1$samples$subjects$`1`[8, 1, ])

fm2 <- readRDS("simulations/multi/multi-1k-bugless.rds")
## multiple population (1000 burn, independent start, garbage)
fm3 <- readRDS("simulations/multi/multi-10k.rds")
source("graphs.R")
plot_spline(fm1)
plot_spline(fm2)
plot_spline(fm3)
plot_spline(fm4)

## multiple population (1000 burn, previous start)
source("main-tpf.R")
system.time(fm4 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 1000, burn = 0))

plot(fm4$samples$population$`1`[1, ])
plot(fm4$samples$population$`1`[2, ])
plot(fm4$samples$population$`1`[3, ])
plot(fm4$samples$population$`1`[4, ])
plot(fm4$samples$population$`1`[5, ])
plot(fm4$samples$population$`1`[6, ])
plot(fm4$samples$population$`1`[7, ])
plot(fm4$samples$population$`1`[8, ])

plot(fm4$samples$subjects$`1`[1, 1, ])
plot(fm4$samples$subjects$`1`[2, 1, ])
plot(fm4$samples$subjects$`1`[3, 1, ])
plot(fm4$samples$subjects$`1`[4, 1, ])
plot(fm4$samples$subjects$`1`[5, 1, ])
plot(fm4$samples$subjects$`1`[6, 1, ])
plot(fm4$samples$subjects$`1`[7, 1, ])
plot(fm4$samples$subjects$`1`[8, 1, ])

source("main-tpf.R")
system.time(fm4 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 1000, burn = 0))

fm5 <- readRDS("simulations/multi/multi-sitka.rds")
system.time(fm5 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 10000, burn = 0))

source("graphs.R")
fm7t <- truncate(fm7, 2000)
plot_spline(fm7t)

fm2t <- truncate(fm2, 300)
plot_spline(fm2t)


## Berkeley growth dataset
fm6 <- readRDS("simulations/multi/multi-growth.rds")
fm7 <- readRDS("simulations/multi/multi-growth-uncon.rds")
growth <- reshape2::melt(fda::growth[-3])
growth <- with(growth, data.frame(x = Var1, y = value, grp.sub = Var2, grp.pop = L1))
growth10 <- subset(growth,
                   grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                  "girl01", "girl02", "girl03", "girl04", "girl05"),
                   drop = TRUE)
growth10$grp.sub <- droplevels(growth10$grp.sub)

growth20 <- subset(growth,
                   grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                  "girl01", "girl02", "girl03", "girl04", "girl05",
                                  "boy06", "boy07", "boy08", "boy09", "boy10",
                                  "girl06", "girl07", "girl08", "girl09", "girl10"),
                   drop = TRUE)
growth20$grp.sub <- droplevels(growth20$grp.sub)

set.seed(1)
source("main-tpf.R")
system.time(fm8 <- SubjectsTpf(growth20, 8, deg = 2, shape = "increasing", size = 15000, burn = 0, verbose = T))

plot(fm8$samples$population[1, ])
plot(fm8$samples$population[2, ])
plot(fm8$samples$population[3, ])
plot(fm8$samples$population[4, ])
plot(fm8$samples$population[5, ])
plot(fm8$samples$population[6, ])
plot(fm8$samples$population[7, ])
plot(fm8$samples$population[8, ])
plot(fm8$samples$population[9, ])
plot(fm8$samples$population[10, ])
plot(fm8$samples$population[11, ])

plot(fm8$samples$subjects[1, 1, ])
plot(fm8$samples$subjects[2, 1, ])
plot(fm8$samples$subjects[3, 1, ])
plot(fm8$samples$subjects[4, 1, ])
plot(fm8$samples$subjects[5, 1, ])
plot(fm8$samples$subjects[6, 1, ])
plot(fm8$samples$subjects[7, 1, ])
plot(fm8$samples$subjects[8, 1, ])
plot(fm8$samples$subjects[9, 1, ])
plot(fm8$samples$subjects[10, 1, ])
plot(fm8$samples$subjects[11, 1, ])

plot(1 / fm8$samples$precision$pop)
plot(1 / fm8$samples$precision$sub)
plot(fm8$samples$precision$eps)

plot_spline(fm8)
plot_spline(truncate_spline(fm8, 10000))


set.seed(1)
source("main-tpf.R")
system.time(fm9 <- SubjectsTpf(growth10, 8, deg = 2, shape = "increasing", size = 40000, burn = 0, verbose = T))

plot(fm9$samples$population[1, ])
plot(fm9$samples$population[2, ])
plot(fm9$samples$population[3, ])
plot(fm9$samples$population[4, ])
plot(fm9$samples$population[5, ])
plot(fm9$samples$population[6, ])
plot(fm9$samples$population[7, ])
plot(fm9$samples$population[8, ])
plot(fm9$samples$population[9, ])
plot(fm9$samples$population[10, ])
plot(fm9$samples$population[11, ])

plot(fm9$samples$subjects[1, 1, ])
plot(fm9$samples$subjects[2, 1, ])
plot(fm9$samples$subjects[3, 1, ])
plot(fm9$samples$subjects[4, 1, ])
plot(fm9$samples$subjects[5, 1, ])
plot(fm9$samples$subjects[6, 1, ])
plot(fm9$samples$subjects[7, 1, ])
plot(fm9$samples$subjects[8, 1, ])
plot(fm9$samples$subjects[9, 1, ])
plot(fm9$samples$subjects[10, 1, ])
plot(fm9$samples$subjects[11, 1, ])

plot(1 / fm9$samples$precision$pop)
plot(1 / fm9$samples$precision$sub)
plot(fm9$samples$precision$eps)

plot_spline(fm9)
plot_spline(truncate_spline(fm9, 10000))

set.seed(1)
source("main-tpf.R")
system.time(fm10 <- SubjectsTpf(growth, 8, deg = 2, shape = "increasing", size = 20000, burn = 0, verbose = T))

plot(fm10$samples$population[1, ])
plot(fm10$samples$population[2, ])
plot(fm10$samples$population[3, ])
plot(fm10$samples$population[4, ])
plot(fm10$samples$population[5, ])
plot(fm10$samples$population[6, ])
plot(fm10$samples$population[7, ])
plot(fm10$samples$population[8, ])


plot(fm10$samples$subjects[1, 1, ])
plot(fm10$samples$subjects[2, 1, ])
plot(fm10$samples$subjects[3, 1, ])
plot(fm10$samples$subjects[4, 1, ])
plot(fm10$samples$subjects[5, 1, ])
plot(fm10$samples$subjects[6, 1, ])
plot(fm10$samples$subjects[7, 1, ])
plot(fm10$samples$subjects[8, 1, ])

plot(1 / fm10$samples$precision$pop)
plot(1 / fm10$samples$precision$sub)
plot(fm10$samples$precision$eps)

plot_spline(fm10)
plot_spline(truncate_spline(fm10, 10000))

system.time(fm11 <- SubjectsTpfMul(growth, 8, deg = 2, shape = "increasing", size = 10000, burn = 0))
plot(fm11$samples$population$hgtf[1, ])
plot(fm11$samples$population$hgtf[2, ])
plot(fm11$samples$population$hgtf[3, ])
plot(fm11$samples$population$hgtf[4, ])
plot(fm11$samples$population$hgtf[5, ])
plot(fm11$samples$population$hgtf[6, ])
plot(fm11$samples$population$hgtf[7, ])
plot(fm11$samples$population$hgtf[8, ])


plot(fm11$samples$subjects$hgtf[1, 1, ])
plot(fm11$samples$subjects$hgtf[2, 1, ])
plot(fm11$samples$subjects$hgtf[3, 1, ])
plot(fm11$samples$subjects$hgtf[4, 1, ])
plot(fm11$samples$subjects$hgtf[5, 1, ])
plot(fm11$samples$subjects$hgtf[6, 1, ])
plot(fm11$samples$subjects$hgtf[7, 1, ])
plot(fm11$samples$subjects$hgtf[8, 1, ])


plot_spline(fm11)
fm11t <- truncate_spline(fm11, 1500)
plot(1 / fm10t$samples$precision$pop)
plot(fm10t$samples$subject$hgtm[3, "boy01", ])

plot_spline(fm11t)
get_array_invs <- function(ary) {
    array(apply(ary, 3, solve), dim(ary))
}

## models without random effects on spline coefficients
source("pop-tpf.R")
fmpop <- pop_tpf(growth10, 8, deg = 2, shape = "increasing", size = 10000, burn = 0, verbose = T)
plot_spline(fmpop)


