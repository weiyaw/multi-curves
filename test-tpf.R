library(nlme)
rm(list = ls())
setwd("~/Dropbox/master/algo/")
source("main-tpf.R")
source("subor.R")
sitka <- read.table("data/sitka10.txt", header = T)
##sitka <- read.table("data/sitka.txt", header = T)

## Initialise dataset and factors
y <- sitka$log.size
## y[length(y)] <- 1                       # Introduce an outlier
id.num <- sitka$id.num
n <- NROW(sitka)
K <- 5
n.subject <- length(unique(sitka$id.num))
n.per.sub <- tapply(sitka$id.num, sitka$id.num, length)
pop.level <- factor(rep(1, n))
sub.level <- factor(sitka$id.num)


## Construct the knots, fixed and random effects design matrix (linear spline)
time.l <- sitka$days / 674
knots.l <- quantile(unique(time.l), seq(0, 1, length = K + 2))[-c(1, K + 2)]
X.l <- model.matrix(y ~ time.l, sitka)
Z.l <- outer(time.l, knots.l, `-`)
Z.l <- Z.l * (Z.l > 0)

## Fit a lme model (linear)
pop.pd.l <- pdIdent(~ Z.l - 1)
sub.pd.l <- pdBlocked(list(pdSymm(~ time.l), pdIdent(~ Z.l - 1)))
lin.fm <- lme(fixed = y ~ time.l,
              random = list(pop.level = pop.pd.l, sub.level = sub.pd.l))

## Construct the knots, fixed and random effects design matrix (linear spline)
time.q <- sitka$days / 67.4
knots.q <- quantile(unique(time.q), seq(0, 1, length = K + 2))[-c(1, K + 2)]
X.q <- model.matrix(y ~ time.q + I(time.q^2), sitka)
Z.q <- outer(time.q, knots.q, `-`)
Z.q <- Z.q^2 * (Z.q > 0)

## Fit a lme model (quadratic)
pop.pd.q <- pdIdent(~ Z.q - 1)
## sub.pd.q <- pdBlocked(list(pdSymm(~ time.q + I(time.q^2)), pdIdent(~ Z.q - 1)))
## sub.pd.q <- pdSymm(~ time.q + I(time.q^2))
sub.pd.q <- pdBlocked(list(pdSymm(~ 1), pdIdent(~ Z.q - 1)))
quad.fm <- lme(fixed = y ~ time.q + I(time.q^2),
               random = list(pop.level = pop.pd.q, sub.level = sub.pd.q))

## TEST LINEAR SPLINE

lin.sample1 <- SubjectLin(y, X.l, Z.l, lin.fm, 100, 10)
lapply(lin.sample1, rowMeans)
## lin.sample2 <- SubjectLin(y, X.l, Z.l, lin.fm, 10000, 100)

## grps.fact <- factor(lin.fm$groups$sub.level,
##                     levels = unique(lin.fm$groups$sub.level))
## PlotLinMean(lin.sample1, knots.l, range(time.l),
##             data.frame(time.l, y, grps.fact))

## PlotLinMean(lin.sample2, knots.l, range(time.l),
##             data.frame(time.l, y, grps.fact))

## uncon.set <- cbind(t(coef(lin.fm)), t(coef(lin.fm, level = 1)))
## colnames(uncon.set) <- with(lin.fm$groups, c(levels(sub.level), "population"))
## PlotLinSpline(uncon.set, knots.l, range(time.l),
##               data.frame(time.l, y, lin.fm$groups$sub.level))


## TEST QUADRATIC SPLINE

post <- SubjectQuad(y, X.q, Z.q, quad.fm, knots.q, range(time.q), 5000, 1)

source("subor.R")
quad.co <- t(coef(quad.fm, level = 2))
PlotQuadSpline(quad.co, knots.q, range(time.q),
               data.frame(time.q, y, quad.fm$groups$sub.level))



## TEST MIXED MODEL LINEAR SPLINE (SubjectsTpf)
rm(list = ls())
sitka <- read.table("data/sitka5.txt", header = T)
sitka <- with(sitka, data.frame(x = days / 674,
                                y = log.size,
                                grps = id.num))
source("main-tpf.R")
system.time(fm2 <- get_tpf_old(sitka, 5, 2, size = 100, burn = 0))
saveRDS(fm2, "tpf-long.rds")
source("graphs.R")
fm1 <- readRDS("simulations/tpf-lin.rds")
fm2 <- readRDS("simulations/tpf-quad.rds")
PlotSpline(fm1, range(sitka$x), sitka)
PlotSpline(fm2, range(sitka$x), sitka)

fm3 <- readRDS("simulations/bspline-lin.rds")
fm4 <- readRDS("simulations/bspline-quad.rds")
PlotSpline(fm3, range(sitka$x), sitka)
PlotSpline(fm4, range(sitka$x), sitka)

## TEST MIXED MODEL LINEAR SPLINE WITH MULTIPLE POPULATION (SubjectsTpfMul)
rm(list = ls())
sitka10 <- read.table("data/sitka10.txt", header = T)
sitka10 <- with(sitka10, data.frame(x = days / 674,
                                    y = log.size,
                                    grp.sub = id.num,
                                    grp.pop = ozone))
sitka10$grp.sub <- factor(sitka10$grp.sub,
                          levels = c("1", "2", "3", "4", "5",
                                     "60", "59", "56", "57", "58"))
source("main-tpf.R")
## single population
system.time(fm1 <- SubjectsTpf(sitka10, 5, deg = 2, shape = "increasing", size = 10))

system.time(fm2 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 100, burn = 0))
## multiple population (1000 burn, independent start)
fm2 <- readRDS("simulations/multi/multi-1k-bugless.rds")
## multiple population (1000 burn, independent start, garbage)
fm3 <- readRDS("simulations/multi/multi-10k.rds")
source("graphs.R")
plot_spline(fm1)
plot_spline(fm2, c(0, 1))
plot_spline(fm3, c(0, 1))
plot_spline(fm4, c(0, 1))

## multiple population (1000 burn, previous start)
system.time(fm4 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 10000, burn = 0))

source("main-tpf.R")

system.time(fm4 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 1000, burn = 0))

fm5 <- readRDS("simulations/multi/multi-sitka.rds")

system.time(fm5 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 1000, burn = 0))

source("graphs.R")
fm7t <- truncate(fm7, 2000)
plot_spline(fm7t)

plot(fm5$samples$population$`1`[1, ])
plot(fm5$samples$population$`1`[2, ])
plot(fm5$samples$population$`1`[3, ])
plot(fm5$samples$population$`1`[4, ])
plot(fm5$samples$population$`1`[5, ])
plot(fm5$samples$population$`1`[6, ])
plot(fm5$samples$population$`1`[7, ])

plot(fm5$samples$subject$`1`[1, "5", ])
plot(fm5$samples$subject$`1`[2, "5", ])
plot(fm5$samples$subject$`1`[3, "5", ])
plot(fm5$samples$subject$`1`[4, "5", ])
plot(fm5$samples$subject$`1`[5, "5", ])
plot(fm5$samples$subject$`1`[6, "5", ])
plot(fm5$samples$subject$`1`[7, "5", ])
plot_spline(fm5)

plot(fm6$samples$precision$pop)
plot(fm6$samples$precision$sub)
plot(fm6$samples$precision$poly[1, 1, 2000:10000])
plot(fm6$samples$precision$poly[2, 2, 2000:10000])
plot(fm6$samples$precision$poly[3, 3, 2000:10000])
plot(fm6$samples$precision$poly[1, 2, 2000:10000])
plot(fm6$samples$precision$poly[1, 3, 2000:10000])

source("graphs.R")
fm2t <- truncate(fm2, 300)
plot_spline(fm2t)


## Berkeley growth dataset
fm6 <- readRDS("simulations/multi/multi-growth.rds")
fm7 <- readRDS("simulations/multi/multi-growth-uncon.rds")

plot(fm7$samples$population$hgtf[1, ])
plot(fm7$samples$population$hgtf[2, ])
plot(fm7$samples$population$hgtf[3, ])
plot(fm7$samples$population$hgtf[4, ])
plot(fm7$samples$population$hgtf[5, ])
plot(fm7$samples$population$hgtf[6, ])
plot(fm7$samples$population$hgtf[7, ])


plot(fm6$samples$subject$hgtm[1, "boy01", ])
plot(fm6$samples$subject$hgtm[2, "boy01", ])
plot(fm6$samples$subject$hgtm[3, "boy01", ])
plot(fm6$samples$subject$hgtm[4, "boy01", ])
plot(fm6$samples$subject$hgtm[5, "boy01", ])
plot(fm6$samples$subject$hgtm[6, "boy01", ])
plot(fm6$samples$subject$hgtm[7, "boy01", ])

