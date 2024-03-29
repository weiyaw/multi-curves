##### B SPLINES #####

library(microbenchmark)
rm(list = ls())
setwd("~/Gdrive/master/algo/")
source("main-bspline.R")
source("subor.R")
source("graphs.R")
source("normgen.R")
sitka <- read.table("data/tsitka.txt", header = T)
## sitka <- sitka[sitka$id.num == 1 | sitka$id.num == 2, ]
## sitka <- read.table("data/sitka.txt", header = T)

## Initialise dataset and factors
y <- sitka$log.size
## y[length(y)] <- 1                       # Introduce an outlier
time <- sitka$days / 674
id.num <- sitka$id.num

dat <- data.frame(x = time, y, grp = factor(id.num))

source("main.R")
## system.time(fm1 <- SubjectsBsPar(dat, 5, size = 100, burn = 1))
system.time(fm1 <- SubjectsBs(dat, 5, size = 10, burn = 1))
PlotSpline(fm1, c(0, 1.5), dat, 500)

fm1.long <- readRDS("fm1-long.rds")
PlotSpline(fm1.long, c(0, 1.5), dat, 500)

source("normgen.R")
tmp <- tempfile()
Rprof(tmp)
for (ii in 1:1000) {
    sub <- tmvtnorm::rtmvnorm2(1, mean = as.vector(fm1[[1]]),
                               sigma = fm1[[2]],
                               D = fm1[[3]],
                               lower = fm1[[4]],
                               upper = rep(Inf, 6),
                               algorithm = "gibbs",
                               start.value = fm1[[5]],
                               burn.in = 100)}
Rprof(NULL)
summaryRprof(tmp)


sub1 <- TruncatedNormal::rmvtgauss.lin(1, arg$mu, arg$sig,
                                       Amat = t(arg[[4]]),
                                       Avec = -1 * arg[[5]],
                                       start = arg[[3]],
                                       burn = 1000)

fm1 <- SubjectsBsPar(dat, 5, size = 500, burn = 1)
rm(list = setdiff(ls(), c("dat", "fm1")))
PlotSpline(fm1, c(0, 1.5), dat, 500)
lapply(fm1$samples, rowMeans)

fm2 <- SubjectsBs(dat, 5, deg = 2)
PlotSpline(fm2, c(0, 1.5), dat, 500)


## Construct the knots, fixed and random effects design matrix (linear spline)
K <- 5
deg <- 1
EPS <- 1e-6

design <- BSplineDesign(time, K, deg)
B <- design$design
knots <- design$knots
rm(design)


pop.level <- rep(1, length(time))
sub.level <- id.num
g <- 1:NCOL(B)
G <- model.matrix(~ g)


## Fit a lme model (linear)
library(nlme)
source("graphs.R")
pop.pd.l <- pdIdent(~ B - 1)
sub.pd.l <- pdIdent(~ B - 1)
fnlme.1 <- lme(fixed = y ~ -1,
               random = list(pop.level = pop.pd.l, sub.level = sub.pd.l))
PlotLinSpline(t(coef(fnlme.1)[1, ]), knots, c(0, 1.5), dat, "bs")
fnlme.2 <- lme(fixed = y ~ -1,
               random = list(sub.level = sub.pd.l))
fnlme.3 <- lme(fixed = y ~ B - 1,
               random = list(sub.level = sub.pd.l))

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



#### BERKELY GROWTH DATASET ####
rm(list = ls())
growth <- reshape2::melt(fda::growth[-3])
growth <- with(growth, data.frame(x = Var1 / max(Var1), y = value / max(value), grp.sub = Var2, grp.pop = L1))
growth10 <- subset(growth,
                   grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                  "girl01", "girl02", "girl03", "girl04", "girl05"),
                   drop = TRUE)
growth10$grp.sub <- droplevels(growth10$grp.sub)

growth10boys <- subset(growth,
                       grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                      "boy06", "boy07", "boy08", "boy09", "boy10"),
                   drop = TRUE)
growth10boys$grp.sub <- droplevels(growth10boys$grp.sub)

data <- growth10boys
deg <- 2
data[, 1:2] <- apply(data[, 1:2], 2, function(x) x / max(x))

set.seed(1)
## set.seed(2)
source("main-bs.R")
system.time(fm8 <- sub_bs(data, 8, deg = 2, shape = "increasing", size = 10000, burn = 0, sampler = "gibbs", verbose = T))

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

source("graphs.R")
plot_spline(fm8, mle = F)
plot_spline(fm9, mle = T)
plot_spline(truncate_spline(fm8, 5000))
plot_spline(truncate_spline(fm9, 5000))






