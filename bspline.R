##### B SPLINES #####

library(TruncatedNormal)
library(nlme)
library(microbenchmark)
rm(list = ls())
setwd("~/Gdrive/master/algo/")
source("main.R")
source("subor.R")
sitka <- read.table("data/tsitka.txt", header = T)
##sitka <- read.table("data/sitka.txt", header = T)

## Initialise dataset and factors
y <- sitka$log.size
## y[length(y)] <- 1                       # Introduce an outlier
time <- sitka$days / 674
id.num <- sitka$id.num

dat <- data.frame(x = time, y, grp = id.num)

fm1 <- BasisLin(dat, 5)
lapply(fm1, rowMeans)
PlotLinSpline(


## Construct the knots, fixed and random effects design matrix (linear spline)

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

