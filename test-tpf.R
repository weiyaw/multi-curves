library(TruncatedNormal)
library(nlme)
library(microbenchmark)
rm(list = ls())
setwd("~/Gdrive/master/algo/")
source("main-tpf.R")
source("subor.R")
sitka <- read.table("data/tsitka.txt", header = T)
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
sitka <- read.table("data/tsitka.txt", header = T)
sitka <- with(sitka, data.frame(x = days / 674,
                                y = log.size,
                                grps = id.num))
source("main-tpf.R")
system.time(fm2 <- SubjectsTpf(sitka, 5, size = 100))
saveRDS(fm2, "tpf-long.rds")
source("graphs.R")
PlotSpline(fm1, range(sitka$x), sitka)






