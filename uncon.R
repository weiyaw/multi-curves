rm(list = ls())
library(nlme)
library(ggplot2)
sitka <- read.table("data/tsitka.txt", header = T)
ggplot(sitka, aes(days, log.size, group = id.num)) + geom_line(aes(col = id.num))

y <- sitka$log.size
time <- sitka$days
subject <- sitka$id.num
nSubject <- length(unique(subject))
## Set up knots
K <- max(5, min(floor(length(unique(time)) / 4), 40))
knots <- quantile(unique(time), seq(0, 1, length = K + 2))[-c(1, K + 2)]

## Set up design matrices
X <- model.matrix(y ~ time)
Z <- outer(time, knots, "-")
Z <- Z * (Z > 0)
n <- length(y)

popSpline <- rep(1, n)
subLinear <- subject
fit1 <- lme(y ~ time, random = list(popSpline = pdIdent(~ Z - 1),
                                    subLinear = pdIdent(~ 1)))
summary(fit1)

fit2 <- lme(y ~ time, random = list(popSpline = pdIdent(~ Z - 1),
                                    subLinear = pdSymm(~ time)))
summary(fit2)

## Subject specific
scalTime <- time / 674
K <- 5
kSubject <- K
knots <- quantile(unique(scalTime), seq(0, 1, length = K + 2))[-c(1, K + 2)]

X <- model.matrix(y ~ scalTime)
Z <- outer(scalTime, knots, "-")
Z <- Z * (Z > 0)
n <- length(y)

knotsSubject <- quantile(unique(scalTime), seq(0, 1, length = kSubject + 2)
                         )[-c(1, kSubject + 2)]
zSubject <- outer(scalTime, knotsSubject, "-")
zSubject <- zSubject * (zSubject > 0)

pop.level <- factor(rep(1, n))
sub.level <- factor(subject)
fit3 <- lme(fixed = y ~ scalTime,
            random = list(pop.level = pdIdent(~ Z - 1),
                          sub.level = pdBlocked(list(pdSymm(~ scalTime),
                                                     pdIdent(~ zSubject - 1)))))

## Subject specific (same population and individuals splines)
fit4 <- lme(fixed = y ~ scalTime,
            random = list(pop.level = pdIdent(~ Z - 1),
                          sub.level = pdBlocked(list(pdSymm(~ scalTime),
                                                     pdIdent(~ Z - 1)))))

## Calculate prediction
PredictLme <- function(X, beta, Z, u) {
    X <- as.matrix(X)
    beta <- as.matrix(beta)
    Z <- as.matrix(Z)
    u <- as.matrix(u)
    y <- X %*% beta + Z %*% u
    colnames(y) <- rownames(y) <- NULL
    return(y)
}

# Calculate subject specific prediction
new.beta <- t(coef(fit4)[, 1:2])
new.u <- t(coef(fit4)[, -(1:2)])
new.scalTime <- c(min(scalTime), knots, max(scalTime))
new.X <- cbind(1, new.scalTime, deparse.level = 0)
new.Z <- outer(new.scalTime, knots, `-`)
new.Z <- new.Z * (new.Z > 0)
new.y.mat <- PredictLme(new.X, new.beta, new.Z, new.u) # each columm corresponds to prediction of each subject
new.data <- data.frame(log.size = as.vector(new.y.mat),
                       scalTime = rep(new.scalTime, nSubject),
                       id.num = rep(seq(1, nSubject), each = (K + 2)))

# Calculate population prediction
new.pop.beta <- t(coef(fit4, level = 1)[1:2])
new.pop.u <- t(coef(fit4, level = 1)[-(1:2)])
new.pop.y.mat <- PredictLme(new.X, new.pop.beta, new.Z, new.pop.u)
new.pop.data <- data.frame(log.size = as.vector(new.pop.y.mat),
                           scalTime = new.scalTime)

ggplot(mapping = aes(x = scalTime, y = log.size)) +
    geom_point(aes(col = id.num), sitka, size = 2) +
    geom_line(aes(group = id.num, col = id.num), new.data) +
    geom_line(data = new.pop.data, col = 'red')

dims <- rev(fit4$dims$ncol)[-(1:2)]
pop.varcov.idx <- seq(1, length.out = dims["pop.level"] + 1)
sub.varcov.idx <- seq(dims["pop.level"] + 2, length.out = dims["sub.level"] + 1)
pop.varcov <- VarCorr(fit4)[pop.varcov.idx, ]
sub.varcov <- VarCorr(fit4)[sub.varcov.idx, ]


?rmvtgauss.lin
