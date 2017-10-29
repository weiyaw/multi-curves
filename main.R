library(TruncatedNormal)
library(nlme)
library(microbenchmark)
library(mvtnorm)
rm(list = ls())
setwd("~/Gdrive/master/algo/")
source("subor.R")
sitka <- read.table("data/tsitka.txt", header = T)

## Initialise dataset and factors
y <- sitka$log.size
y[length(y)] <- 4
scal.time <- sitka$days / 674
id.num <- sitka$id.num
n <- NROW(sitka)
K <- 5
n.subject <- length(unique(sitka$id.num))
n.per.sub <- tapply(sitka$id.num, sitka$id.num, length)
n.terms <- K + 2
pop.level <- factor(rep(1, n))
sub.level <- factor(sitka$id.num)

## Construct the knots, fixed and random effects design matrix
knots <- quantile(unique(scal.time), seq(0, 1, length = K + 2))[-c(1, K + 2)]
X <- model.matrix(y ~ scal.time, sitka)
Z <- outer(scal.time, knots, `-`)
Z <- Z * (Z > 0)

## Fit a lme model
fit1 <- lme(fixed = y ~ scal.time,
            random = list(pop.level = pdIdent(~ Z - 1),
                          sub.level = pdBlocked(list(pdSymm(~ scal.time),
                                                     pdIdent(~ Z - 1)))))

fit2 <- lme(fixed = y ~ scal.time,
            random = list(pop.level = pdIdent(~ Z - 1),
                          sub.level = pdBlocked(list(pdLogChol(~ scal.time),
                                                     pdIdent(~ Z - 1)))))
## Extract EBLUPS
pop.coef <- as.matrix(coef(fit1, level = 1))
sub.coef <- as.matrix(coef(fit1, level = 2)) -
    matrix(rep(pop.coef, n.subject), n.subject, byrow = TRUE)
crude.zeta <- c(pop.coef, t(sub.coef))

## Extract variance covariance matrix
varcov.list <- MultiVarCov(fit1$modelStruct$reStruct, fit1$sigma^2)
pop.varcov.list <- list(diag(1, 2), varcov.list$pop.level)
sub.varcov.list <- rep(varcov.list["sub.level"], n.subject)
crude.varcov <- DiagMat(c(pop.varcov.list, sub.varcov.list))

## Construct the constraint matrix A
n.fixef <- fit1$dims$ncol[3]
n.sub.ranef <- fit1$dims$ncol[1]
single.A <- matrix(1, K + 1, K + 1)
single.A[upper.tri(single.A)] <- 0
single.A <- cbind(0, single.A)


A <- DiagMat(single.A, n.subject)
A <- rbind(matrix(0, NROW(single.A), NCOL(A)), A)
A <- cbind(do.call("rbind", rep(list(single.A), n.subject + 1)), A)

## rmvtgauss.lin(500, crude.zeta, crude.varcov, Amat = t(A), Avec = rep(0, NROW(A)))

zeta.dim <- NCOL(crude.zeta)


ApproxGibbs <- function(y, linear.mat, spline.mat, lme.obj, size = 1000) {
    ## y : response
    ## linear.mat: model matrix for straight line
    ## spline.mat: model matrix for splines
    ## lme.obj: unconstrained lme object from nlme
    ## size: number of samples from the posterior distribution

    ## Initialise
    varcov.list <- lapply(lme.obj$modelStruct$reStruct, as.matrix)
    varcov.list <- lapply(varcov.list, function (x) { x * lme.obj$sigma^2 })
    n.per.sub <- tapply(lme.obj$groups[, -1], lme.obj$groups[, -1], length)
    n.subject <- lme.obj$dims$ngrps[1]
    n.fixed <- lme.obj$dims$ncol[length(lme.obj$dims$ncol) - 1]
    n.terms <- lme.obj$dims$qvec[1]
    n.splines <- n.terms - n.fixed

    ## EBLUPS
    pop.coef <- as.matrix(coef(lme.obj, level = 1))
    sub.coef <- as.matrix(coef(lme.obj, level = 2)) -
        VecToMat(pop.coef, n.subject, FALSE)
    crude.mean <- c(pop.coef, t(sub.coef))

    ## REML variance covariance
    pop.varcov <- list(diag(1, n.fixed), varcov.list[[2]])
    sub.varcov <- rep(varcov.list[1], n.subject)
    varcov <- DiagMat(c(pop.varcov, sub.varcov))

    ## Constraint matrix
    single.A <- matrix(1, n.splines + 1, n.splines + 1)
    single.A[upper.tri(single.A)] <- 0
    single.A <- cbind(0, single.A)

    A <- DiagMat(single.A, n.subject)
    A <- rbind(matrix(0, NROW(single.A), NCOL(A)), A)
    A <- cbind(do.call("rbind", rep(list(single.A), n.subject + 1)), A)

    ## Candidates from independent proposal
    cand <- rmvtgauss.lin(size, crude.mean, varcov, Amat = t(A),
                          Avec = rep(0, NROW(A)))
    unif.rv <- runif(size, 0, 1)

    ## Model matrix
    model.mat <- cbind(linear.mat, spline.mat, deparse.level = 0)

    ## Initialise density functions
    Prop <- ProposalLogFac(crude.mean, varcov, n.subject, n.terms)
    Like <- LikelihoodLogFac(y, model.mat, n.terms, n.per.sub)
    Prior <- PriorLogFac(n.terms, n.subject)

    ## Fixed variance of random components
    eps.prec <- 1 / lme.obj$sigma^2
    u.prec <- 1 / pop.varcov[[2]][1, 1]
    v.prec <- 1 / sub.varcov[[1]][3, 3]
    b.prec <- solve(sub.varcov[[1]][1:2, 1:2])

    ## MH step
    zeta.prev <- crude.mean
    prop.prev <- Prop(crude.mean)
    like.prev <- Like(crude.mean, eps.prec)
    prior.prev <- Prior(crude.mean, u.prec, v.prec, b.prec)

    res <- matrix(NA, length(crude.mean), size)

    for (i in seq_len(size)) {
        zeta.curr <- cand[, i]
        prop.curr <- Prop(zeta.curr)
        like.curr <- Like(zeta.curr, eps.prec)
        prior.curr <- Prior(zeta.curr, u.prec, v.prec, b.prec)
        browser()
        accpt.prob <- exp((prop.prev + like.curr + prior.curr) -
                          (prop.curr + like.prev + prior.prev))

        if (unif.rv[i] < accpt.prob) {
            res[, i] <- zeta.curr
            zeta.prev <- zeta.curr
            like.prev <- like.curr
            prior.prev <- prior.curr
        } else {
            res[, i] <- zeta.prev
        }
        if (i %% 100 == 0) {
            cat(i, " ")
        }
    }
    cat("\n")
    return(res)
}

post <- ApproxGibbs(y, X, Z, fit2, 5)
hist(post[3,])


ExactGibbs <- function(y, linear.mat, spline.mat, lme.obj, size = 100,
                       burn = size / 10) {
    ## y : response (better be sorted on subjects indices)
    ## linear.mat: model matrix for straight line
    ## spline.mat: model matrix for splines
    ## lme.obj: unconstrained lme object from nlme
    ## size: number of samples from the posterior distribution

    ## EBLUPS
    pop.coef <- as.matrix(coef(lme.obj, level = 1))

    ## REML variance covariance
    var.error <- lme.obj$sigma^2
    varcov.list <- lapply(lme.obj$modelStruct$reStruct, as.matrix)
    varcov.list <- lapply(varcov.list, function (x) { x * var.error })
    var.pop <- varcov.list[[2]][1, 1]
    varcov.sub <- varcov.list[[1]]

    ## Initialise various lengths
    n <- length(y)
    n.subject <- lme.obj$dims$ngrps[1]
    n.fixed <- lme.obj$dims$ncol[length(lme.obj$dims$ncol) - 1]
    n.terms <- lme.obj$dims$qvec[1]
    n.splines <- n.terms - n.fixed

    ## Population and subject share the same constraint matrix
    constrt.mat <- matrix(1, n.splines + 1, n.splines + 1)
    constrt.mat[lower.tri(constrt.mat)] <- 0
    constrt.mat <- rbind(0, constrt.mat)

    ## Indices in the dataset for each subject
    idx.sub <- tapply(seq_len(n), lme.obj$groups[, -1], function(x) x)

    ## Model matrix
    model.mat <- cbind(linear.mat, spline.mat)

    ## Variance of the conditional SUBJECT posterior
    ## Factor in front of the mean of the conditional SUBJECT posterior
    M.inv.sub <- rep(list(solve(varcov.sub) * var.error), n.subject)
    cond.fact.sub <- list()             # factor of the conditional mean
    cond.var.sub <- list()              # varcov of the conditional dist.
    for (i in seq_len(n.subject)) {
        M.inv.sub[[i]] <- solve(crossprod(model.mat[idx.sub[[i]], ]) +
                                M.inv.sub[[i]])
        cond.fact.sub[[i]] <- tcrossprod(M.inv.sub[[i]],
                                         model.mat[idx.sub[[i]], ])
        cond.var.sub[[i]] <- var.error * M.inv.sub[[i]]
    }

    ## Variance of the conditional POPULATION posterior
    ## Factor in front of the mean of the conditional POPULATION posterior
    M.inv.pop <- diag(c(rep(0, n.fixed), rep(var.error / var.pop, n.splines)))
    M.inv.pop <- solve(crossprod(model.mat) + M.inv.pop)
    cond.fact.pop <- tcrossprod(M.inv.pop, model.mat) # factor of the cond. mean
    cond.var.pop <- var.error * M.inv.pop       # varcov of the cond. dist.

    ## Initialise current estimates, EBLUPS as initial population response curve
    curr.pop <- as.vector(pop.coef)
    curr.sub <- matrix(NA, n.terms, n.subject)

    ## Record current individual curves contribution to the prediction.
    ## Essentially: model.mat %*% coef.sub
    curr.pred.sub <- rep(NA, n)

    if (burn < 1) {
        stop("Must burn at least 1 sample.")
    }

    ## Sequence of the subject (indices for the "for" loop)
    seq.subject <- seq_len(n.subject)

    ## Burning period. Results are discarded
    for (i in seq_len(burn)) {

        ## Simulate subject posterior (conditioned on previous population)
        ## Lower-bound vector for the subject posterior
        avec.sub <- -1 * cumsum(curr.pop[-1])

        ## Generate an individual curve for each subject
        for (j in seq.subject) {
            idx <- idx.sub[[j]]
            model.mat.sub <- model.mat[idx, , drop = FALSE]
            mu.sub <- cond.fact.sub[[j]] %*%
                (y[idx] -  model.mat.sub %*% curr.pop)
            curr.sub[, j] <- rmvtgauss.lin(1, mu.sub, cond.var.sub[[j]],
                                           Amat = constrt.mat,
                                           Avec = avec.sub)
            curr.pred.sub[idx] <- model.mat.sub %*% curr.sub[, j]
        }

        if (i >  burn) {
            break
        }

        ## Simulated population posterior (conditioned on current subject)
        ## Lower-bound vector for the population posterior
        avec.pop <- apply(curr.sub[-1, , drop = FALSE], 2, cumsum)
        avec.pop <- -1 * apply(avec.pop, 1, min)
        avec.pop[avec.pop < 0] <- 0

        ## Generate a population curve
        mu.pop <- cond.fact.pop %*% (y - curr.pred.sub)
        start.gauss <- c(mu.pop[1], avec.pop[1] + 1, diff(avec.pop))
        curr.pop <- new.rmvtgauss.lin(1, mu.pop, cond.var.pop,
                                      Amat = constrt.mat,
                                      Avec = avec.pop,
                                      start = start.gauss,
                                      burnin = 50)
    }

    ## Initialise the output list
    ## "res[[i]]" to access i th subject curve
    ## "res$pop" to access population curve
    res <- rep(list(matrix(NA, n.terms, size)), n.subject)
    res$pop <- matrix(NA, n.terms, size)

    ## Generate posterior. Results are recorded
    for (i in seq_len(size)) {

        ## Simulated population posterior (conditioned on current subject)
        ## Lower-bound vector for the population posterior
        avec.pop <- apply(curr.sub[-1, , drop = FALSE], 2, cumsum)
        avec.pop <- -1 * apply(avec.pop, 1, min)
        avec.pop[avec.pop < 0] <- 0

        ## Generate a population curve
        mu.pop <- cond.fact.pop %*% (y - curr.pred.sub)
        start <- c(mu.pop[1], avec.pop[1] + 1, diff(avec.pop))
        curr.pop <- new.rmvtgauss.lin(1, mu.pop, cond.var.pop,
                                      Amat = constrt.mat,
                                      Avec = avec.pop,
                                      start = start,
                                      burnin = 50)
        res$pop[, i] <- curr.pop

        ## Simulate subject posterior (conditioned on previous population)
        ## Lower-bound vector for the subject posterior
        avec.sub <- -1 * cumsum(curr.pop[-1])

        ## Generate an individual curve for each subject
        for (j in seq.subject) {
            idx <- idx.sub[[j]]
            model.mat.sub <- model.mat[idx, , drop = FALSE]
            mu.sub <- cond.fact.sub[[j]] %*%
                (y[idx] -  model.mat.sub %*% curr.pop)
            curr.sub[, j] <- rmvtgauss.lin(1, mu.sub, cond.var.sub[[j]],
                                           Amat = constrt.mat,
                                           Avec = avec.sub)
            curr.pred.sub[idx] <- model.mat.sub %*% curr.sub[, j]
            res[[j]][, i] <- curr.sub[, j]
        }

        if (i %% 1000 == 0) {
            cat(i, " samples generated.\n")
        }
    }
    return(res)
}



linear.mat <- X
spline.mat <- Z
lme.obj <- fit2
post <- ExactGibbs(y, linear.mat, spline.mat, lme.obj, 100, 10)
post2 <- ExactGibbs(y, linear.mat, spline.mat, lme.obj, 10000, 100)



PlotLinSpline <- function(coefs, knots, limits, data) {

    if (is.vector(coefs)) {
        coefs <- as.matrix(coefs)
    } else if (!is.matrix(coefs)) {
        stop("coefs must be a matrix")
    }

    if (!is.vector(knots) || !is.vector(limits)) {
        stop("knots and limits must be vectors")
    }

    rownames(coefs) <- NULL
    names(knots) <- NULL
    time <- c(seq(min(limits), max(limits), length.out = 100), knots)
    time <- time[order(time)]
    basis <- outer(time, knots, `-`)
    basis <- basis * (basis > 0)
    model.mat <- cbind(1, time, basis, deparse.level = 0)
    y <- model.mat %*% coefs
    plot.data <- melt(y)
    plot.data$time <- rep(time, times = NCOL(coefs))

    if (missing(data)) {
        ggplot() +
            geom_line(aes(time, value, group = X2, col = factor(X2)), plot.data)
    } else if (is.data.frame(data)) {
        ggplot() +
            geom_line(aes(time, value, group = X2, col = factor(X2)), plot.data) +
            geom_point(aes(data[[1]], data[[2]], col = factor(data[[3]])))
    }
}

j <- 1000
one.set <- matrix(NA, n.terms, n.subject + 1)
colnames(one.set) <- with(fit2$groups, c(levels(sub.level), "population"))
one.set[, NCOL(one.set)] <- post2$pop[, j]
for (i in seq_len(n.subject)) {
    one.set[, i] <- post2$pop[, j] + post2[[i]][, j]
}
PlotLinSpline(one.set, knots, range(scal.time),
              data.frame(scal.time, y, fit2$groups$sub.level))

uncon.set <- cbind(t(coef(fit2)), t(coef(fit2, level = 1)))
colnames(uncon.set) <- with(fit2$groups, c(levels(sub.level), "population"))
PlotLinSpline(uncon.set, knots, range(scal.time),
              data.frame(scal.time, y, fit2$groups$sub.level))

hist(post[[1]][1, -(1:100)])

