SubjectLin <- function(y, linear.mat, spline.mat, lme.obj, size = 100,
                       burn = size / 10) {
    ## y : response used in lme.obj (better be sorted according to subjects
    ##     indices)
    ## linear.mat: model matrix for straight line
    ## spline.mat: model matrix for splines
    ## lme.obj: unconstrained lme object from nlme
    ## size: number of samples from the posterior distribution

    if (burn < 1) {
        stop("Must burn at least 1 sample.")
    }

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
    constrt <- matrix(1, n.splines + 1, n.splines + 1)
    constrt[lower.tri(constrt)] <- 0
    constrt <- rbind(0, constrt)
    constrt.inv <- solve(constrt[-1, , drop = FALSE])

    ## Indices of datapoints corresponding to each subject
    fact.sub <- factor(lme.obj$groups[, -1],
                       levels = unique(as.character(lme.obj$groups[, -1])))
    idx.sub <- tapply(seq_len(n), fact.sub, function(x) x)

    ## Model matrix
    model.mat <- cbind(linear.mat, spline.mat)

    ## Variance of the conditional SUBJECT posterior
    ## Factor in front of the mean of the conditional SUBJECT posterior
    M.sub <- rep(list(solve(varcov.sub) * var.error), n.subject)
    cond.fact.sub <- list()             # factor of the conditional mean
    cond.var.sub <- list()              # varcov of the conditional dist.
    for (i in seq_len(n.subject)) {
        M.sub[[i]] <- solve(crossprod(model.mat[idx.sub[[i]], ]) +
                                M.sub[[i]])
        cond.fact.sub[[i]] <- tcrossprod(M.sub[[i]],
                                         model.mat[idx.sub[[i]], ])
        cond.var.sub[[i]] <- var.error * M.sub[[i]]
    }

    ## Variance of the conditional POPULATION posterior
    ## Factor in front of the mean of the conditional POPULATION posterior
    M.pop <- diag(c(rep(0, n.fixed), rep(var.error / var.pop, n.splines)))
    M.pop <- solve(crossprod(model.mat) + M.pop)
    cond.fact.pop <- tcrossprod(M.pop, model.mat) # factor of the cond. mean
    cond.var.pop <- var.error * M.pop       # varcov of the cond. dist.

    ## Initialise current estimates, EBLUPS as initial population response curve
    ## If EBLUPS is not monotone, make it monotone whilst retaining its shape AMAP.
    ## grad.pop = gradients of the population curve

    curr.pop <- as.vector(pop.coef)
    grad.pop <- crossprod(constrt, curr.pop)
    grad.pop[grad.pop < 0.01] <- 0.01
    curr.pop <- c(curr.pop[1], crossprod(constrt.inv, grad.pop))

    curr.sub <- matrix(NA, n.terms, n.subject)

    ## Record current individual curves contribution to the prediction.
    ## Essentially: model.mat %*% coef.sub
    curr.pred.sub <- rep(NA, n)

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
                                           Amat = constrt,
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
                                      Amat = constrt,
                                      Avec = avec.pop,
                                      start = start.gauss,
                                      burnin = 50)

        if (i %% 1000 == 0) {
            cat(i, " samples burned.\n")
        }

    }

    ## Initialise the output list
    ## "samples[[i]]" to access i th subject curve
    ## "samples$pop" to access population curve
    samples <- rep(list(matrix(NA, n.terms, size)), n.subject)
    names(samples) <- names(idx.sub)
    samples$population <- matrix(NA, n.terms, size)

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
                                      Amat = constrt,
                                      Avec = avec.pop,
                                      start = start,
                                      burnin = 50)
        samples$population[, i] <- curr.pop

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
                                           Amat = constrt,
                                           Avec = avec.sub)
            curr.pred.sub[idx] <- model.mat.sub %*% curr.sub[, j]
            samples[[j]][, i] <- curr.sub[, j]
        }

        if (i %% 1000 == 0) {
            cat(i, " samples generated.\n")
        }
    }
    means <- lapply(samples, rowMeans)
    basis <- list(type = "tpf", knots = NA, degree = 1)
    res <- list(means = means, samples = samples, basis = basis)
    return(res)
}





SubjectQuad <- function(y, quad.mat, spline.mat, lme.obj, knots, limits,
                        size = 100, burn = size / 10) {
    ## y : response used in lme.obj (better be sorted according to subjects
    ##     indices)
    ## quad.mat: model matrix for quadratic polynomial
    ## spline.mat: model matrix for quadratic splines
    ## lme.obj: unconstrained lme object from nlme
    ## knots : the knots locations
    ## limits: the range on which the monotonicity constraint is applied
    ## size: number of samples from the posterior distribution

    if (burn < 1) {
        stop("Must burn at least 1 sample.")
    }

    ## EBLUPS
    pop.coef <- as.matrix(coef(lme.obj, level = 1))
    kappa <- c(min(limits), knots, max(limits))
    names(kappa) <- NULL

    ## REML variance covariance
    var.error <- lme.obj$sigma^2
    varcov.list <- lapply(lme.obj$modelStruct$reStruct, as.matrix)
    varcov.list <- lapply(varcov.list, function (x) { x * var.error })
    var.pop <- varcov.list[[2]][1, 1]
    varcov.sub <- varcov.list[[1]]

    ## Initialise various lengths
    n <- length(y)
    n.subject <- lme.obj$dims$ngrps[1]
    n.terms <- lme.obj$dims$qvec[1]
    n.fixed <- lme.obj$dims$ncol[3]
    n.splines <- lme.obj$dims$qvec[2]
    n.knots <- length(knots)

    if ((n.fixed + n.splines) != n.terms) {
        stop("Number of fixed or spline terms incorrect.")
    }

    ## Population and subject share the same constraint matrix
    constrt <- -2 * outer(knots, c(knots[-1], max(limits)), `-`)
    constrt[lower.tri(constrt)] <- 0
    constrt <- cbind(0, 0, constrt, deparse.level = 0)
    constrt <- rbind(0, 1, 2 * kappa, constrt, deparse.level = 0)
    constrt.inv <- solve(constrt[-1, , drop = FALSE])

    colnames(constrt) <- NULL
    rownames(constrt) <- NULL

    ## Indices of datapoints corresponding to each subject
    fact.sub <- factor(lme.obj$groups[, -1],
                       levels = unique(as.character(lme.obj$groups[, -1])))
    idx.sub <- tapply(seq_len(n), fact.sub, function(x) x)

    ## Model matrix
    model.mat <- cbind(quad.mat, spline.mat)

    ## Variance of the conditional SUBJECT posterior
    ## Factor in front of the mean of the conditional SUBJECT posterior
    M.sub <- rep(list(solve(varcov.sub) * var.error), n.subject)
    cond.fact.sub <- list()             # factor of the conditional mean
    cond.var.sub <- list()              # varcov of the conditional dist.
    for (i in seq_len(n.subject)) {
        M.sub[[i]] <- solve(crossprod(model.mat[idx.sub[[i]], ]) +
                                M.sub[[i]])
        cond.fact.sub[[i]] <- tcrossprod(M.sub[[i]],
                                         model.mat[idx.sub[[i]], ])
        cond.var.sub[[i]] <- var.error * M.sub[[i]]
    }

    ## Variance of the conditional POPULATION posterior
    ## Factor in front of the mean of the conditional POPULATION posterior
    M.pop <- diag(c(rep(0, n.fixed), rep(var.error / var.pop, n.splines)))
    M.pop <- solve(crossprod(model.mat) + M.pop)
    cond.fact.pop <- tcrossprod(M.pop, model.mat) # factor of the cond. mean
    cond.var.pop <- var.error * M.pop       # varcov of the cond. dist.

    ## Initialise current estimates, EBLUPS as initial population response curve
    ## If EBLUPS is not monotone, make it monotone whilst retaining its shape AMAP.
    curr.pop <- as.vector(pop.coef)
    grad.pop <- crossprod(constrt, curr.pop)
    grad.pop[grad.pop < 0.01] <- 0.01
    curr.pop <- c(curr.pop[1], crossprod(constrt.inv, grad.pop))

    curr.sub <- matrix(NA, n.terms, n.subject)

    ## Record current individual curves contribution to the prediction.
    ## Essentially: model.mat %*% coef.sub
    curr.pred.sub <- rep(NA, n)

    ## Sequence of the subject (indices for the "for" loop)
    seq.subject <- seq_len(n.subject)

    ## Starting values for rmvtgauss.lin when generating SUBJECTS
    start.sub <- rep(0, n.terms)

    ## Burning period. Results are discarded
    for (k in seq_len(burn)) {

        ## Simulate subject posterior (conditioned on previous population)
        ## Lower-bound vector for the subject posterior
        lower.sub <- -1 * crossprod(constrt, curr.pop)

        ## Generate an individual curve for each subject
        for (i in seq.subject) {
            idx <- idx.sub[[i]]
            X.i <- model.mat[idx, , drop = FALSE]
            mu.sub <- cond.fact.sub[[i]] %*% (y[idx] -  X.i %*% curr.pop)
            curr.sub[, i] <- rmvtgauss.lin(1, mu.sub, cond.var.sub[[i]],
                                           Amat = constrt,
                                           Avec = lower.sub,
                                           start = start.sub,
                                           burnin = 50)
            curr.pred.sub[idx] <- X.i %*% curr.sub[, i]
        }


        if (k >  burn) {
            break
        }

        ## Simulated population posterior (conditioned on current subject)
        ## Lower-bound vector for the population posterior
        lower.pop <- crossprod(constrt, curr.sub)
        lower.pop <- -1 * apply(lower.pop, 1, min)
        lower.pop[lower.pop < 0] <- 0

        ## Generate a population curve
        mu.pop <- cond.fact.pop %*% (y - curr.pred.sub)
        start.pop <- c(mu.pop[1], crossprod(constrt.inv, lower.pop + 0.01))
        curr.pop <- new.rmvtgauss.lin(1, mu.pop, cond.var.pop,
                                      Amat = constrt,
                                      Avec = lower.pop,
                                      start = start.pop,
                                      burnin = 50)

        if (k %% 1000 == 0) {
            cat(k, " samples burned.\n")
        }

    }

    ## Initialise the output list
    ## "res[[i]]" to access i th subject curve
    ## "res$pop" to access population curve
    res <- rep(list(matrix(NA, n.terms, size)), n.subject)
    names(res) <- names(idx.sub)
    res$population <- matrix(NA, n.terms, size)

    ## Generate posterior. Results are recorded
    for (k in seq_len(size)) {

        ## Simulated population posterior (conditioned on current subject)
        ## Lower-bound vector for the population posterior
        lower.pop <- crossprod(constrt, curr.sub)
        lower.pop <- -1 * apply(lower.pop, 1, min)
        lower.pop[lower.pop < 0] <- 0

        ## Generate a population curve
        mu.pop <- cond.fact.pop %*% (y - curr.pred.sub)
        start.pop <- c(mu.pop[1], crossprod(constrt.inv, lower.pop + 0.01))
        curr.pop <- new.rmvtgauss.lin(1, mu.pop, cond.var.pop,
                                      Amat = constrt,
                                      Avec = lower.pop,
                                      start = start.pop,
                                      burnin = 50)
        res$population[, k] <- curr.pop

        ## Simulate subject posterior (conditioned on previous population)
        ## Lower-bound vector for the subject posterior
        lower.sub <- -1 * crossprod(constrt, curr.pop)

        ## Generate an individual curve for each subject
        for (i in seq.subject) {
            idx <- idx.sub[[i]]
            X.i <- model.mat[idx, , drop = FALSE]
            mu.sub <- cond.fact.sub[[i]] %*%
                (y[idx] -  X.i %*% curr.pop)
            curr.sub[, i] <- rmvtgauss.lin(1, mu.sub, cond.var.sub[[i]],
                                           Amat = constrt,
                                           Avec = lower.sub,
                                           start = start.sub,
                                           burnin = 50)
            curr.pred.sub[idx] <- X.i %*% curr.sub[, i]
            res[[i]][, k] <- curr.sub[, i]
        }

        if (k %% 1000 == 0) {
            cat(k, " samples generated.\n")
        }
    }
    return(res)
}


#' Fit a monotonic linear B-spline model.
#'
#' \code{BasisLin} fit a monotonic linear B-spline model which gives population
#' and subject-specific curves.
#'
#' This model is fitted in a Bayesian framework, that is, this routine gives the
#' samples drawing from the posterior distribution associated to each regression
#' coefficients.
#'
#' @param data First and second columns correspond to explanatory (x) and
#'     response (y) data respectively. The third column specifies the subject
#'     (group) to each datapoint corresponds to. Required.
#'
#' @param K The number of (inner) knots. The knots are equally spaced between
#'     \code{range(data[[1]])}. Required.
#'
#' @param shape Specify the shape of the response curves. Default is (monotonic)
#'     "increasing". Possible values are "increasing", "decreasing" ...
#'
#' @param size The numbre of samples to be drawn from the posterior.
#'
#' @param burn The number of samples to burn before recording.
#'
#' @return A list of (K + 1) matrices containing samples from the
#'     posteriors. The dimension of each matrix is \code{K + 2} \times
#'     \code{size}.
SubjectsBs <- function(data, K, shape = "increasing", size = 100,
                     burn = size / 10) {

    deg <- 1

    x <- data[[1]]
    y <- data[[2]]
    grp <- data[[3]]

    idx.sub <- tapply(seq_len(nrow(data)), grp, function(x) x)
    n.subject <- length(idx.sub)
    n.terms <- K + 2


    ## Construct design matrix and its cross-products
    dist <- diff(range(x)) / (K + 1)
    knots <- seq(min(x) - deg*dist, max(x) + deg*dist, len = K + 2*(deg + 1))
    B.pop <- splines::splineDesign(knots, x, ord = deg + 1)
    B.pop.sq <- crossprod(B.pop)
    B.sub.sq <- list()
    for (i in seq_len(n.subject)) {
        B.sub.sq[[i]] <- crossprod(B.pop[idx.sub[[i]], ])
    }
    names(B.sub.sq) <- names(idx.sub)

    ## Construct the constraint matrix and its cross-product
    D <- matrix(0, n.terms - deg, n.terms)
    if (shape == "increasing") {
        D[col(D) == row(D)] <- -1
        D[col(D) == (row(D) + 1)] <- 1
    } else if (shape == "decreasing") {
        D[col(D) == row(D)] <- 1
        D[col(D) == (row(D) + 1)] <- -1
    }
    D.sq <- crossprod(D)
    D.t <- t(D)

    ## Initial (current) values
    c.pop <- seq(0, 1, len = n.terms)
    c.sub <- matrix(0, n.terms, n.subject)
    c.veps <- 0.077^2
    c.vpop <- 4.29^2
    c.vsub <- 1.74^2

    ## Initial prediction contribution by subject specific effects
    pred.sub <- rep(0, nrow(data))
    for (i in seq_len(n.subject)) {
        idx <- idx.sub[[i]]
        pred.sub[idx] <- B.pop[idx, ] %*% c.sub[, i]
    }


    ## Initialise the output list, ordered by the order of levels(grp)
    ## "samples[[i]]" to access the i^th subject (in levels(grp)) curve
    ## "samples$pop" to access the population curve
    samples <- rep(list(matrix(NA, n.terms, size)), n.subject)
    names(samples) <- names(idx.sub)
    samples$population <- matrix(NA, n.terms, size)

    ## Burnin step
    for (k in seq_len(burn + size)) {
        ## Update population estimates
        M.pop <- solve(B.pop.sq + (c.veps / c.vpop) * D.sq)
        mu.pop <- tcrossprod(M.pop, B.pop) %*% (y - pred.sub)
        sig.pop <- c.veps * M.pop

        lower.pop <- D %*% c.sub
        lower.pop <- -1 * apply(lower.pop, 1, min)
        lower.pop[lower.pop < 0] <- 0

        c.pop <- new.rmvtgauss.lin(1, mu.pop, sig.pop,
                                   Amat = D.t,
                                   Avec = lower.pop,
                                   start = c.pop,
                                   burnin = 50)
        if (k > burn) {
            samples$population[, k - burn] <- c.pop
        }

        ## Update prediction contribution by the population curve.
        pred.pop <- B.pop %*% c.pop
        y.diff.pop <- y - pred.pop

        ## Update subject specific estimates
        lower.sub <- -1 * (D %*% c.pop)

        for (i in seq_len(n.subject)) {
            idx <- idx.sub[[i]]
            M.sub <- solve(B.sub.sq[[i]] + diag(c.veps / c.vsub, n.terms))
            mu.sub <- tcrossprod(M.sub, B.pop[idx, ]) %*% y.diff.pop[idx]
            sig.sub <- c.veps * M.sub

            c.sub[, i] <- new.rmvtgauss.lin(1, mu.sub, sig.sub,
                                            Amat = D.t,
                                            Avec = lower.sub,
                                            start = c.sub[, i],
                                            burnin = 50)

            ## Update prediction contribution by subject curves
            pred.sub[idx] <- B.pop[idx, ] %*% c.sub[, i]

            if (k > burn) {
                samples[[i]][, k - burn] <- c.sub[, i]
            }
        }
    }

    ## Return posterior mean, samples, and information regarding the basis
    ## functions.
    means <- lapply(samples, rowMeans)
    basis <- list(type = "bs", knots = knots, degree = deg)
    res <- list(means = means, samples = samples, basis = basis)
    return(res)
}

