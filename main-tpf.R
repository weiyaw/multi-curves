source("subor.R")
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


## Truncated power function
## data : 1st col x, 2nd col y, 3rd col groups
## K : number of quantile (inner) knots, or a vector of inner knots
## deg : degree of spline polynomial
SubjectsTpf <- function(data, K, deg = 1, shape = "increasing", size = 100,
                        burn = size / 10) {

    ## if (deg != 1 && deg != 2) {
    ##     stop("Invalid spline degree. Must be 1 or 2.")
    ## }

    x <- data[[1]]
    y <- data[[2]]
    grp <- data[[3]]

    n.terms <- K + deg + 1
    n.sample <- NROW(data)
    idx.sub <- tapply(seq_len(n.sample), grp, function(x) x)
    n.subject <- length(idx.sub)

    ## Construct design matrix and its cross-products
    design <- TpfDesign(x, K, deg)
    knots <- design$knots               # all knots (with extrema)
    n.spline <- length(knots) - 2       # number of inner knots (w/o extrema)
    X.pop <- design$design
    X.pop.sq <- crossprod(X.pop)
    rm(design)

    X.sub.sq <- list()
    for (i in seq_len(n.subject)) {
        X.sub.sq[[i]] <- crossprod(X.pop[idx.sub[[i]], ])
    }
    names(X.sub.sq) <- names(idx.sub)

    ## An idempotent matrix to extract the spline terms
    Kmat <- diag(c(rep(0, deg + 1), rep(1, n.spline)))
    idx.poly <- seq_len(deg + 1)        # index of polynomial terms

    ## Construct the constraint matrix and its cross-product
    A <- TpfConstMat(knots, shape, deg)
    A.t <- t(A)

    ## Calculate an inverse of A to easily produce feasible states
    A.inv <- diag(NCOL(A))
    A.inv[row(A.inv) > diff(dim(A))] <- A
    A.inv <- solve(A.inv)

    ## Initial (current) values
    c.pop <- seq(0, 1, len = n.terms)
    c.sub <- matrix(0, n.terms, n.subject)

    ## Initial prediction contribution by population effects
    pred.pop <- X.pop %*% c.pop

    ## Initial prediction contribution by subject specific effects
    pred.sub <- rep(0, nrow(data))
    for (i in seq_len(n.subject)) {
        idx <- idx.sub[[i]]
        pred.sub[idx] <- X.pop[idx, ] %*% c.sub[, i]
    }


    ## Initialise the output list, ordered by the order of levels(grp)
    ## "samples[[i]]" to access the i^th subject (in levels(grp)) curve
    ## "samples$pop" to access the population curve
    samples <- rep(list(matrix(NA, n.terms, size)), n.subject)
    names(samples) <- names(idx.sub)
    samples$population <- matrix(NA, n.terms, size)
    samples$var.pop <- rep(NA, size)
    samples$prec.poly <- array(NA, c(deg + 1, deg + 1, size))
    samples$var.sub <- rep(NA, size)
    samples$var.eps <- rep(NA, size)

    ## hyperparameters of priors
    ig.a <- 0.001
    ig.b <- 0.001
    wi.df <- 1
    wi.sig <- diag(0.001, deg + 1)

    ## Burnin step
    for (k in seq_len(burn + size)) {
        ## Update variances
        c.vpop <- 1 / rgamma(1, shape = n.spline / 2 + ig.a,
                             scale = 0.5 * crossprod(Kmat %*% c.pop) + ig.b)
        c.vsub <- 1 / rgamma(1, shape = (n.subject * n.spline) / 2 + ig.a,
                             scale = 0.5 * crossprod(c(Kmat %*% c.sub)) + ig.b)
        c.ppoly <- rWishart(1, df = wi.df + n.subject,
                            solve(wi.sig + tcrossprod(c.sub[idx.poly, ])))
        resid.vec <- y - pred.pop - pred.sub
        c.veps <- 1 / rgamma(1, shape = 0.5 * n.sample + ig.a,
                             scale = 0.5 * crossprod(resid.vec) + ig.b)


        ## Update population estimates
        M.pop <- solve(X.pop.sq + (c.veps / c.vpop) * Kmat)
        mu.pop <- tcrossprod(M.pop, X.pop) %*% (y - pred.sub)
        sig.pop <- c.veps * M.pop

        lower.pop <- A %*% c.sub
        lower.pop <- -1 * apply(lower.pop, 1, min)
        lower.pop[lower.pop < 0] <- 0

        ## Initialise the starting values of the truncated normal sampler
        start.pop <- A.inv %*% c(mu.pop[1], (lower.pop + 0.1))

        c.pop <- TruncatedNormal::rmvtgauss.lin(1, mu.pop, sig.pop,
                                                Amat = A.t,
                                                Avec = lower.pop,
                                                start = start.pop,
                                                burnin = 1000)
        if (k > burn) {
            rec.idx <- k - burn
            if (rec.idx %% 1000 == 0) {
                cat(rec.idx, "samples generated. \n")
            }
            samples$population[, rec.idx] <- c.pop
            samples$var.pop[rec.idx] <- c.vpop
            samples$prec.poly[, , rec.idx] <- c.ppoly
            samples$var.sub[rec.idx] <- c.vsub
            samples$var.eps[rec.idx] <- c.veps
        }

        ## Update prediction contribution by the population curve.
        pred.pop <- X.pop %*% c.pop
        y.diff.pop <- y - pred.pop

        ## Update subject specific estimates
        lower.sub <- -1 * (A %*% c.pop)

        ## Initialise the starting values of the truncated normal sampler
        ## start.sub <- rep(0, n.terms)
        zeros <- rep(0, n.terms)

        ## Calculate the precision matrix term
        G.term <- diag(1 / c.vsub, n.terms)
        G.term[1:NROW(c.ppoly), 1:NCOL(c.ppoly)] <- c.ppoly
        G.term <- c.veps * G.term

        for (i in seq_len(n.subject)) {
            idx <- idx.sub[[i]]
            M.sub <- solve(X.sub.sq[[i]] + G.term)
            mu.sub <- tcrossprod(M.sub, X.pop[idx, ]) %*% y.diff.pop[idx]
            sig.sub <- c.veps * M.sub

            c.sub[, i] <- TruncatedNormal::rmvtgauss.lin(1, mu.sub, sig.sub,
                                                         Amat = A.t,
                                                         Avec = lower.sub,
                                                         start = zeros,
                                                         burnin = 1000)

            ## Update prediction contribution by subject curves
            pred.sub[idx] <- X.pop[idx, ] %*% c.sub[, i]

            if (k > burn) {
                samples[[i]][, k - burn] <- c.sub[, i]
            }
        }
    }

    ## Return posterior mean, samples, and information regarding the basis
    ## functions.
    means <- lapply(samples[seq_len(n.subject + 1)], rowMeans)
    basis <- list(type = "tpf", knots = knots, degree = deg)
    info <- list(n.subject = n.subject, n.terms = n.terms)
    data <- data.frame(x = x, y = y, grps = factor(grp))
    res <- list(means = means, samples = samples, basis = basis, info = info,
                data = data)
    return(res)
}

