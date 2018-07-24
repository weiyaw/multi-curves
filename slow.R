ExactGibbsTmg <- function(y, linear.mat, spline.mat, lme.obj, size = 100,
                          burn = size / 10) {

    ## Code that use HMC to simulate truncated normal distribution

    ## y : response used in lme.obj (better be sorted according to subjects
    ##     indices)
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
    constrt.mat[upper.tri(constrt.mat)] <- 0
    constrt.mat <- cbind(0, constrt.mat)

    ## Indices of datapoints corresponding to each subject
    fact.sub <- factor(lme.obj$groups[, -1],
                       levels = unique(as.character(lme.obj$groups[, -1])))
    idx.sub <- tapply(seq_len(n), fact.sub, function(x) x)

    ## Model matrix
    model.mat <- cbind(linear.mat, spline.mat)

    ## Precision of the conditional posterior for SUBJECT
    prec.sub <- rep(list(solve(varcov.sub)), n.subject)
    for (i in seq_len(n.subject)) {
        prec.sub[[i]] <- (crossprod(model.mat[idx.sub[[i]], ]) / var.error) +
            prec.sub[[i]]
    }

    ## Precision of the conditional posterior for POPULATION
    prec.pop <- diag(c(rep(0, n.fixed), rep(1 / var.pop, n.splines)))
    prec.pop <- (crossprod(model.mat) / var.error) + prec.pop

    ## Initialise current estimates, EBLUPS as initial population response curve
    ## If the EBLUPS is not monotone, make it monotone.
    curr.pop <- as.vector(pop.coef)
    need <- -cumsum(curr.pop[-1])
    need[need < 0.01] <- 0
    curr.pop <- curr.pop + c(0, need * 1.1)

    curr.sub <- matrix(NA, n.terms, n.subject)

    ## Record current individual curves contribution to the prediction.
    ## Essentially: model.mat %*% coef.sub
    curr.pred.sub <- rep(NA, n)

    if (burn < 1) {
        stop("Must burn at least 1 sample.")
    }

    ## Sequence of the subject (indices for the "for" loop)
    seq.subject <- seq_len(n.subject)

    ## Initial value when simulating SUBJECT posteriors
    init.sub <- rep(0, n.terms)

    ## Burning period. Results are discarded
    for (k in seq_len(burn)) {

        ## Simulate subject posterior (conditioned on previous population)
        ## Lower-bound vector for the subject posterior
        avec.sub <- cumsum(curr.pop[-1])

        ## Generate an individual curve for each subject
        for (i in seq.subject) {
            idx <- idx.sub[[i]]
            X.i <- model.mat[idx, , drop = FALSE]
            r.i <- crossprod(X.i %*% curr.pop - y[idx], X.i) / -var.error
            curr.sub[, i] <- rtmg(1, M = prec.sub[[i]],
                                  r = as.vector(r.i),
                                  initial = init.sub,
                                  f = constrt.mat,
                                  g = avec.sub)
            if (is.null(curr.sub[1, i])) {
                browser()
            }
            curr.pred.sub[idx] <- X.i %*% curr.sub[, i]
        }

        if (k >  burn) {
            break
        }

        ## Simulated population posterior (conditioned on current subject)
        ## Lower-bound vector for the population posterior
        avec.pop <- apply(curr.sub[-1, , drop = FALSE], 2, cumsum)
        avec.pop <- apply(avec.pop, 1, min)
        avec.pop[avec.pop > 0] <- 0

        ## Generate a population curve
        r <- crossprod(curr.pred.sub - y, model.mat) / -var.error
        init.pop <- c(curr.sub[1, 1], 1 - avec.pop[1], -diff(avec.pop))
        curr.pop <- as.vector(rtmg(1, M = prec.pop,
                                   r = as.vector(r),
                                   initial = init.pop,
                                   f = constrt.mat,
                                   g = avec.pop))

        if (is.null(curr.pop[1])) {
            browser()
        }


        if (k %% 1000 == 0) {
            cat(i, " samples burned.\n")
        }

    }

    ## Initialise the output list
    ## "res[[i]]" to access i th subject curve
    ## "res$pop" to access population curve
    res <- rep(list(matrix(NA, n.terms, size)), n.subject)
    names(res) <- names(idx.sub)
    res$population <- matrix(NA, n.terms, size)

    browser()
    ## Generate posterior. Results are recorded
    for (k in seq_len(size)) {

        ## Simulated population posterior (conditioned on current subject)
        ## Lower-bound vector for the population posterior
        avec.pop <- apply(curr.sub[-1, , drop = FALSE], 2, cumsum)
        avec.pop <- apply(avec.pop, 1, min)
        avec.pop[avec.pop > 0] <- 0

        ## Generate a population curve
        r <- crossprod(curr.pred.sub - y, model.mat) / -var.error
        init.pop <- c(curr.sub[1, 1], 1 - avec.pop[1], -diff(avec.pop))
        curr.pop <- as.vector(rtmg(1, M = prec.pop,
                                   r = as.vector(r),
                                   initial = init.pop,
                                   f = constrt.mat,
                                   g = avec.pop)
                              res$population[, k] <- curr.pop)

        ## Simulate subject posterior (conditioned on previous population)
        ## Lower-bound vector for the subject posterior
        avec.sub <- cumsum(curr.pop[-1])

        ## Generate an individual curve for each subject
        for (i in seq.subject) {
            idx <- idx.sub[[i]]
            X.i <- model.mat[idx, , drop = FALSE]
            r.i <- crossprod(X.i %*% curr.pop - y[idx], X.i) / -var.error
            curr.sub[, i] <- rtmg(1, M = prec.sub[[i]], r = as.vector(r.i),
                                  initial = init.sub,
                                  f = constrt.mat,
                                  g = avec.sub)
            curr.pred.sub[idx] <- X.i %*% curr.sub[, i]
            res[[i]][, k] <- curr.sub[, i]
        }

        if (k %% 1000 == 0) {
            cat(k, " samples generated.\n")
        }
    }
    return(res)
}


ApproxGibbs <- function(y, linear.mat, spline.mat, lme.obj, size = 1000) {
    ## y : response
    ## linear.mat: model matrix for straight line
    ## spline.mat: model matrix for splines
    ## lme.obj: unconstrained lme object from nlme
    ## size: number of samples from the posterior distribution

    ## Initialise
    varcov.list <- lapply(lme.obj$modelStruct$reStruct, as.matrix)
    varcov.list <- lapply(varcov.list, function (x) { x * lme.obj$sigma^2 })
    ## WARNING THE FACTOR OF TAPPLY NEED TO BE FIXED
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

## post <- ApproxGibbs(y, X, Z, fit2, 5)
## hist(post[3,])

SubjectLin <- function(y, linear.mat, spline.mat, lme.obj, size = 100,
                       burn = size / 10) {
    ## y : response used in lme.obj (better be sorted according to subjects
    ##     indices)
    ## linear.mat: model matrix for straight line
    ## spline.mat: model matrix for splines
    ## lme.obj: unconstrained lme object from nlme
    ## size: number of samples from the posterior distribution

    if (burn < 1) stop("Must burn at least 1 sample.")

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

