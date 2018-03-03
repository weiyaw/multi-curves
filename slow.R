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

