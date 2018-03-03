source("subor.R")
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
#' @param deg The degree of spline polynomial
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
SubjectsBs <- function(data, K, deg = 1, shape = "increasing", size = 100,
                     burn = size / 10) {
    EPS <- 1e-6

    if (deg > K + 1) {
        stop("Number of knots not sufficient for spline order.")
    }

    x <- data[[1]]
    y <- data[[2]]
    grp <- data[[3]]

    idx.sub <- tapply(seq_len(nrow(data)), grp, function(x) x)
    n.subject <- length(idx.sub)
    n.terms <- K + deg + 1

    ## Construct design matrix and its cross-products
    design <- BSplineDesign(x, K, deg, EPS)
    knots <- design$knots
    B.pop <- design$design
    B.pop.sq <- crossprod(B.pop)
    rm(design)

    B.sub.sq <- list()
    for (i in seq_len(n.subject)) {
        B.sub.sq[[i]] <- crossprod(B.pop[idx.sub[[i]], ])
    }
    names(B.sub.sq) <- names(idx.sub)

    ## Construct the difference matrix and its cross-product and transposition
    D <- DiffMat(n.terms, deg)
    D.sq <- crossprod(D)

    ## Construct the constraint matrix
    ## An extra row is added on top of the A matrix to make it invertible.
    A <- BSplineConstMat(n.terms, shape)
    A.t <- t(A)

    ## Calculate an inverse of A to easily produce feasible states
    ## A.inv <- solve(rbind(c(1, rep(0, n.terms - 1)), A))
    A.inv <- diag(NCOL(A))
    A.inv[row(A.inv) > diff(dim(A))] <- A
    A.inv <- solve(A.inv)

    ## Initial (current) values
    c.pop <- seq(0, 1, len = n.terms)
    c.sub <- matrix(0, n.terms, n.subject)

    ## Initial prediction contribution by population effects
    pred.pop <- B.pop %*% c.pop

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
    samples$var.pop <- rep(NA, size)
    samples$var.sub <- rep(NA, size)
    samples$var.eps <- rep(NA, size)

    ## Burnin step
    for (k in seq_len(burn + size)) {
        ## Update variances
        c.vpop <- 1 / rgamma(1, shape = NROW(D) / 2 + 0.001,
                             scale = 0.5 * crossprod(D %*% c.pop) + 0.001)
        c.vsub <- 1 / rgamma(1, shape = (n.subject * n.terms) / 2 + 0.001,
                             scale = 0.5 * crossprod(c(c.sub)) + 0.001)
        resid.vec <- y - pred.pop - pred.sub
        c.veps <- 1 / rgamma(1, shape = 0.5 * length(y) + 0.001,
                             scale = 0.5 * crossprod(resid.vec) + 0.001)

        if (k == 1) {
            c.vpop <- 1
            c.vsub <- 1
            c.veps <- 1
        }

        ## Update population estimates
        M.pop <- solve(B.pop.sq + (c.veps / c.vpop) * D.sq)
        mu.pop <- tcrossprod(M.pop, B.pop) %*% (y - pred.sub)
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
        ## c.pop <- tmvtnorm:::rtmvnorm2(1, mean = as.vector(mu.pop),
        ##                              sigma = sig.pop + diag(1, ncol(sig.pop)),
        ##                              D = A,
        ##                              lower = lower.pop,
        ##                              upper = rep(Inf, NROW(D)),
        ##                              start.value = start.pop,
        ##                              burn.in = 1000)

        if (k > burn) {
            rec.idx <- k - burn
            if (rec.idx %% 1000 == 0) {
                cat(rec.idx, "samples generated. \n")
            }
            samples$population[, rec.idx] <- c.pop
            samples$var.pop[rec.idx] <- c.vpop
            samples$var.sub[rec.idx] <- c.vsub
            samples$var.eps[rec.idx] <- c.veps
        }

        ## Update prediction contribution by the population curve.
        pred.pop <- B.pop %*% c.pop
        y.diff.pop <- y - pred.pop

        ## Update lower bound
        lower.sub <- -1 * (A %*% c.pop)

        ## Initialise the starting values of the truncated normal sampler
        ## start.sub <- rep(0, n.terms)
        zeros <- rep(0, n.terms)

        for (i in seq_len(n.subject)) {
            idx <- idx.sub[[i]]
            M.sub <- solve(B.sub.sq[[i]] + diag(c.veps / c.vsub, n.terms))
            mu.sub <- tcrossprod(M.sub, B.pop[idx, ]) %*% y.diff.pop[idx]
            sig.sub <- c.veps * M.sub

            ## ## for rtmg
            ## Minv.sub <- B.sub.sq[[i]] + diag(c.veps / c.vsub, n.terms)
            ## prec.sub <- Minv.sub / c.veps
            ## X.i <- B.pop[idx, , drop = FALSE]
            ## r.i <- crossprod(X.i %*% c.pop - y[idx], X.i) / -c.veps

            ## for (ii in 1:1000) {
            c.sub[, i] <- TruncatedNormal::rmvtgauss.lin(1, mu.sub, sig.sub,
                                                         Amat = A.t,
                                                         Avec = lower.sub,
                                                         start = zeros,
                                                         burnin = 1000)

            ## c.sub[, i] <- tmg::rtmg(1, M = prec.sub,
            ##                         r = as.vector(r.i),
            ##                         initial = rep(0, n.terms),
            ##                         f = A,
            ##                         g = as.vector(-1 * lower.sub),
            ##                         burn.in = 30)

            ## c.sub[, i] <- tmvtnorm::rtmvnorm2(1, mean = as.vector(mu.sub),
            ##                                   sigma = sig.sub,
            ##                                   D = A,
            ##                                   lower = lower.sub,
            ##                                   upper = rep(Inf, NROW(A)),
            ##                                   algorithm = "gibbs",
            ##                                   start.value = as.vector(start.sub),
            ##                                   burn.in = 1000)

            ## Update prediction contribution by subject curves
            pred.sub[idx] <- B.pop[idx, ] %*% c.sub[, i]

            if (k > burn) {
                samples[[i]][, k - burn] <- c.sub[, i]
            }
        }
    }

    ## Return posterior mean, samples, and information regarding the basis
    ## functions.
    means <- lapply(samples[seq_len(n.subject + 1)], rowMeans)
    basis <- list(type = "bs", knots = knots, degree = deg)
    info <- list(n.subject = n.subject, n.terms = n.terms)
    res <- list(means = means, samples = samples, basis = basis, info = info)
    return(res)
}


### Parallel ###

SubjectsBsPar <- function(data, K, deg = 1, shape = "increasing", size = 100,
                          burn = size / 10) {
    library(parallel)
    EPS <- 1e-6

    if (deg > K + 1) {
        stop("Number of knots not sufficient for spline order.")
    }

    x <- data[[1]]
    y <- data[[2]]
    grp <- data[[3]]

    idx.sub <- tapply(seq_len(nrow(data)), grp, function(x) x)
    n.subject <- length(idx.sub)
    n.terms <- K + deg + 1


    ## Construct design matrix and its cross-products
    dist <- diff(range(x)) / (K + 1)
    knots <- seq(min(x) - deg*dist - EPS, max(x) + deg*dist + EPS,
                 len = K + 2*(deg + 1))
    B.pop <- splines::splineDesign(knots, x, ord = deg + 1)
    B.pop.sq <- crossprod(B.pop)
    B.sub.sq <- list()
    for (i in seq_len(n.subject)) {
        B.sub.sq[[i]] <- crossprod(B.pop[idx.sub[[i]], ])
    }
    names(B.sub.sq) <- names(idx.sub)

    ## Construct the difference matrix and its cross-product and transposition
    D <- diag(1, n.terms)
    D[(col(D) + 1) == row(D)] <- -1
    D <- D[-1, ]
    for (i in seq_len(deg - 1)) {
        D <- D[-1, -1] %*% D
    }

    D.sq <- crossprod(D)

    ## Construct the constraint matrix
    ## An extra row is added on top of the A matrix to make it invertible.
    if (shape == "increasing") {
        A <- diag(1, n.terms)
        A[(col(A) + 1) == row(A)] <- -1
        A.inv <- solve(A)
        A <- A[-1, ]
    } else if (shape == "decreasing") {
        A <- diag(-1, n.terms)
        A[1, 1] <- 1
        A[(col(A) + 1) == row(A)] <- 1
        A.inv <- solve(A)
        A <- A[-1, ]
    }
    A.t <- t(A)

    ## Initial (current) values
    c.pop <- seq(0, 1, len = n.terms)
    c.sub <- matrix(0, n.terms, n.subject)
    ## c.veps <- 0.077^2
    ## c.vpop <- 4.29^2
    ## c.vsub <- 1.74^2


    ## Initial prediction contribution by population effects
    pred.pop <- B.pop %*% c.pop

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
    samples$var.pop <- rep(NA, size)
    samples$var.sub <- rep(NA, size)
    samples$var.eps <- rep(NA, size)

    ## Burnin step
    for (k in seq_len(burn + size)) {
        ## Update variances
        c.vpop <- 1 / rgamma(1, shape = NROW(D) / 2 + 0.001,
                             scale = 0.5 * crossprod(D %*% c.pop) + 0.001)
        c.vsub <- 1 / rgamma(1, shape = (n.subject * n.terms) / 2 + 0.001,
                             scale = 0.5 * crossprod(c(c.sub)) + 0.001)
        resid.vec <- y - pred.pop - pred.sub
        c.veps <- 1 / rgamma(1, shape = 0.5 * length(y) + 1,
                             scale = 0.5 * crossprod(resid.vec) + 0.001)

        ## Update population estimates
        M.pop <- solve(B.pop.sq + (c.veps / c.vpop) * D.sq)
        mu.pop <- tcrossprod(M.pop, B.pop) %*% (y - pred.sub)
        sig.pop <- c.veps * M.pop

        lower.pop <- A %*% c.sub
        lower.pop <- -1 * apply(lower.pop, 1, min)
        lower.pop[lower.pop < 0] <- 0

        ## Initialise the starting values of the truncated normal sampler
        start.pop <- A.inv %*% c(mu.pop[1], (lower.pop + 0.1))

        c.pop <- new.rmvtgauss.lin(1, mu.pop, sig.pop,
                                   Amat = A.t,
                                   Avec = lower.pop,
                                   start = start.pop,
                                   burnin = 100)
        ## c.pop <- rtmvnorm2(1, mean = as.vector(mu.pop),
        ##                    sigma = sig.pop,
        ##                    D = A,
        ##                    lower = lower.pop,
        ##                    upper = rep(Inf, NROW(D)),
        ##                    algorithm = "gibbs",
        ##                    start.value = start.pop,
        ##                    burn.in = 50)

        if (k > burn) {
            rec.idx <- k - burn
            samples$population[, rec.idx] <- c.pop
            samples$var.pop[rec.idx] <- c.vpop
            samples$var.sub[rec.idx] <- c.vsub
            samples$var.eps[rec.idx] <- c.veps
        }

        ## Update prediction contribution by the population curve.
        pred.pop <- B.pop %*% c.pop
        y.diff.pop <- y - pred.pop

        ## Update lower bound
        lower.sub <- -1 * (A %*% c.pop)

        ## Initialise the starting values of the truncated normal sampler
        start.sub <- A.inv %*% c(mu.pop[1], (lower.sub + 0.1))


        ## Helper function to facilitate parallel computing. "idx" and
        ## "B.sub.sq" are lists.
        genSubSamps <- function(idx, B.sub.sq, B.pop, y.diff.pop,
                                c.vsub, c.veps, n.terms, A.t,
                                lower.sub, start.sub) {

            M.sub <- solve(B.sub.sq + diag(c.veps / c.vsub, n.terms))
            mu.sub <- tcrossprod(M.sub, B.pop[idx, ]) %*% y.diff.pop[idx]
            sig.sub <- c.veps * M.sub

            new.rmvtgauss.lin(1, mu.sub, sig.sub,
                              Amat = A.t, Avec = lower.sub,
                              start = start.sub, burnin = 100)

        }
        ## Update subject-specific estimates
        c.sub <- mapply(genSubSamps, idx.sub, B.sub.sq,
                          MoreArgs = list(B.pop = B.pop,
                                          y.diff.pop = y.diff.pop,
                                          c.vsub = c.vsub,
                                          c.veps = c.veps,
                                          n.terms = n.terms,
                                          A.t = A.t,
                                          lower.sub = lower.sub,
                                          start.sub = start.sub))
                          ## mc.cores = detectCores())

        for (i in seq_len(n.subject)) {
            idx <- idx.sub[[i]]
        ##     M.sub <- solve(B.sub.sq[[i]] + diag(c.veps / c.vsub, n.terms))
        ##     mu.sub <- tcrossprod(M.sub, B.pop[idx, ]) %*% y.diff.pop[idx]
        ##     sig.sub <- c.veps * M.sub

        ##     c.sub[, i] <- new.rmvtgauss.lin(1, mu.sub, sig.sub,
        ##                                     Amat = A.t,
        ##                                     Avec = lower.sub,
        ##                                     start = start.sub,
        ##                                     burnin = 100)

        ## c.sub[, i] <- rtmvnorm2(1, mean = as.vector(mu.sub),
        ##                         sigma = sig.sub,
        ##                         D = A,
        ##                         lower = lower.sub,
        ##                         upper = rep(Inf, NROW(D)),
        ##                         algorithm = "gibbs",
        ##                         start.value = start.sub,
        ##                         burn.in = 50)

        ## Update prediction contribution by subject curves
            pred.sub[idx] <- B.pop[idx, ] %*% c.sub[, i]

            if (k > burn) {
                samples[[i]][, k - burn] <- c.sub[, i]
            }
        }
    }


    ## Return posterior mean, samples, and information regarding the basis
    ## functions.
    means <- lapply(samples[seq_len(n.subject + 1)], rowMeans)
    basis <- list(type = "bs", knots = knots, degree = deg)
    info <- list(n.subject = n.subject, n.terms = n.terms)
    data <- data.frame(x = x, y = y, grps = factor(grp))
    res <- list(means = means, samples = samples, basis = basis, info = info,
                data = data)
    return(res)
}
