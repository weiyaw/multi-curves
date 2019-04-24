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
sub_bs <- function(data, K, deg = 1, shape = "increasing", size = 100,
                   burn = size / 10, sampler = "gibbs", verbose) {
    EPS <- 1e-6

    if (deg > K + 1) {
        stop("Number of knots not sufficient for spline order.")
    }

    if (burn < 0) {
        stop("Negative burn.")
    } else if (burn == 0) {
        warning("Not burning any samples.")
    }

    x <- data[[1]]
    y <- data[[2]]
    ## convert the group variable into a factor
    if (is.factor(data[[3]])) {
        grp <- droplevels(data[[3]])
    } else {
        grp <- factor(data[[3]], levels = unique(data[[3]]))
    }

    n_terms <- K + deg + 1              # number of basis functions
    n_samples <- NROW(data)             # number of samples

    lvl_sub <- levels(grp)
    idx_sub <- tapply(seq_len(n_samples), grp, function(x) x)
    n_subs <- length(idx_sub)


    ## construct design matrix and its cross-products
    design_ls <- get_design_bs(x, K, deg, EPS)
    knots <- design_ls$knots
    B_pop <- design_ls$design

    ## get rid of the design_ls to save space
    rm(design_ls)

    ## construct cross-products of the design matrix
    B_pop_sq <- crossprod(B_pop)
    B_sub_sq <- array(NA, c(n_terms, n_terms, n_subs),
                      list(NULL, NULL, lvl_sub))
    for (j in lvl_sub) {
        B_sub_sq[, , j] <- crossprod(B_pop[idx_sub[[j]], ])
    }

    ## construct the difference matrix and its cross-product
    D <- get_diff_mat(n_terms, deg)

    ## construct the constraint matrix
    A <- get_constmat_bs(n_terms, shape)
    A_t <- t(A)

    ## calculate an inverse of A to easily produce feasible states
    ## an extra row is added on top of the A matrix to make it invertible.
    ## A.inv <- solve(rbind(c(1, rep(0, n.terms - 1)), A))
    A_inv <- diag(NCOL(A))
    A_inv[row(A_inv) > diff(dim(A))] <- A
    A_inv <- solve(A_inv)

    ## initialise population coefs and subjects deviations
    kcoef_pop <- get_ols(y, B_pop)
    kcoef_sub <- matrix(0, n_terms, n_subs, dimnames = list(NULL, lvl_sub))

    ## initialise prediction contribution by population coefs
    kpred_pop <- B_pop %*% kcoef_pop

    ## initialise prediction contribution by subjects deviations
    kpred_sub <- rep(NA, n_samples)
    for (j in lvl_sub) {
        idx <- idx_sub[[j]]
        kpred_sub[idx] <- B_pop[idx, ] %*% kcoef_sub[, j]
    }

    ## initialise the output list, by the order of lvl_sub
    samples <- list(population = matrix(NA, n_terms, size),
                    subjects = array(NA, c(n_terms, n_subs, size),
                                     dimnames = list(NULL, lvl_sub)),
                    precision = list())
    ## samples$precision: precisions (matrix) (num vec/mat list)
    ## pop: population coefs; poly: subject polynomial terms;
    ## sub: subject deviations; eps: error terms.
    samples$precision$pop <- rep(NA, size)
    samples$precision$sub <- rep(NA, size)
    samples$precision$eps <- rep(NA, size)
    loglike_ls <- rep(NA, size)

    ## Burnin step
    for (k in seq.int(-burn + 1, size)) {
        ## Update variances
        kprecs <- get_cov_bs(list(kcoef_pop), list(kcoef_sub), list(kpred_pop),
                             list(kpred_sub), list(y), D, 1,
                             n_subs, n_terms, n_samples)

        if (k > 0) {
            ## store the precisions
            samples$precision$pop[k] <- kprecs$pop
            samples$precision$poly[, , k] <- kprecs$poly
            samples$precision$sub[k] <- kprecs$sub
            samples$precision$eps[k] <- kprecs$eps
        }

        kcoefs <- get_coefs_gibbs(kcoef_pop, kcoef_sub, kpred_sub,
                                  B_pop, B_pop_sq, B_sub_sq,
                                  lvl_sub, idx_sub, y,
                                  kprecs$eps, kprecs$pop, kprecs$sub,
                                  A, A_t, A_inv, D, n_terms)

        ## for the ease of reading
        kcoef_pop <- kcoefs$coef_pop
        kcoef_sub <- kcoefs$coef_sub
        kpred_pop <- kcoefs$pred_pop
        kpred_sub <- kcoefs$pred_sub
        if (k > 0) {
            ## store the results
            samples$precision$pop[k] <- kprecs$pop
            samples$precision$poly[, , k] <- kprecs$poly
            samples$precision$sub[k] <- kprecs$sub
            samples$precision$eps[k] <- kprecs$eps
            samples$population[, k] <- kcoef_pop
            samples$subjects[, , k] <- kcoef_sub

            ## store the log-likelihood
            loglike_ls[k] <- mvtnorm::dmvnorm(y, kpred_pop + kpred_sub,
                                              diag(1/kprecs$eps, n_samples),
                                              log = TRUE)
        }

        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }
    }

    ## return posterior means (coefs only), samples, and information regarding
    ## the basis functions.
    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    basis <- list(type = "bs", knots = knots, degree = deg)
    info <- list(lvl_pop = NULL, lvl_sub = lvl_sub, n_terms = n_terms)
    data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
    mle_idx <- which.max(loglike_ls)
    mle <- list(population = samples$population[, mle_idx],
                subjects = samples$subjects[, , mle_idx],
                loglike = loglike_ls[mle_idx])

    list(means = means, samples = samples, basis = basis, info = info,
         data = data, mle = mle)
}


## get a sample from the coefs posterior (one population)

## different in each ITERATION and POPULATION
## coef_pop: previous population estimate (num vec)
## coef_sub: previous individual deviations (num mat)
## pred_sub: prediction contribution by the sub deviations (num vec)

## different in each POPULATION
## X_pop: design matrix (num mat)
## X_pop_sq: crossprod of the whole design matrix (num mat)
## X_sub_sq: crossprod of the individual design matrix (num 3d ary)
## lvl_sub: names of each subject (str vec)
## idx_sub: indices of each subject in X_pop (num vec list)
## y_pop: response data (num vec list)

## different in each ITERATAION
## prc_eps: precision of the Gaussian noise (num)
## prc_pop: precision of the population spline terms (num)
## prc_sub: precision of the individual spline terms (num)

## constants
## A: constraint matrix (num mat)
## A_t: t(A), for the purpose of using Berwin's sampler (num mat)
## A_inv: pseudo-inverse of A for generating feasible starting values (num mat)
## Kmat: the penalty matrix (num mat)
## n_terms: number of parameters to descrive a curve (num)

## RETURN
## coef_pop: population coefs (num vec)
## coef_sub: individual deviations (num mat)
## pred_pop: prediction contribution by the pop coefs (num vec)
## pred_sub: prediction contribution by the sub deviations (num vec)
get_coefs_gibbs <- function(coef_pop, coef_sub, pred_sub, X_pop, X_pop_sq, X_sub_sq,
                            lvl_sub, idx_sub, y_pop,
                            prc_eps, prc_pop, prc_sub, A,
                            A_t, A_inv, Kmat, n_terms, DEBUG = FALSE) {

    M_pop <- solve(X_pop_sq + (prc_pop / prc_eps) * crossprod(Kmat))
    mu_pop <- tcrossprod(M_pop, X_pop) %*% (y_pop - pred_sub)
    sig_pop <- M_pop / prc_eps

    lower_pop <- A %*% coef_sub
    lower_pop <- -1 * apply(lower_pop, 1, min)
    lower_pop[lower_pop < 0] <- 0

    ## Initialise the starting values of the truncated normal sampler
    start_pop <- A_inv %*% c(mu_pop[1], (lower_pop + 0.1))

    coef_pop <- start_pop
    if (DEBUG) {browser()}

    ## unconstrained sampler
    ## coef_pop <- t(mvtnorm::rmvnorm(1, mu_pop, sig_pop))

    ## Berwin's constrained sampler
    ## coef_pop <- TruncatedNormal::rmvtgauss.lin(1, mu_pop, sig_pop,
    ##                                         Amat = A_t,
    ##                                         Avec = lower_pop,
    ##                                         ## start = start_pop,
    ##                                         ## burnin = 1000)
    ##                                         start = coef_pop,
    ##                                         burnin = 500)

    ## modified HMC constrained sampler
    coef_pop <- t(tnorm::rmvtnorm(1, mu_pop, sig_pop,
                                F = A,
                                g = -lower_pop,
                                ## start = start_pop,
                                ## burnin = 1000)
                                initial = coef_pop,
                                burn = 20))


    ## Update prediction contribution by the population curve.
    pred_pop <- X_pop %*% coef_pop
    y_diff_pop <- y_pop - pred_pop

    ## Update subject specific estimates
    lower_sub <- -1 * (A %*% coef_pop)

    ## initialise the starting values for the truncated normal sampler, and
    ## the coef_sub matrix
    zeros <- rep(0, n_terms)
    coef_sub <- matrix(0, n_terms, length(lvl_sub),
                       dimnames = list(NULL, lvl_sub))

    ## Calculate the precision matrix term
    half_N <- diag(prc_sub, n_terms)
    half_N <- half_N / prc_eps

    for (j in lvl_sub) {
        idx <- idx_sub[[j]]
        M_sub <- solve(X_sub_sq[, , j] + half_N)
        mu_sub <- tcrossprod(M_sub, X_pop[idx, ]) %*% y_diff_pop[idx]
        sig_sub <- M_sub / prc_eps

        ## unconstrained sampler
        ## coef_sub[, j] <- mvtnorm::rmvnorm(1, mu_sub, sig_sub)

        ## Berwin's constrained sampler
        ## tmp <- TruncatedNormal::rmvtgauss.lin(500, mu_sub, sig_sub,
        ##                                       Amat = A_t,
        ##                                       Avec = lower_sub,
        ##                                       ## start = zeros,
        ##                                       ## burnin = 1000)
        ##                                       start = coef_sub[, j],
        ##                                       burnin = 0)

        ## modified HMC constrained sampler
        coef_sub[, j] <- tnorm::rmvtnorm(1, mu_sub, sig_sub,
                                         F = A,
                                         g = -lower_sub,
                                         ## start = zeros,
                                         ## burnin = 1000)
                                         initial = coef_sub[, j],
                                         burn = 20)


        ## Update prediction contribution by subject curves
        pred_sub[idx] <- X_pop[idx, ] %*% coef_sub[, j]
    }

    list(coef_pop = coef_pop, coef_sub = coef_sub, pred_pop = pred_pop,
         pred_sub = pred_sub, loglike = NA)
}



## get a sample from the precision posterior (one/multiple population)

## different in each ITERATION and POPULATION
## coef_pop: population coefs (num vec list)
## coef_sub: individual deviations (num mat list)
## pred_pop: prediction contribution by the pop coefs (num vec list)
## pred_sub: prediction contribution by the sub deviations (num vec list)

## different in each POPULATION
## y_pop: response data (num vec list)

## constants
## Kmat: the penalty matrix (num mat)
## n_pops: number of populations (num)
## n_subs: number of subjects (num vec)
## n_spline: number of spline terms/knots (num)
## n_samples: total sample size (num)

## RETURN
## prc_eps: precision of the Gaussian noise (num)
## prc_pop: precision of the population terms (num)
## prc_sub: precision of the individual terms (num)
get_cov_bs <- function(coef_pop, coef_sub, pred_pop, pred_sub, y_pop,
                       Kmat, n_pops, n_subs, n_spline, n_samples) {

    ## hyperparameters of priors
    ig_a <- 0.5
    ig_b <- 0.0001

    n_pen <- NROW(Kmat)                 # number of penalised terms

    if (!all(is.list(coef_pop), is.list(coef_sub))) {
        stop("coefs must be lists.")
    }
    if (!all(is.list(pred_pop), is.list(pred_sub))) {
        stop("preds must be lists.")
    }
    if (!all(is.list(y_pop))) {
        stop("response must be a list.")
    }

    ## crossprod of the spline terms (assuming independent coefs)
    xspl_pop <- sum(vapply(coef_pop, function(x) crossprod(Kmat %*% x), 0))
    xspl_sub <- sum(vapply(coef_sub, function(x) crossprod(c(Kmat %*% x)), 0))

    ## precision of population terms
    shp_pop <- 0.5 * n_pops * n_pen + ig_a
    scl_pop <- 0.5 * xspl_pop + ig_b
    prc_pop <- rgamma(1, shape = shp_pop, rate = scl_pop)

    ## precision of individual terms (assuming independent coefs)
    shp_sub <- 0.5 * sum(n_subs) * n_spline + ig_a
    scl_sub <- 0.5 * xspl_sub + ig_b
    ## if (scl_sub < 0.00001) {scl_sub <- 1}
    prc_sub <- rgamma(1, shape = shp_sub, rate = scl_sub)

    ## precision of residuals
    resid_vec <- mapply(function(x, y, z) x - y - z,
                        y_pop, pred_pop, pred_sub,
                        SIMPLIFY = FALSE)
    shp_eps <- 0.5 * n_samples + ig_a
    scl_eps <- 0.5 * crossprod(unlist(resid_vec)) + ig_b
    prc_eps <- rgamma(1, shape = shp_eps, rate = scl_eps)

    list(pop = prc_pop, sub = prc_sub, eps = prc_eps)
}










### Parallel ### outdated ###

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

## problematic marginal covariance posterior
marginal_cov <- function(n, burn = 100,
                         data = list(Bmat = NULL, y = NULL, K = NULL),
                         start = list(pop = 0.1, eps = 0.1)) {

    if (any(is.null(data$Bmat), is.null(data$y), is.null(data$K))) {
        stop("Data incomplete.")
    }

    if (NCOL(data$K) == NCOL(data$Bmat)) {
        n_terms <- NCOL(data$K)
    } else {
        stop("Dim of K and Bmat mismatch.")
    }

    if (NROW(data$y) == NROW(data$Bmat)) {
        n_samples <- NROW(data$y)
    } else  {
        stop("Dim of y and Bmat mismatch.")
    }

    xBmat <- crossprod(data$Bmat)
    xK <- crossprod(data$K)
    xy <- crossprod(y)
    Bmatxy <- crossprod(data$Bmat, y)

    ## starting values
    kvar <- start

    ## log of the prior distribution
    log_prior <- function(kvar) {
        pop <- kvar$pop
        eps <- kvar$eps
        a <- -1/2                       # use a flat prior IG(-1/2, 0)
        b <- 0
        if (all(pop > 0, eps > 0)) {
            log_pop <- (-a-1) * log(pop) + (-2 * b / pop) # IG(a, b)
            log_eps <- (-a-1) * log(eps) + (-2 * b / eps)
            log_pop + log_eps
        } else {
            -Inf
        }
    }
    debug <- F
    ## log-likelihood
    log_like <- function(kvar) {
        pop <- kvar$pop
        eps <- kvar$eps
        L_inv <- solve(xK * pop + xBmat * eps)

        if (all(pop > 0, eps > 0)) {
            -0.5 * (n_terms * log(pop) + n_samples * log(eps)) -
                0.5 / eps * (crossprod(Bmatxy, L_inv %*% Bmatxy) - xy)
        } else {
            stop()
            -Inf
        }
    }

    ## independent proposal TN(0, 1)
    log_q <- function(kvar) {
        sig2 <- 0.1
        -0.5 / sig2 * crossprod(as.numeric(kvar))
    }

    samples <- list(pop = rep(NA, n), eps = rep(NA, n))

    ## MH algorithm, independent TN(0, 1)
    for (k in seq.int(-burn + 1, n)) {
        ## cat(k)
        ## if (k == 8) {debug <- T; browser()}
        kvar_tmp <- tmvtnorm::rtmvnorm(1, sigma = diag(0.1, 2), lower = c(0, 0))
        kvar_tmp <- list(pop = kvar_tmp[1], eps = kvar_tmp[2])

        prob <- exp(log_prior(kvar_tmp) + log_like(kvar_tmp) -
                    log_prior(kvar) - log_like(kvar) +
                    log_q(kvar) - log_q(kvar_tmp))

        if (prob > runif(1)) {
            kvar <- kvar_tmp
        }
        if (k > 0) {
            samples$pop[k] <- kvar$pop
            samples$eps[k] <- kvar$eps
        }
    }
    samples
}
