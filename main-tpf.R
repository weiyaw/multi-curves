source("subor.R")

## Truncated power function
## data : 1st col x, 2nd col y, 3rd col groups
## K : number of quantile (inner) knots, or a vector of inner knots
## deg : degree of spline polynomial
sub_tpf <- function(data, K, deg = 1, penalty = TRUE, shape = "increasing", size = 100,
                    burn = size / 10, verbose = FALSE) {

    if (deg != 1 && deg != 2) {
        stop("Invalid spline degree. Must be 1 or 2.")
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

    n_terms <- K + deg + 1
    n_samples <- NROW(data)

    lvl_sub <- levels(grp)
    idx_sub <- tapply(seq_len(n_samples), grp, function(x) x)
    n_subs <- length(idx_sub)

    ## construct the design matrix and knots
    design_ls <- get_design_tpf(x, K, deg)
    knots <- design_ls$knots            # all knots (with boundaries)
    n_spline <- length(knots) - 2       # number of inner knots (w/o boundaries)
    X_pop <- design_ls$design
    browser()
    ## get rid of the design_ls to save space
    rm(design_ls)

    ## construct cross-products of the design matrix
    X_pop_sq <- crossprod(X_pop)
    X_sub_sq <- array(NA, c(n_terms, n_terms, n_subs),
                      list(NULL, NULL, lvl_sub))
    for (j in lvl_sub) {
        X_sub_sq[, , j] <- crossprod(X_pop[idx_sub[[j]], ])
    }

    ## An idempotent matrix to extract the spline terms
    Kmat <- diag(c(rep(0, deg + 1), rep(1, n_spline)))
    idx_poly <- seq_len(deg + 1)        # index of polynomial terms

    ## Construct the constraint matrix and its cross-product
    A <- get_constmat_tpf(knots, shape, deg)
    A_t <- t(A)

    ## Calculate an inverse of A to easily produce feasible states
    A_inv <- diag(NCOL(A))
    A_inv[row(A_inv) > diff(dim(A))] <- A
    A_inv <- solve(A_inv)

    ## initialise population coefs and subjects deviations
    kcoef_pop <- get_ols(y, X_pop)
    kcoef_sub <- matrix(0, n_terms, n_subs, dimnames = list(NULL, lvl_sub))

    ## initialise prediction contribution by population coefs
    kpred_pop <- X_pop %*% kcoef_pop

    ## initialise prediction contribution by subjects deviations
    kpred_sub <- rep(NA, n_samples)
    for (j in lvl_sub) {
        idx <- idx_sub[[j]]
        kpred_sub[idx] <- X_pop[idx, ] %*% kcoef_sub[, j]
    }

    ## remove the dummy variables used in the for loop
    rm(j, idx)

    ## initialise the output list, by the order of lvl_sub
    samples <- list(population = matrix(NA, n_terms, size),
                    subjects = array(NA, c(n_terms, n_subs, size),
                                     dimnames = list(NULL, lvl_sub)),
                    precision = list())
    ## samples$precision: precisions (matrix) (num vec/mat list)
    ## pop: population coefs; poly: subject polynomial terms;
    ## sub: subject deviations; eps: error terms.
    samples$precision$pop <- rep(NA, size)
    samples$precision$poly <- array(NA, c(deg + 1, deg + 1, size))
    samples$precision$sub <- rep(NA, size)
    samples$precision$eps <- rep(NA, size)
    loglike_ls <- rep(NA, size)

    ## burnin followed by actual sampling
    for (k in seq.int(-burn + 1, size)) {
        ## get the precisions (varariance-covariance matrices)
        kprecs <- get_cov_tpf(list(kcoef_pop), list(kcoef_sub), list(kpred_pop),
                              list(kpred_sub), list(y), idx_poly, Kmat, 1,
                              n_subs, n_spline, n_samples)

        ## No penalty on the population spline terms
        if (!penalty) {kprecs$pop <- 0}

        ## get the coefs and deviations
        kcoefs <- get_coefs_hmc(kcoef_pop, kcoef_sub, kpred_sub,
                                X_pop, X_pop_sq, X_sub_sq,
                                lvl_sub, idx_sub, y,
                                kprecs$eps, kprecs$pop, kprecs$poly, kprecs$sub,
                                A, A_t, A_inv, Kmat, n_terms)

        ## for the ease of reading
        kcoef_pop <- kcoefs$coef_pop
        kcoef_sub <- kcoefs$coef_sub
        kpred_pop <- kcoefs$pred_pop
        kpred_sub <- kcoefs$pred_sub

        if (k > 0) {
            ## store the precisions
            samples$precision$pop[k] <- kprecs$pop
            samples$precision$poly[, , k] <- kprecs$poly
            samples$precision$sub[k] <- kprecs$sub
            samples$precision$eps[k] <- kprecs$eps

            ## store the coefs and deviations
            samples$population[, k] <- kcoef_pop
            samples$subjects[, , k] <- kcoef_sub

            ## store the log-likelihood
            loglike_ls[k] <- kcoefs$loglike
        }

        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }
    }

    ## return posterior means (coefs only), samples, and information regarding
    ## the basis functions.
    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    basis <- list(type = "tpf", knots = knots, degree = deg)
    info <- list(lvl_pop = NULL, lvl_sub = lvl_sub, n_terms = n_terms)
    data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
    mle_idx <- which.max(loglike_ls)
    mle <- list(population = samples$population[, mle_idx],
                subjects = samples$subjects[, , mle_idx],
                loglike = loglike_ls[mle_idx])

    list(means = means, samples = samples, basis = basis, info = info,
         data = data, mle = mle)
}


source("subor.R")


## Subjects model with multiple population curves
## data: 1st predictor, 2nd response, 3rd subjects, 4th populations
## K: number of quantile (inner) knots, or a vector of inner knots
## deg: degree of spline polynomial
multisub_tpf <- function(data, K, deg = 1, shape = "increasing", size = 100,
                         burn = size / 10, verbose = FALSE) {
    ## if (deg != 1 && deg != 2) {
    ##     stop("Invalid spline degree. Must be 1 or 2.")
    ## }
    DEBUG <- FALSE

    n_terms <- K + deg + 1              # number of coefs to describe a curve
    n_samples <- NROW(data)             # number of provided datapoints

    x <- data[[1]]
    y <- data[[2]]

    ## convert groups variables to factors
    if (is.factor(data[[3]])) {
        grp_sub <- droplevels(data[[3]])
    } else {
        grp_sub <- factor(data[[3]], levels = unique(data[[3]]))
    }
    if (is.factor(data[[4]])) {
        grp_pop <- droplevels(data[[4]])
    } else {
        grp_pop <- factor(data[[4]], levels = unique(data[[4]]))
    }

    ## calculate the population indices
    ## these are used to extract from the input (main) dataframe
    lvl_pop <- levels(grp_pop)
    idx_pop <- tapply(seq_len(n_samples), grp_pop, function(x) x)
    n_pops <- length(idx_pop)

    ## calculate the subject indices within each population
    ## these are used to extract from each population dataframe (P_i)
    lvl_sub <- list()
    idx_sub <- list()
    for (i in lvl_pop) {
        grp_sub_i <- grp_sub[grp_pop == i, drop = TRUE]
        lvl_sub[[i]] <- levels(grp_sub_i)
        idx_sub[[i]] <- tapply(seq_along(grp_sub_i), grp_sub_i , function(x) x)
    }
    n_subs <- vapply(idx_sub, length, 0L) # number of subjects in each pop

    ## remove the dummy variables used in the for loop
    rm(i, grp_sub_i)

    ## construct design matrix and knots
    design_ls <- get_design_tpf(x, K, deg)
    knots <- design_ls$knots            # all knots (with extrema)
    n_spline <- length(knots) - 2       # number of inner knots (w/o extrema)

    ## construct cross-products of the design matrix
    ## also, seperate the design matrix for each population (P_i)
    ## X_pop_sq = crossprod(P_i), X_sub_sq = crossprod(X_ij)
    ## X_pop_sq is a 3D array, X_sub_sq is a list of 3D arrays.
    X_pop <- list()
    X_pop_sq <- list()
    X_sub_sq <- list()
    for (i in lvl_pop) {
        X_pop[[i]] <- design_ls$design[idx_pop[[i]], ]
        X_pop_sq[[i]] <- crossprod(X_pop[[i]])
        X_sub_sq[[i]] <- array(NA, c(n_terms, n_terms, n_subs[[i]]),
                               list(NULL, NULL, lvl_sub[[i]]))
        for (j in lvl_sub[[i]]) {
            X_sub_sq[[i]][, , j] <- crossprod(X_pop[[i]][idx_sub[[i]][[j]], ])
        }
    }

    ## get rid of the design_ls to save space
    rm(design_ls)

    ## get y for each population (y_i)
    y_pop <- tapply(y, grp_pop, function(x) x)

    ## an idempotent matrix to extract the spline terms
    Kmat <- diag(c(rep(0, deg + 1), rep(1, n_spline)))
    idx_poly <- seq_len(deg + 1)        # index of polynomial terms

    ## construct the constraint matrix and its cross-product
    A <- get_constmat_tpf(knots, shape, deg)
    A_t <- t(A)

    ## calculate an inverse of A to easily produce feasible states
    A_inv <- diag(NCOL(A))
    A_inv[row(A_inv) > diff(dim(A))] <- A
    A_inv <- solve(A_inv)

    ## initialise the current estimates
    kcoef_pop <- list()
    kcoef_sub <- list()
    kpred_pop <- list()
    kpred_sub <- list()


    for (i in lvl_pop) {
        ## initialise population coefs and subjects deviations
        ## OLS as starting values
        ## kcoef_pop[[i]] <- get_ols(y_pop[[i]], X_pop[[i]])
        ## Zero vector as starting values
        kcoef_pop[[i]] <- rep(0, n_terms)
        kcoef_sub[[i]] <- matrix(0, n_terms, n_subs[[i]],
                                 dimnames = list(NULL, lvl_sub[[i]]))

        ## initialise prediction contribution by population coefs
        kpred_pop[[i]] <- X_pop[[i]] %*% kcoef_pop[[i]]

        ## initialise prediction contribution by subjects deviations
        kpred_sub[[i]] <- rep(NA, length(idx_pop[[i]]))
        for (j in lvl_sub[[i]]) {
            idx <- idx_sub[[i]][[j]]
            kpred_sub[[i]][idx] <- X_pop[[i]][idx, ] %*% kcoef_sub[[i]][, j]
        }
    }

    ## remove the dummy variables used in the for loop
    rm(i, j, idx)

    ## initialise the output list, by the order of lvl_pop and lvl_sub
    samples <- list(population = list(), subjects = list(), precision = list())
    for (i in lvl_pop) {
        ## samples$population: population coefs (num mat list)
        ## samples$subjects: subject deviations (num 3D ary list)
        samples$population[[i]] <- matrix(NA, n_terms, size)
        samples$subjects[[i]] <- array(NA, c(n_terms, n_subs[[i]], size),
                                       dimnames = list(NULL, lvl_sub[[i]]))
    }
    ## samples$precision: precisions (matrix) (num vec/mat list)
    ## pop: population coefs; poly: subject polynomial terms;
    ## sub: subject deviations; eps: error terms.
    samples$precision$pop <- rep(NA, size)
    samples$precision$poly <- array(NA, c(deg + 1, deg + 1, size))
    samples$precision$sub <- rep(NA, size)
    samples$precision$eps <- rep(NA, size)

    ## remove the dummy variables used in the for loop
    rm(i)

    ## burnin followed by actual sampling
    for (k in seq.int(-burn + 1, size)) {
        ## get the precisions (varariance-covariance matrices)
        ## if (epm_bayes) {
        ##  ##    kprecs <- list(coef_pop = coef_pop, coef_sub = coef_sub, pred_pop = pred_pop,
        ##  ## pred_sub = pred_sub)

        ## } else {
            kprecs <- get_cov_tpf(kcoef_pop, kcoef_sub, kpred_pop, kpred_sub,
                                  y_pop, idx_poly, Kmat, n_pops, n_subs,
                                  n_spline, n_samples)
        ## }
        if (k > 0) {
            ## store the precisions
            samples$precision$pop[k] <- kprecs$pop
            samples$precision$poly[, , k] <- kprecs$poly
            samples$precision$sub[k] <- kprecs$sub
            samples$precision$eps[k] <- kprecs$eps
        }

        ## get the coefs and deviations
        ## DEBUG
        ## if (k == 67) {DEBUG <- TRUE}
        ## if (k == 68) {DEBUG <- FALSE}
        kcoefs <- mapply(get_coefs_tpf, kcoef_pop, kcoef_sub, kpred_sub, X_pop,
                         X_pop_sq, X_sub_sq, lvl_sub, idx_sub,
                         y_pop,
                         MoreArgs = list(kprecs$eps, kprecs$pop,
                                         kprecs$poly, kprecs$sub, A,
                                         A_t, A_inv, Kmat, n_terms, DEBUG),
                         SIMPLIFY = FALSE)

        for (i in lvl_pop) {
            ## convert mapply output into a suitable format
            kcoef_pop[[i]] <- kcoefs[[i]]$coef_pop
            kcoef_sub[[i]] <- kcoefs[[i]]$coef_sub
            kpred_pop[[i]] <- kcoefs[[i]]$pred_pop
            kpred_sub[[i]] <- kcoefs[[i]]$pred_sub

            if (k > 0) {
                ## store the coefs and deviations
                samples$population[[i]][, k] <- kcoef_pop[[i]]
                samples$subjects[[i]][, , k] <- kcoef_sub[[i]]
            }
        }

        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }
    }

    ## return posterior means (coefs only), samples, and information regarding
    ## the basis functions.
    means <- list(population = lapply(samples$population, rowMeans),
                  subjects = lapply(samples$subjects,
                                    function(x) rowMeans(x, dims = 2)))
    basis <- list(type = "tpf", knots = knots, degree = deg)
    info <- list(lvl_pop = lvl_pop, lvl_sub = lvl_sub, n_terms = n_terms)
    data <- data.frame(x = x, y = y, grp_sub = grp_sub, grp_pop = grp_pop)

    list(means = means, samples = samples, basis = basis, info = info,
         data = data)

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
## prc_poly: precision of the individual polynomial term (num mat)
## prc_sub: precision of the individual spline terms (num)

## constants
## A: constraint matrix (num mat)
## A_t: t(A), for the purpose of using Berwin's sampler (num mat)
## A_inv: pseudo-inverse of A for generating feasible starting values (num mat)
## Kmat: to extract spline terms (num mat)
## n_terms: number of parameters to descrive a curve (num)

## RETURN
## coef_pop: population coefs (num vec)
## coef_sub: individual deviations (num mat)
## pred_pop: prediction contribution by the pop coefs (num vec)
## pred_sub: prediction contribution by the sub deviations (num vec)
get_coefs_tpf <- function(coef_pop, coef_sub, pred_sub, X_pop, X_pop_sq, X_sub_sq,
                          lvl_sub, idx_sub, y_pop,
                          prc_eps, prc_pop, prc_poly, prc_sub, A,
                          A_t, A_inv, Kmat, n_terms, DEBUG = FALSE) {

    M_pop <- solve(X_pop_sq + (prc_pop / prc_eps) * Kmat)
    mu_pop <- tcrossprod(M_pop, X_pop) %*% (y_pop - pred_sub)
    sig_pop <- M_pop / prc_eps

    lower_pop <- A %*% coef_sub
    lower_pop <- -1 * apply(lower_pop, 1, min)
    lower_pop[lower_pop < 0] <- 0

    ## Initialise the starting values of the truncated normal sampler
    start_pop <- A_inv %*% c(mu_pop[1], (lower_pop + 0.1))

    coef_pop <- start_pop
    if (DEBUG) {browser()}

    ## TRUNCATED NORMAL WITH BERWIN'S SAMPLER
    ## coef_pop <- TruncatedNormal::rmvtgauss.lin(1, mu_pop, sig_pop,
    ##                                         Amat = A_t,
    ##                                         Avec = lower_pop,
    ##                                         ## start = start_pop,
    ##                                         ## burnin = 1000)
    ##                                         start = coef_pop,
    ##                                         burnin = 500)

    ## ORDINARY NORMAL WITH BUILT IN SAMPLER
    ## coef_pop <- t(mvtnorm::rmvnorm(1, mu_pop, sig_pop))

    ## TRUNCATED NORMAL WITH MODIFIED HMC
    coef_pop <- t(tnorm::rmvtnorm(1, mu_pop, sig_pop,
                                F = A,
                                g = -lower_pop,
                                ## start = start_pop,
                                ## burnin = 1000)
                                initial = coef_pop,
                                burn = 20))

    ## tmp <- TruncatedNormal::rmvtgauss.lin(500, mu_pop, sig_pop,
    ##                                         Amat = A_t,
    ##                                         Avec = lower_pop,
    ##                                         ## start = start_pop,
    ##                                         ## burnin = 1000)
    ##                                         start = coef_pop,
    ##                                         burnin = 0)

    ## tmp <- t(tnorm::rmvtnorm(500, mu_pop, sig_pop,
    ##                             F = A,
    ##                             g = -lower_pop,
    ##                             ## start = start_pop,
    ##                             ## burnin = 1000)
    ##                             initial = coef_pop,
    ##                             burn = 0))


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
    half_N[1:NROW(prc_poly), 1:NCOL(prc_poly)] <- prc_poly
    half_N <- half_N / prc_eps

    for (j in lvl_sub) {
        idx <- idx_sub[[j]]
        M_sub <- solve(X_sub_sq[, , j] + half_N)
        mu_sub <- tcrossprod(M_sub, X_pop[idx, ]) %*% y_diff_pop[idx]
        sig_sub <- M_sub / prc_eps

        ## TRUNCATED NORMAL WITH BERWIN'S SAMPLER
        ## tmp <- TruncatedNormal::rmvtgauss.lin(500, mu_sub, sig_sub,
        ##                                       Amat = A_t,
        ##                                       Avec = lower_sub,
        ##                                       ## start = zeros,
        ##                                       ## burnin = 1000)
        ##                                       start = coef_sub[, j],
        ##                                       burnin = 0)

        ## ORDINARY NORMAL WITH BUILT IN SAMPLER
        ## coef_sub[, j] <- mvtnorm::rmvnorm(1, mu_sub, sig_sub)

        ## TRUNCATED NORMAL WITH MODIFIED HMC
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
         pred_sub = pred_sub)
}

## HMC version of get_coefs
get_coefs_hmc <- function(coef_pop, coef_sub, pred_sub, X_pop, X_pop_sq, X_sub_sq,
                          lvl_sub, idx_sub, y_pop,
                          prc_eps, prc_pop, prc_poly, prc_sub, A,
                          A_t, A_inv, Kmat, n_terms, DEBUG = FALSE) {

    ## design matrix
    C_diag <- lapply(idx_sub, function(x) X_pop[x, ])
    C_diag <- DiagMat(C_diag)
    C <- cbind(X_pop, C_diag)

    ## penalty variance
    D_top <- Kmat * prc_pop
    D_btm <- DiagMat(DiagMat(list(prc_poly, diag(prc_sub, sum(Kmat)))), length(idx_sub))
    D_inv <- prc_eps * DiagMat(list(D_top, D_btm))

    ## constraint matrix
    A_big <- DiagMat(A, length(lvl_sub) + 1)
    A_big[, 1:NCOL(A)] <- matrix(rep(t(A), length(lvl_sub) + 1), ncol = NCOL(A), byrow = TRUE)


    ## lower bound of the constraint
    lower <- rep(0, times = NROW(A_big))


    ## initialise the starting values of the truncated normal sampler
    if (any(A %*% coef_pop <= 0)) {
        coefs <- rep(A_inv %*% (lower[1:NCOL(A)] + 1), length(lvl_sub) + 1)
        ## print("no satisfaction.")
    } else {
        coefs <- c(coef_pop, coef_sub)
    }

    ## calculate "posterior" mean and variance
    M <- solve(crossprod(C) + D_inv)
    mu <- tcrossprod(M, C) %*% y_pop
    sig <- M / prc_eps

    ## TRUNCATED NORMAL WITH MODIFIED HMC

    ## coefs <- t(tnorm::rmvtnorm(1, mu, sig, F = NULL, g = NULL,
    ##                           ## start = start_pop,
    ##                           ## burnin = 1000)
    ##                           initial = coefs,
    ##                           burn = 20))

    coefs <- t(tnorm::rmvtnorm(1, mu, sig, F = A_big, g = -lower,
                              ## start = start_pop,
                              ## burnin = 1000)
                              initial = coefs,
                              burn = 20))

    coef_pop <- coefs[1:n_terms]
    coef_sub <- matrix(coefs[-(1:n_terms)], n_terms, length(lvl_sub),
                        dimnames = list(NULL, lvl_sub))


    ## Update prediction contribution by the population curve.
    pred_pop <- X_pop %*% coef_pop
    y_diff_pop <- y_pop - pred_pop

    for (j in lvl_sub) {
        idx <- idx_sub[[j]]

        ## Update prediction contribution by subject curves
        pred_sub[idx] <- X_pop[idx, ] %*% coef_sub[, j]
    }

    loglike <- mvtnorm::dmvnorm(t(coefs), mu, sig, log = TRUE)

    list(coef_pop = coef_pop, coef_sub = coef_sub, pred_pop = pred_pop,
         pred_sub = pred_sub, loglike = loglike)
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
## idx_poly: indices to extract polynomial terms (num vec)
## Kmat: to extract spline terms (num mat)
## n_pops: number of populations (num)
## n_subs: number of subjects (num vec)
## n_spline: number of spline terms/knots (num)
## n_samples: total sample size (num)

## RETURN
## prc_eps: precision of the Gaussian noise (num)
## prc_pop: precision of the population spline terms (num)
## prc_poly: precision of the individual polynomial term (num mat)
## prc_sub: precision of the individual spline terms (num)
get_cov_tpf <- function(coef_pop, coef_sub, pred_pop, pred_sub, y_pop,
                        idx_poly, Kmat, n_pops, n_subs, n_spline, n_samples) {

    ## hyperparameters of priors
    ig_a <- 0.001
    ig_b <- 0.001
    wi_df <- length(idx_poly)
    wi_sig <- diag(0.0011, length(idx_poly))

    ## if dealing with a single population model, convert everything to lists.
    if (typeof(coef_pop) == "list") coef_pop else list(coef_pop)
    if (typeof(coef_sub) == "list") coef_sub else list(coef_sub)
    if (typeof(pred_pop) == "list") pred_pop else list(pred_pop)
    if (typeof(pred_sub) == "list") pred_sub else list(pred_sub)

    ## crossprod of the spline terms
    xspl_pop <- sum(vapply(coef_pop, function(x) crossprod(Kmat %*% x), 0))
    xspl_sub <- sum(vapply(coef_sub, function(x) crossprod(c(Kmat %*% x)), 0))

    ## precision of population spline terms
    shp_pop <- (n_pops * n_spline) / 2 + ig_a
    scl_pop <- 0.5 * xspl_pop + ig_b
    prc_pop <- rgamma(1, shape = shp_pop, rate = scl_pop)

    ## precision of individual spline terms
    shp_sub <- (sum(n_subs) * n_spline) / 2 + ig_a
    scl_sub <- 0.5 * xspl_sub + ig_b
    prc_sub <- rgamma(1, shape = shp_sub, rate = scl_sub)

    ## precision of residuals
    resid_vec <- mapply(function(x, y, z) x - y - z,
                        y_pop, pred_pop, pred_sub,
                        SIMPLIFY = FALSE)
    shp_eps <- 0.5 * n_samples + ig_a
    scl_eps <- 0.5 * crossprod(unlist(resid_vec)) + ig_b
    prc_eps <- rgamma(1, shape = shp_eps, rate = scl_eps)

    ## tcrossprod of the polynomial terms, presented as a 3D array
    txpoly <- vapply(coef_sub, function(x) tcrossprod(x[idx_poly, , drop = FALSE]),
                     diag(as.double(idx_poly)))

    ## if there is only one row in coef_sub (i.e. model only with a random intercept)
    if (!is.array(txpoly)) {
        dim(txpoly) <- c(1, 1, 1)
    }

    ## precision of individual polynomial terms
    df_poly <- wi_df + sum(n_subs)
    Sigma_poly <- solve(wi_sig + rowSums(txpoly, dims = 2))
    prc_poly <- rWishart(1, df = df_poly, Sigma = Sigma_poly)[, , 1]

    list(pop = prc_pop, sub = prc_sub, poly = prc_poly, eps = prc_eps)
}


## Fit a penalised spline with nlme. A random effect is added to each population
## coefficient.
## data : 1st col x, 2nd col y, 3rd col groups
## K : number of quantile (inner) knots
## deg : degree of spline polynomial
lme_tpf <- function(data, K, deg = 1) {
    x <- data[[1]]
    y <- data[[2]]

    ## convert the group variable into a factor
    if (is.factor(data[[3]])) {
        grp <- droplevels(data[[3]])
    } else {
        grp <- factor(data[[3]], levels = unique(data[[3]]))
    }

    ## dummy variable of ones
    ones <- rep(1, length(grp))

    ## design matrices
    design_ls <- get_design_tpf(x, K, deg)
    X <- design_ls$design[, 1:(deg + 1)]
    Z <- design_ls$design[, (deg + 2):(deg + 1 + K)]

    ## covariance structures of random effects
    pop_pd <- nlme::pdIdent(~ Z - 1)
    sub_pd <- nlme::pdBlocked(list(nlme::pdSymm(~ X - 1), nlme::pdIdent(~ Z - 1)))

    fm <- nlme::lme(fixed = y ~ X - 1, random = list(ones = pop_pd, grp = sub_pd),
                    control = list(maxIter = 50, msMaxIter = 150, niterEM = 150))

    ## extract coefficients
    pop_coef <- as.numeric(coef(fm, level = 1))
    sub_coef <- t(as.matrix(coef(fm, level = 2))) - pop_coef
    colnames(sub_coef) <- levels(grp)
    rownames(sub_coef) <- NULL

    means <- list(population = pop_coef, subjects = sub_coef)
    basis <- list(type = "tpf", knots = design_ls$knots, degree = deg)
    info <- list(lvl_pop = NULL, lvl_sub = levels(grp), n_terms = deg + 1 + K)
    data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)

    list(means = means, samples = NULL, basis = basis, info = info,
         data = data, model = fm)
}


## Fit a penalised spline with nlme. A random effect is added to each population
## coefficient. Different knots on population and subject-specific levels.
## data : 1st col x, 2nd col y, 3rd col groups
## K_pop : number of quantile (inner) knots
## K_sub : number of quantile (inner) knots
## deg : degree of spline polynomial
lme_tpf <- function(data, K, deg = 1) {
    x <- data[[1]]
    y <- data[[2]]

    ## convert the group variable into a factor
    if (is.factor(data[[3]])) {
        grp <- droplevels(data[[3]])
    } else {
        grp <- factor(data[[3]], levels = unique(data[[3]]))
    }

    ## dummy variable of ones
    ones <- rep(1, length(grp))

    ## design matrices for the population curve
    design_ls <- get_design_tpf(x, K, deg)
    X <- design_ls$design[, 1:(deg + 1)]
    Z <- design_ls$design[, (deg + 2):(deg + 1 + K)]

    ## covariance structures of random effects
    pop_pd <- nlme::pdIdent(~ Z - 1)
    sub_pd <- nlme::pdBlocked(list(nlme::pdSymm(~ X - 1), nlme::pdIdent(~ Z - 1)))

    fm <- nlme::lme(fixed = y ~ X - 1, random = list(ones = pop_pd, grp = sub_pd),
                    control = list(maxIter = 50, msMaxIter = 150, niterEM = 150))

    ## extract coefficients
    pop_coef <- as.numeric(coef(fm, level = 1))
    sub_coef <- t(as.matrix(coef(fm, level = 2))) - pop_coef
    colnames(sub_coef) <- levels(grp)
    rownames(sub_coef) <- NULL

    means <- list(population = pop_coef, subjects = sub_coef)
    basis <- list(type = "tpf", knots = design_ls$knots, degree = deg)
    info <- list(lvl_pop = NULL, lvl_sub = levels(grp), n_terms = deg + 1 + K)
    data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)

    list(means = means, samples = NULL, basis = basis, info = info,
         data = data, model = fm)
}






## fit_
##     pop_pd <- nlme::pdIdent(~ Z.q - 1)
##     ## sub.pd.q <- pdBlocked(list(pdSymm(~ time.q + I(time.q^2)), pdIdent(~ Z.q - 1)))
##     ## sub.pd.q <- pdSymm(~ time.q + I(time.q^2))
##     sub_pd <- pdBlocked(list(pdSymm(~ 1), pdIdent(~ Z.q - 1)))
##     quad_fm <- lme(fixed = y ~ time.q + I(time.q^2),
##                    random = list(pop.level = pop.pd.q, sub.level = sub.pd.q))


