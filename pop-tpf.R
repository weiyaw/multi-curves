source("subor.R")

## Penalised splines, random effects on polynomials only, truncated power function
## data : 1st col x, 2nd col y, 3rd col groups
## K : number of quantile (inner) knots, or a vector of inner knots
## deg : degree of spline polynomial
pop_tpf <- function(data, K, deg = 1, random = deg + 1, shape = "increasing",
                    size = 100, burn = size / 10, verbose = FALSE) {

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

    n_samples <- NROW(data)
    n_terms_pop <- K + deg + 1

    if (!(random > 0 || random < deg + 2)) {
        stop("Unsupported number of random effects.")
    }
    n_terms_sub <- random

    lvl_sub <- levels(grp)
    idx_sub <- tapply(seq_len(n_samples), grp, function(x) x)
    n_subs <- length(idx_sub)

    ## construct the design matrix and knots
    design_ls <- get_design_tpf(x, K, deg)
    knots <- design_ls$knots            # all knots (with extrema)
    n_spline <- length(knots) - 2       # number of inner knots (w/o extrema)
    X_pop <- design_ls$design
    X_sub <- design_ls$design[, 1:n_terms_sub, drop = FALSE]

    ## get rid of the design_ls to save space
    rm(design_ls)

    ## construct cross-products of the design matrix
    X_pop_sq <- crossprod(X_pop)
    X_sub_sq <- array(NA, c(n_terms_sub, n_terms_sub, n_subs),
                      list(NULL, NULL, lvl_sub))
    for (j in lvl_sub) {
        X_sub_sq[, , j] <- crossprod(X_sub[idx_sub[[j]], ])
    }

    ## An idempotent matrix to extract the spline terms
    Kmat <- diag(c(rep(0, deg + 1), rep(1, n_spline)))
    idx_poly <- seq_len(n_terms_sub)        # index of polynomial terms

    ## Construct the constraint matrix and its cross-product
    A <- get_constmat_tpf(knots, shape, deg)
    A_t <- t(A)

    ## Calculate an inverse of A to easily produce feasible states
    A_inv <- diag(NCOL(A))
    A_inv[row(A_inv) > diff(dim(A))] <- A
    A_inv <- solve(A_inv)

    ## initialise population coefs and subjects deviations
    kcoef_pop <- get_ols(y, X_pop)
    kcoef_sub <- matrix(0, n_terms_sub, n_subs, dimnames = list(NULL, lvl_sub))

    ## initialise prediction contribution by population coefs
    kpred_pop <- X_pop %*% kcoef_pop

    ## initialise prediction contribution by subjects deviations
    kpred_sub <- rep(NA, n_samples)
    for (j in lvl_sub) {
        idx <- idx_sub[[j]]
        kpred_sub[idx] <- X_sub[idx, , drop = FALSE] %*% kcoef_sub[, j]
    }

    ## remove the dummy variables used in the for loop
    rm(j, idx)

    ## initialise the output list, by the order of lvl_sub
    samples <- list(population = matrix(NA, n_terms_pop, size),
                    subjects = array(NA, c(n_terms_pop, n_subs, size),
                                     dimnames = list(NULL, lvl_sub)),
                    precision = list())
    ## samples$precision: precisions (matrix) (num vec/mat list)
    ## pop: population coefs; poly: subject polynomial terms;
    ## sub: subject deviations; eps: error terms.
    samples$precision$pop <- rep(NA, size)
    samples$precision$poly <- array(NA, c(deg + 1, deg + 1, size))
    samples$precision$eps <- rep(NA, size)
    ## samples$precision$sub <- rep(NA, size)

    ## burnin followed by actual sampling
    for (k in seq.int(-burn + 1, size)) {
        ## get the precisions (varariance-covariance matrices)
        kprecs <- get_cov_tpf(list(kcoef_pop), list(kcoef_sub), list(kpred_pop),
                              list(kpred_sub), list(y), idx_poly, Kmat, 1,
                              n_subs, n_spline, n_samples)

        ## get the coefs and deviations
        kcoefs <- get_coefs_tpf(kcoef_sub, kpred_sub, X_pop, X_sub, X_pop_sq, X_sub_sq,
                                lvl_sub, idx_sub, y, kprecs$eps, kprecs$pop,
                                kprecs$poly, kprecs$sub, A, A_t, A_inv, Kmat,
                                n_terms_pop, n_terms_sub)

        ## for the ease of reading
        kcoef_pop <- kcoefs$coef_pop
        kcoef_sub <- kcoefs$coef_sub
        kpred_pop <- kcoefs$pred_pop
        kpred_sub <- kcoefs$pred_sub

        if (k > 0) {
            ## store the precisions
            samples$precision$pop[k] <- kprecs$pop
            samples$precision$poly[, , k] <- kprecs$poly
            ## samples$precision$sub[k] <- NULL
            samples$precision$eps[k] <- kprecs$eps

            ## store the coefs and deviations
            samples$population[, k] <- kcoef_pop
            samples$subjects[1:n_terms_sub, , k] <- kcoef_sub
        }

        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }
    }

    samples$subjects[(n_terms_sub + 1):n_terms_pop, , ] <- 0

    ## return posterior means (coefs only), samples, and information regarding
    ## the basis functions.
    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    basis <- list(type = "tpf", knots = knots, degree = deg)
    info <- list(lvl_pop = NULL, lvl_sub = lvl_sub, n_terms = n_terms_pop)
    data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)

    cat((A %*% means$population) > 0, '\n')
    list(means = means, samples = samples, basis = basis, info = info,
         data = data)
}

## get a sample from the coefs posterior (one population)

## different in each ITERATION and POPULATION
## coef_sub: individual deviations (num mat)
## pred_sub: prediction contribution by the sub deviations (num vec)

## different in each POPULATION
## X_pop: design matrix for population coefficients (num mat)
## X_sub: design matrix for individual coefficients (num mat)
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
get_coefs_tpf <- function(coef_sub, pred_sub, X_pop, X_sub, X_pop_sq, X_sub_sq,
                          lvl_sub, idx_sub, y_pop,
                          prc_eps, prc_pop, prc_poly, prc_sub, A,
                          A_t, A_inv, Kmat, n_terms_pop, n_terms_sub,
                          DEBUG = FALSE) {

    M_pop <- solve(X_pop_sq + (prc_pop / prc_eps) * Kmat)
    mu_pop <- tcrossprod(M_pop, X_pop) %*% (y_pop - pred_sub)
    sig_pop <- M_pop / prc_eps

    lower_pop <- A[, 1:n_terms_sub, drop = FALSE] %*% coef_sub
    lower_pop <- -1 * apply(lower_pop, 1, min)
    lower_pop[lower_pop < 0] <- 0

    ## Initialise the starting values of the truncated normal sampler
    start_pop <- A_inv %*% c(mu_pop[1], (lower_pop + 0.1))

    coef_pop <- start_pop

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
    ## zeros <- rep(0, n_terms_sub)
    coef_sub <- matrix(0, n_terms_sub, length(lvl_sub),
                       dimnames = list(NULL, lvl_sub))


    ## Calculate the precision matrix term
    ## half_N <- diag(prc_sub, n_terms)
    half_N <- prc_poly
    half_N <- half_N / prc_eps

    for (j in lvl_sub) {
        idx <- idx_sub[[j]]
        M_sub <- solve(X_sub_sq[, , j] + half_N)
        mu_sub <- tcrossprod(M_sub, X_sub[idx, , drop = FALSE]) %*% y_diff_pop[idx]
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
                                         F = A[, 1:n_terms_sub, drop = FALSE],
                                         g = -lower_sub,
                                         ## start = zeros,
                                         ## burnin = 1000)
                                         initial = coef_sub[, j],
                                         burn = 20)


        ## Update prediction contribution by subject curves
        pred_sub[idx] <- X_sub[idx, , drop = FALSE] %*% coef_sub[, j]
    }

    list(coef_pop = coef_pop, coef_sub = coef_sub, pred_pop = pred_pop,
         pred_sub = pred_sub)
}

## get a sample from the precision posterior (one/multiple populations)

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
    wi_df <- 1
    wi_sig <- diag(0.001, length(idx_poly))

    ## if dealing with a single population model, convert everything to lists.
    if (typeof(coef_pop) == "list") coef_pop else list(coef_pop)
    if (typeof(coef_sub) == "list") coef_sub else list(coef_sub)
    if (typeof(pred_pop) == "list") pred_pop else list(pred_pop)
    if (typeof(pred_sub) == "list") pred_sub else list(pred_sub)

    ## crossprod of the spline terms
    xspl_pop <- sum(vapply(coef_pop, function(x) crossprod(Kmat %*% x), 0))
    ## xspl_sub <- sum(vapply(coef_sub, function(x) crossprod(c(Kmat %*% x)), 0))

    ## precision of population spline terms
    shp_pop <- (n_pops * n_spline) / 2 + ig_a
    scl_pop <- 0.5 * xspl_pop + ig_b
    prc_pop <- rgamma(1, shape = shp_pop, rate = scl_pop)

    ## precision of individual spline terms
    ## shp_sub <- (sum(n_subs) * n_spline) / 2 + ig_a
    ## scl_sub <- 0.5 * xspl_sub + ig_b
    ## prc_sub <- rgamma(1, shape = shp_sub, rate = scl_sub)

    ## precision of residuals
    resid_vec <- mapply(function(x, y, z) x - y - z,
                        y_pop, pred_pop, pred_sub,
                        SIMPLIFY = FALSE)
    shp_eps <- 0.5 * n_samples + ig_a
    scl_eps <- 0.5 * crossprod(unlist(resid_vec)) + ig_b
    prc_eps <- rgamma(1, shape = shp_eps, rate = scl_eps)

    ## tcrossprod of the polynomial terms, presented as a 3D array
    txpoly <- vapply(coef_sub, function(x) tcrossprod(x[idx_poly, , drop = FALSE]),
                     diag(length(idx_poly)))

    ## if there is only one row in coef_sub (i.e. model only with a random intercept)
    if (!is.array(txpoly)) {
        dim(txpoly) <- c(1, 1, 1)
    }

    ## precision of individual polynomial terms
    df_poly <- wi_df + sum(n_subs)
    Sigma_poly <- solve(wi_sig + rowSums(txpoly, dims = 2))
    prc_poly <- rWishart(1, df = df_poly, Sigma = Sigma_poly)[, , 1]

    list(pop = prc_pop, sub = NULL, poly = prc_poly, eps = prc_eps)
}

## Get an ordinary least squares estimate with a given response and design
## matrix.
get_ols <- function(response, design) {
    tcrossprod(solve(crossprod(design)), design) %*% as.vector(response)
}


