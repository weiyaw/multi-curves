## source("subor.R")
## This is a Gibbs sampler for Bayesian ridge; see thesis.
## Requirements: Bmat, y, Kmat
## Algorithms paremeters: burn, size
## Extras: verbose
bayes_ridge <- function(y, Bmat, Kmat, burn, size, start = NULL, verbose = TRUE) {

    ## hyperparemeters for priors
    ig_a_pop <- 0.001                        # flat prior on standard deviations
    ig_b_pop <- 0.001
    ig_a_eps <- 0.001                        # flat prior on standard deviations
    ig_b_eps <- 0.001

    ## some precalculation
    rank_K <- NROW(Kmat)
    n_terms <- NCOL(Bmat)
    n_samples <- NROW(Bmat)
    xBmat <- crossprod(Bmat)
    xKmat <- crossprod(Kmat)
    Bxy <- crossprod(Bmat, y)

    ## initialise the output list
    samples <- list(population = matrix(NA, n_terms, size),
                    precision = list())

    if (is.null(start)) {
        ## initialise theta using ols
        kcoef_pop <- tcrossprod(solve(crossprod(Bmat)), Bmat) %*% as.vector(y)
    } else {
        kcoef_pop <- as.vector(start)
    }

    for (k in seq.int(-burn + 1, size)) {

        ## update sigma^2_theta
        shape_pop <- 0.5 * rank_K + ig_a_pop
        rate_pop <- 0.5 * crossprod(Kmat %*% kcoef_pop) + ig_b_pop
        kprec_pop <- rgamma(1, shape = shape_pop, rate = rate_pop)

        ## update sigma^2_epsilon
        shape_eps <- 0.5 * n_samples + ig_a_eps
        rate_eps <- 0.5 * crossprod(y - Bmat %*% kcoef_pop) + ig_b_eps
        kprec_eps <- rgamma(1, shape = shape_eps, rate = rate_eps)

        ## update theta
        M <- chol2inv(chol(xBmat + kprec_pop / kprec_eps * xKmat))
        mu <- M %*% Bxy
        sig <- 1 / kprec_eps * M
        kcoef_pop <- t(mvtnorm::rmvnorm(1, mu, sig))

        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }

        if (k > 0) {
            samples$precision$pop[k] <- kprec_pop
            samples$precision$eps[k] <- kprec_eps
            samples$population[, k] <- kcoef_pop
        }
    }
    means <- list(population = rowMeans(samples$population))
    list(means = means, samples = samples)
}



## This is a Gibbs sampler for longitudinal Bayesian ridge; see thesis.
## Requirements: Bmat, y, grp, Kmat, dim_sub1
## Algorithms paremeters: burn, size
## Extras: verbose
bayes_ridge_sub <- function(y, grp, Bmat, Kmat, dim_sub1, burn, size, start_pop = NULL,
                            verbose = TRUE) {

    ## hyperparemeters for priors
    ig_a_pop <- 0.001
    ig_b_pop <- 0.001
    ig_a_sub2 <- 0.001
    ig_b_sub2 <- 0.001
    ig_a_eps <- 0.001
    ig_b_eps <- 0.001
    iw_v <- dim_sub1 + 1
    iw_lambda <- diag(dim_sub1)          # assuming symmetric

    ## clean up grp variable
    if (is.factor(grp)) {
        grp <- droplevels(grp)
    } else {
        grp <- factor(grp, levels = unique(grp))
    }

    ## some precalculation
    rank_K <- NROW(Kmat)
    n_terms <- NCOL(Bmat)
    n_samples <- NROW(Bmat)
    n_subs <- length(unique(grp))
    idx <- tapply(seq_len(n_samples), grp, function(x) x)

    ## some crossproducts precalculation
    xBmat <- crossprod(Bmat)
    xBmat_sub <- array(NA, c(n_terms, n_terms, n_subs), list(NULL, NULL, levels(grp)))
    for (i in levels(grp)) {
        xBmat_sub[, , i] <- crossprod(Bmat[idx[[i]], ])
    }
    xKmat <- crossprod(Kmat)

    ## initialise the output list
    samples <- list(population = matrix(NA, n_terms, size),
                    subjects = array(NA, c(n_terms, n_subs, size),
                                     dimnames = list(NULL, levels(grp))),
                    precision = list(pop = rep(NA, size),
                                     sub1 = array(NA, c(dim_sub1, dim_sub1, size)),
                                     sub2 = rep(NA, size),
                                     eps = rep(NA, size)))

    ## initialise theta with an LS estimate, and delta with rnorms with small sd
    if (is.null(start_pop)) {
        kcoef_pop <- tcrossprod(solve(crossprod(Bmat)), Bmat) %*% as.vector(y)
    } else {
        kcoef_pop <- as.vector(start_pop)
    }
    kcoef_sub <- matrix(rnorm(n_terms * n_subs) * 0.01, n_terms, n_subs,
                        dimnames = list(NULL, levels(grp)))

    ## initialise prediction contribution by population coefs and subjects deviations
    kcontrib_pop <- Bmat %*% kcoef_pop
    kcontrib_sub <- rep(NA, n_samples)
    for (i in levels(grp)) {
        kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef_sub[, i]
    }

    for (k in seq.int(-burn + 1, size)) {

        ## update sigma^2_theta
        shape_pop <- 0.5 * rank_K + ig_a_pop
        rate_pop <- 0.5 * crossprod(Kmat %*% kcoef_pop) + ig_b_pop
        kprec_pop <- rgamma(1, shape = shape_pop, rate = rate_pop)

        ## update sigma^2_epsilon
        shape_eps <- 0.5 * n_samples + ig_a_eps
        rate_eps <- 0.5 * crossprod(y - kcontrib_pop - kcontrib_sub) + ig_b_eps
        kprec_eps <- rgamma(1, shape = shape_eps, rate = rate_eps)

        ## update Sigma_dev1
        df_sub1 <- iw_v + n_subs
        scale_sub1 <- iw_lambda + tcrossprod(kcoef_sub[1:dim_sub1, , drop = FALSE])
        inv_scale_sub1 <- chol2inv(chol(scale_sub1))
        kprec_sub1 <- rWishart(1, df = df_sub1, Sigma = inv_scale_sub1)[, , 1]

        ## update sigma^2_dev2
        shape_sub2 <- 0.5 * n_subs * (n_terms - dim_sub1) + ig_a_sub2
        rate_sub2 <- 0.5 * crossprod(c(kcoef_sub[-(1:dim_sub1), ])) + ig_b_sub2
        kprec_sub2 <- rgamma(1, shape = shape_sub2, rate = rate_sub2)

        ## update theta
        M_pop <- chol2inv(chol(xBmat + kprec_pop / kprec_eps * xKmat))
        mu <- M_pop %*% crossprod(Bmat, y - kcontrib_sub)
        sig <- 1 / kprec_eps * M_pop
        kcoef_pop <- t(mvtnorm::rmvnorm(1, mu, sig))
        kcontrib_pop <- Bmat %*% kcoef_pop

        ## update delta
        kprec_sub <- diag(kprec_sub2, n_terms)
        kprec_sub[1:NROW(kprec_sub1), 1:NCOL(kprec_sub1)] <- kprec_sub1

        for (i in levels(grp)) {
            M_sub <- chol2inv(chol(xBmat_sub[, , i] + kprec_sub / kprec_eps))
            y_star <- y[idx[[i]]] - kcontrib_pop[idx[[i]]]
            mu <- M_sub %*% crossprod(Bmat[idx[[i]], ], y_star)
            sig <- M_sub / kprec_eps
            kcoef_sub[, i] <- t(mvtnorm::rmvnorm(1, mu, sig))
            kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef_sub[, i]
        }

        ## print progress
        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }

        ## store samples after burn-in iterations
        if (k > 0) {
            samples$precision$pop[k] <- kprec_pop
            samples$precision$sub1[, , k] <- kprec_sub1
            samples$precision$sub2[k] <- kprec_sub2
            samples$precision$eps[k] <- kprec_eps
            samples$population[, k] <- kcoef_pop
            samples$subjects[, , k] <- kcoef_sub
        }
    }

    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    list(means = means, samples = samples)
}






