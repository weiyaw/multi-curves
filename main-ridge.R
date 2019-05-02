block_diag <- function(..., size = NULL) {

    ## Construct a big matrix with its diagonal elements the matrices provided
    ## in "A". If "A" is a matrix, the function returns a matrix with "size" of
    ## "A" in the diagonal elements. If "A" is a list, it returns a matrix with
    ## diagonal elements as all the matrices contained inside the list.

    A <- list(...)
    if (length(A) == 1 && is.matrix(A[[1]])) {
        if (is.null(size)) {
            stop("Size of the resulting matrix not supplied.")
        }
        n_row <- NROW(A[[1]])
        n_col <- NCOL(A[[1]])
        row_idx <- seq(1, n_row)
        col_idx <- seq(1, n_col)
        res <- matrix(0, size * n_row, size * n_col)
        for (i in seq(0, size - 1)) {
            res[i * n_row + row_idx, i * n_col + col_idx] <- A[[1]]
        }
        return(res)

   } else if (length(A) > 1) {
       if (!all(vapply(A, is.matrix, TRUE))) {
            stop("The list contains non-matrix objects.")
       }
        dims <- vapply(A, dim, c(1, 1))
        total_dims <- rowSums(dims)
        res <- matrix(0, total_dims[1], total_dims[2])
        row_rolling <- col_rolling <- 0
        for (i in seq(1, NCOL(dims))) {
            row_idx <- row_rolling + seq(1, dims[1, i])
            col_idx <- col_rolling + seq(1, dims[2, i])
            res[row_idx, col_idx] <- A[[i]]
            row_rolling <- row_rolling + dims[1, i]
            col_rolling <- col_rolling + dims[2, i]
        }
        return(res)

   } else {
       warning("Non-matrix or list object supplied")
   }
}

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
        kcoef_pop <- tcrossprod(solve(crossprod(Bmat) + xKmat), Bmat) %*% as.vector(y)
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
        kprec_sub <- block_diag(kprec_sub1, diag(kprec_sub2, n_terms - dim_sub1))

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


## This is a Gibbs sampler for constrained longitudinal Bayesian ridge; see thesis.
## A * theta >= 0 ; A * (theta + delta_i) >= 0
## The vector of constriant has to be zero.
## Assumption: Full row rank for Amat.
## Requirements: Bmat, y, grp, Kmat, dim_sub1, Amat
## Algorithms paremeters: burn, size
## Extras: verbose
bayes_ridge_cons_sub <- function(y, grp, Bmat, Kmat, dim_sub1, Amat, burn, size,
                                 start_pop = NULL, verbose = TRUE) {

    ## hyperparemeters for priors
    ig_a_pop <- 0.001
    ig_b_pop <- 0.001
    ig_a_sub1 <- 0.001
    ig_b_sub1 <- 0.001
    ig_a_sub2 <- 0.001
    ig_b_sub2 <- 0.001
    ig_a_eps <- 0.001
    ig_b_eps <- 0.001
    omega_1 <- diag(dim_sub1)              # positive definite
    omega_2 <- diag(NCOL(Bmat) - dim_sub1) # positive definite

    if (NCOL(Amat) != NCOL(Bmat)) {
        stop("Dims of Amat and Bmat inconsistent.")
    }

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
    inv_omega_1 <- chol2inv(chol(omega_1))
    inv_omega_2 <- chol2inv(chol(omega_2))

    ## some crossproducts precalculation
    xBmat <- crossprod(Bmat)
    xBmat_sub <- array(NA, c(n_terms, n_terms, n_subs), list(NULL, NULL, levels(grp)))
    for (i in levels(grp)) {
        xBmat_sub[, , i] <- crossprod(Bmat[idx[[i]], ])
    }
    xKmat <- crossprod(Kmat)

    ## a matrix to generate feasible starting points quickly
    Ainv <- diag(NCOL(Amat))
    Ainv[row(Ainv) > diff(dim(Amat))] <- Amat
    Ainv <- solve(Ainv)
    gen_init <- function(lower, mu) {
        Ainv %*% c(mu[1:(NCOL(Ainv) - length(lower))], lower + 1)
    }

    ## initialise the output list
    samples <- list(population = matrix(NA, n_terms, size),
                    subjects = array(NA, c(n_terms, n_subs, size),
                                     dimnames = list(NULL, levels(grp))),
                    precision = list(pop = rep(NA, size),
                                     sub1 = array(NA, c(dim_sub1, dim_sub1, size)),
                                     sub2 = rep(NA, size),
                                     eps = rep(NA, size)))

    ## initialise theta and delta with tnorms with small sd
    if (is.null(start_pop)) {
        kcoef_pop <- Ainv %*% rep(1, ncol(Amat))
    } else {
        kcoef_pop <- as.vector(start_pop)
    }
    kcoef_sub <- t(tnorm::rmvtnorm(n_subs, kcoef_pop, diag(n_terms) * 0.01,
                                   initial = gen_init(-Amat %*% kcoef_pop, kcoef_pop),
                                   F = Amat,
                                   g = Amat %*% kcoef_pop))
    colnames(kcoef_sub) <- levels(grp)

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

        ## update sigma^2_dev1
        shape_sub1 <- 0.5 * n_subs * dim_sub1 + ig_a_sub1
        quad_b_omega <- apply(kcoef_sub[1:dim_sub1, ], 2,
                              function(x) crossprod(x, inv_omega_1 %*% x))
        rate_sub1 <- 0.5 * sum(quad_b_omega) + ig_b_sub1
        kprec_sub1 <- rgamma(1, shape = shape_sub1, rate = rate_sub1)

        ## update sigma^2_dev2
        shape_sub2 <- 0.5 * n_subs * (n_terms - dim_sub1) + ig_a_sub2
        quad_v_omega <- apply(kcoef_sub[-(1:dim_sub1), ], 2,
                              function(x) crossprod(x, inv_omega_2 %*% x))
        rate_sub2 <- 0.5 * sum(quad_v_omega) + ig_b_sub2
        kprec_sub2 <- rgamma(1, shape = shape_sub2, rate = rate_sub2)

        ## update theta
        lower_pop <- apply(-Amat %*% kcoef_sub, 1, max)
        lower_pop[lower_pop < 0] <- 0
        M_pop <- chol2inv(chol(xBmat + kprec_pop / kprec_eps * xKmat))
        mu <- M_pop %*% crossprod(Bmat, y - kcontrib_sub)
        sig <- 1 / kprec_eps * M_pop
        ## kcoef_pop <- t(mvtnorm::rmvnorm(1, mu, sig))
        kcoef_pop <- t(tnorm::rmvtnorm(1, mu, sig,
                                       ## initial = gen_init(lower_pop, mu),
                                       initial = kcoef_pop,
                                       F = Amat,
                                       g = -lower_pop))
        kcontrib_pop <- Bmat %*% kcoef_pop

        ## update delta
        lower_sub <- -Amat %*% kcoef_pop
        kprec_sub <- block_diag(kprec_sub1 * inv_omega_1, kprec_sub2 * inv_omega_2)

        for (i in levels(grp)) {
            M_sub <- chol2inv(chol(xBmat_sub[, , i] + kprec_sub / kprec_eps))
            y_star <- y[idx[[i]]] - kcontrib_pop[idx[[i]]]
            mu <- M_sub %*% crossprod(Bmat[idx[[i]], ], y_star)
            sig <- M_sub / kprec_eps
            ## kcoef_sub[, i] <- t(mvtnorm::rmvnorm(1, mu, sig))
            kcoef_sub[, i] <- t(tnorm::rmvtnorm(1, mu, sig,
                                                ## initial = gen_init(lower_sub, mu),
                                                initial = kcoef_sub[, i],
                                                F = Amat,
                                                g = -lower_sub))
            kcontrib_sub[idx[[i]]] <- Bmat[idx[[i]], ] %*% kcoef_sub[, i]
        }

        ## print progress
        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }

        ## store samples after burn-in iterations
        if (k > 0) {
            samples$precision$pop[k] <- kprec_pop
            samples$precision$sub1[, , k] <- kprec_sub1 * inv_omega_1
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

