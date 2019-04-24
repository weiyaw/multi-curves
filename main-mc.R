## This is a MC sampler that takes precision samples and gives coefs
## Requirements: Bmat, y, grp, Kmat, prec
## Algorithms paremeters: burn, size
## Extras: verbose
mc_sub <- function(y, grp, Bmat, Kmat, prec, size, verbose = TRUE) {
    prec_pop <- prec$pop
    prec_sub1 <- prec$sub1
    prec_sub2 <- prec$sub2
    prec_eps <- prec$eps

    n_terms <- NCOL(Bmat)
    n_samples <- NROW(Bmat)
    n_subs <- length(unique(grp))
    dim_sub1 <- NROW(prec_sub1[, , 1])
    idx <- tapply(seq_len(n_samples), grp, function(x) x)

    ## some crossproducts precalculation
    xBmat <- crossprod(Bmat)
    xBmat_sub <- array(NA, c(n_terms, n_terms, n_subs), list(NULL, NULL, levels(grp)))
    Bxy_sub <- matrix(NA, n_terms, n_subs, dimnames = list(NULL, levels(grp)))
    for (i in levels(grp)) {
        xBmat_sub[, , i] <- crossprod(Bmat[idx[[i]], ])
        Bxy_sub[, i] <- crossprod(Bmat[idx[[i]], ], y[idx[[i]]])
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

    ## initialise some intermediate output
    BMB <- array(NA, c(n_terms, n_terms, n_subs), list(NULL, NULL, levels(grp)))
    BMy <- matrix(NA, n_terms, n_subs, dimnames = list(NULL, levels(grp)))
    kcoef_sub <- matrix(NA, n_terms, n_subs, dimnames = list(NULL, levels(grp)))

    ## get indices of the resamped variances
    resampling_idx <- sample.int(length(prec_pop), size, TRUE)
    ## resampling_idx <- (1+100):(size + 100)

    for (k in 1:size) {

        ## get precisions
        r <- resampling_idx[k]
        kprec_pop <- prec_pop[r]
        kprec_eps <- prec_eps[r]
        kprec_sub <- diag(prec_sub2[r], n_terms)
        kprec_sub[1:dim_sub1, 1:dim_sub1] <- prec_sub1[, , r]

        ## get theta
        for (i in levels(grp)) {
            ## for numerical stability, these steps are simplified
            xBmat_i <- xBmat_sub[, , i]
            Li <- xBmat_i + kprec_sub / kprec_eps
            inv_Li <- chol2inv(chol(Li))
            BMB[, , i] <- kprec_eps * xBmat_i - xBmat_i %*% inv_Li %*% xBmat_i

            Bxy_i <- Bxy_sub[, i]
            BMy[, i] <- kprec_eps * Bxy_i - xBmat_i %*% inv_Li %*% Bxy_i
        }
        Phi <- kprec_pop * xKmat + rowSums(BMB, dims = 2)
        inv_Phi <- chol2inv(chol(Phi))
        kcoef_pop <- t(mvtnorm::rmvnorm(1, inv_Phi %*% rowSums(BMy), inv_Phi))
        kcontrib_pop <- Bmat %*% kcoef_pop

        ## get delta
        for (i in levels(grp)) {
            M_sub <- chol2inv(chol(xBmat_sub[, , i] + kprec_sub / kprec_eps))
            y_star <- y[idx[[i]]] - kcontrib_pop[idx[[i]]]
            mu <- M_sub %*% crossprod(Bmat[idx[[i]], ], y_star)
            sig <- M_sub / kprec_eps
            kcoef_sub[, i] <- t(mvtnorm::rmvnorm(1, mu, sig))
        }

        ## print progress
        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }

        ## store samples
        samples$precision$pop[k] <- kprec_pop
        samples$precision$sub1[, , k] <- prec_sub1[, , r]
        samples$precision$sub2[k] <- prec_sub2[r]
        samples$precision$eps[k] <- kprec_eps
        samples$population[, k] <- kcoef_pop
        samples$subjects[, , k] <- kcoef_sub
    }
    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    list(means = means, samples = samples)
}
