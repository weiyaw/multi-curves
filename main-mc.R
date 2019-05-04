block_diag <- function(..., size = NULL) {

    ## Construct a big matrix with its diagonal elements the matrices provided
    ## in "...". If there is only one argument in "...", the function returns a
    ## matrix with multiple ("size") of "..." in the diagonal elements. If there
    ## are more than one argument in "...", it returns a matrix with diagonal
    ## elements as all the matrices contained in "...".

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


## This is a MC sampler that takes precision samples and gives coefs
## The size of coefs equals to the size of prec
## Requirements: Bmat, y, grp, Kmat, prec
## Algorithms paremeters: burn
## Extras: verbose
mc_sub <- function(y, grp, Bmat, Kmat, prec, verbose = TRUE) {
    prec_pop <- prec$pop
    prec_sub1 <- prec$sub1
    prec_sub2 <- prec$sub2
    prec_eps <- prec$eps

    ## clean up grp variable
    if (is.factor(grp)) {
        grp <- droplevels(grp)
    } else {
        grp <- factor(grp, levels = unique(grp))
    }

    size <- length(prec_pop)
    n_terms <- NCOL(Bmat)
    n_samples <- NROW(Bmat)
    n_subs <- length(unique(grp))
    dim_sub1 <- NROW(prec_sub1[, , 1])
    idx <- tapply(seq_len(n_samples), grp, function(x) x, simplify = FALSE)

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
    ## resampling_idx <- sample.int(length(prec_pop), size, TRUE)
    resampling_idx <- 1:size

    for (k in 1:size) {

        ## get precisions
        r <- resampling_idx[k]
        kprec_pop <- prec_pop[r]
        kprec_eps <- prec_eps[r]
        kprec_sub1 <- prec_sub1[, , r]
        kprec_sub2 <- prec_sub2[r]
        kprec_sub <- block_diag(kprec_sub1, diag(kprec_sub2, n_terms - dim_sub1))

        ## kprec_sub <- diag(prec_sub2[r], n_terms)
        ## kprec_sub[1:dim_sub1, 1:dim_sub1] <- prec_sub1[, , r]

        ## get theta
        for (i in levels(grp)) {
            ## for numerical stability, these steps are simplified
            xBmat_i <- xBmat_sub[, , i]
            Li <- xBmat_i + kprec_sub / kprec_eps
            inv_Li <- chol2inv(chol(Li))
            BMB[, , i] <- kprec_eps * (diag(n_terms) - xBmat_i %*% inv_Li) %*% xBmat_i
            ## BMB[, , i] <- kprec_eps * xBmat_i - xBmat_i %*% inv_Li %*% xBmat_i

            Bxy_i <- Bxy_sub[, i]
            BMy[, i] <- kprec_eps * (Bxy_i - xBmat_i %*% inv_Li %*% Bxy_i)
            ## BMy[, i] <- kprec_eps * Bxy_i - xBmat_i %*% inv_Li %*% Bxy_i
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



## This is a MC sampler that takes precision samples and gives constrained coefs
## A * theta >= lower ; A * (theta + delta_i) >= lower
## The size of coefs equals to the size of prec
## Requirements: Bmat, y, grp, Kmat, Amat, prec
## Algorithms paremeters: burn
## Extras: verbose
mc_cons_sub <- function(y, grp, Bmat, Kmat, Amat, lower, prec, size,
                        verbose = TRUE) {
    prec_pop <- prec$pop
    prec_sub1 <- prec$sub1
    prec_sub2 <- prec$sub2
    prec_eps <- prec$eps

    ## clean up grp variable
    if (is.factor(grp)) {
        grp <- droplevels(grp)
    } else {
        grp <- factor(grp, levels = unique(grp))
    }

    size <- length(prec_pop)
    n_terms <- NCOL(Bmat)
    n_samples <- NROW(Bmat)
    n_subs <- length(unique(grp))
    dim_sub1 <- NROW(prec_sub1[, , 1])
    n_knots <- NCOL(Bmat) -  NROW(prec_sub1[, , 1]) # n_terms - dim_sub1
    idx <- tapply(seq_len(n_samples), grp, function(x) x, simplify = FALSE)

    ## construct full constraint matrix
    fAmat_left <- t(matrix(t(Amat), NCOL(Amat), NROW(Amat) * (n_subs + 1)))
    fAmat_topright <- matrix(0, NROW(Amat), n_terms * n_subs)
    fAmat_btmright <- block_diag(Amat, size = n_subs)
    fAmat <- cbind(fAmat_left, rbind(fAmat_topright, fAmat_btmright))
    flower <- rep(lower, n_subs + 1)

    ## generate an initial value for the hmc sampler
    Ainv <- diag(NCOL(Amat))
    Ainv[row(Ainv) > diff(dim(Amat))] <- Amat
    fAinv_left <- t(matrix(t(Ainv), NCOL(Ainv), NROW(Ainv) * (n_subs + 1)))
    fAinv_topright <- matrix(0, NROW(Ainv), n_terms * n_subs)
    fAinv_btmright <- block_diag(Ainv, size = n_subs)
    fAinv <- cbind(fAinv_left, rbind(fAinv_topright, fAinv_btmright))
    fAinv <- solve(fAinv)
    kcoef <- fAinv %*% rep(c(rep(3, diff(dim(Amat))), lower + 1), times = n_subs + 1)

    ## construct full design matrix
    fBmat_right <- do.call(block_diag, lapply(idx, function(x) Bmat[x, ]))
    fBmat <- cbind(Bmat, fBmat_right)

    ## some crossproducts precalculation
    xfBmat <- crossprod(fBmat)
    fBxy <- crossprod(fBmat, y)

    ## initialise the output list
    samples <- list(population = matrix(NA, n_terms, size),
                    subjects = array(NA, c(n_terms, n_subs, size),
                                     dimnames = list(NULL, levels(grp))),
                    precision = list(pop = rep(NA, size),
                                     sub1 = array(NA, c(dim_sub1, dim_sub1, size)),
                                     sub2 = rep(NA, size),
                                     eps = rep(NA, size)))

    ## get indices of the resamped variances
    ## resampling_idx <- sample.int(length(prec_pop), size, TRUE)
    resampling_idx <- 1:size

    for (k in 1:size) {

        ## get precisions
        r <- resampling_idx[k]
        kprec_pop <- prec_pop[r]
        kprec_eps <- prec_eps[r]
        kprec_sub <- block_diag(prec_sub1[, , r], diag(prec_sub2[r], n_knots))
        kprec <- block_diag(kprec_sub, size = n_subs + 1)
        kprec[1:n_terms, 1:n_terms] <- diag(kprec_pop, n_terms)

        ## get theta and delta
        M <- chol2inv(chol(xfBmat + kprec / kprec_eps))
        mu <- M %*% fBxy
        sig <- 1 / kprec_eps * M
        kcoef <- t(tnorm::rmvtnorm(1, mu, sig, kcoef, fAmat, -flower))
        kcoef_mat <- matrix(kcoef, n_terms)

        ## print progress
        if (verbose && (k %% 1000 == 0)) {
            cat(k, " samples generated.\n")
        }

        ## store samples
        samples$precision$pop[k] <- kprec_pop
        samples$precision$sub1[, , k] <- prec_sub1[, , r]
        samples$precision$sub2[k] <- prec_sub2[r]
        samples$precision$eps[k] <- kprec_eps
        samples$population[, k] <- kcoef_mat[, 1]
        samples$subjects[, , k] <- kcoef_mat[, -1]
    }
    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    list(means = means, samples = samples)
}


