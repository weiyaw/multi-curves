
MultiVarCov <- function(restruct, resid.var = NA) {

    ## Return a list of variance-covariance matrices from a "reStruct" object.

    if (!class(restruct) == "reStruct") {
        stop("reStruct object not supplied")
    }
    if (is.na(resid.var)) {
        stop("Variance for scaling not provided")
    }
    lapply(lapply(restruct, as.matrix), function(x) {x * resid.var})
}

DiagMat <- function(A, size = NA) {

    ## Construct a big matrix with its diagonal elements the matrices provided
    ## in "A". If "A" is a matrix, the function returns a matrix with "size" of
    ## "A" in the diagonal elements. If "A" is a list, it returns a matrix with
    ## diagonal elements as all the matrices contained inside the list.


    if (is.matrix(A)) {
        if (is.na(size)) {
            stop("Size of the resulting matrix not supplied")
        }
        n.row <- NROW(A)
        n.col <- NCOL(A)
        row.idx <- seq(1, n.row)
        col.idx <- seq(1, n.col)
        res <- matrix(0, size * n.row, size * n.col)
        for (i in seq(0, size - 1)) {
            res[i * n.row + row.idx, i * n.col + col.idx] <- A
        }
        return(res)

   } else if (is.list(A)) {
        if (!all(sapply(A, is.matrix))) {
            stop("The list contains non-matrix objects.")
        }
        dims <- sapply(A, dim)
        total.dims <- rowSums(dims)
        res <- matrix(0, total.dims[1], total.dims[2])
        row.rolling <- col.rolling <- 0
        for (i in seq(1, NCOL(dims))) {
            row.idx <- row.rolling + seq(1, dims[1, i])
            col.idx <- col.rolling + seq(1, dims[2, i])
            res[row.idx, col.idx] <- A[[i]]
            row.rolling <- row.rolling + dims[1, i]
            col.rolling <- col.rolling + dims[2, i]
        }
        return(res)

   } else {
       warning("Non-matrix or list object supplied")
   }
}

VecToMat <- function(vec, n, colwise = TRUE) {

    ## Create a matrix with identical columns (or rows)

    if (colwise) {
        matrix(rep(vec, times = n), ncol = n, nrow = length(vec))
    } else {
        matrix(rep(vec, times = n), nrow = n, ncol = length(vec), byrow = TRUE)
    }
}

NormKernel <- function(x, mu, inv.sigma, log = TRUE) {
    ## Kernel of a normal distribution

    if (is.matrix(inv.sigma)) {
        res <- -0.5 * (t(x - mu) %*% inv.sigma %*% (x - mu))
    } else if (is.numeric(inv.sigma) && length(inv.sigma) == 1) {
        res <- -0.5 * (crossprod(x - mu) * inv.sigma)
    } else {
        stop("Invalid inverse var-cov matrix")
    }

    if (log) {
        res
    } else {
        exp(res)
    }
}

ProposalLogFac <- function(mu, G, n.subject, n.terms) {
    ## Logarithm of independent normal proposal ratio
    ## G is the varcov matrix supplied in rmvtgauss.lin
    ## mu is the EBLUPS estimate
    ## n.subject is the number of subjects

    if (NROW(G) != NCOL(G)) {
        stop("G is not a square matrix!")
    }
    if (n.terms != (NROW(G) / (n.subject + 1))) {
        stop("Number of terms disagree with dimension of G")
    }
    terms.idx <- seq(1, n.terms)

    ## Calculate list of precision matrix and mean vector
    idx.list <- list()
    P.list <- list()
    mu.list <- list()
    for (i in seq(1, n.subject + 1)) {
        idx <- (i - 1) * n.terms + terms.idx
        idx.list[[i]] <- idx
        P.list[[i]] <- solve(G[idx, idx, drop = FALSE])
        mu.list[[i]] <- mu[idx]
    }
    rm(list = setdiff(ls(), c("idx.list", "P.list", "mu.list")))

    function(zeta) {
        res <- 0
        for (i in seq_along(mu.list)) {
            idx <- idx.list[[i]]
            mu <- mu.list[[i]]
            P <- P.list[[i]]
            res <- res + NormKernel(zeta[idx], mu, P, TRUE)
        }
        return(res)
    }
}

## cand <- rmvtgauss.lin(5, crude.zeta, crude.varcov, Amat = t(A), Avec = rep(0, NROW(A)))
## prop <- ProposalLogFac(crude.zeta, crude.varcov, n.subject, n.terms)
## mu <- environment(prop)$mu.list
## P <- environment(prop)$P.list
## prop(cand[,1]) - prop(cand[,2])
## dmvnorm(cand[,1], crude.zeta, crude.varcov, TRUE) -
##     dmvnorm(cand[,2], crude.zeta, crude.varcov, TRUE)

LikelihoodLogFac <- function(y, model.mat, n.terms, n.per.sub) {
    ## model.mat is the model matrix with linear terms preceding splines terms.

    n <- sum(n.per.sub)
    n.subject <- length(n.per.sub)

    if (!(is.vector(n.per.sub) || is.array(n.per.sub))) {
        stop("Incorrect format of n.per.sub.")
    }
    if (NROW(model.mat) != n) {
        stop("Rows of model matrix mismatch with n.per.sub.")
    }
    if (NCOL(model.mat) != n.terms) {
        stop("Columns of model matrix mismatch with n.terms.")
    }

    idx.end <- cumsum(n.per.sub)
    idx <- list()

    for (i in seq_along(n.per.sub)) {
        idx[[i]] <- seq(to = idx.end[i], len = n.per.sub[i])
    }

    rm(list = setdiff(ls(), c("y", "model.mat", "n.terms", "n.subject", "idx")))

    function(zeta, eps.prec) {
        reg.coef <- VecToMat(zeta[1:n.terms], n.subject) +
            matrix(zeta[-(1:n.terms)], n.terms, n.subject)
        y.predict <- rep(NA, n)

        for (i in seq(1, n.subject)) {
            y.predict[idx[[i]]] <- model.mat[idx[[i]], ] %*% reg.coef[, i]
        }
        ## return(y.predict)
        return(NormKernel(y, y.predict, eps.prec))
    }
}

## Seems correct?
## likelihood <- with(sitka, LikelihoodLogFac(y, cbind(X, Z), n.terms, n.per.sub))
## likelihood(cand[, 1], 1 / fit1$sigma^2)
## likelihood(cand[, 2], 1 / fit1$sigma^2)
## NormKernel(y, predict(fit2), 1 / fit1$sigma^2)

## library(ggplot2)
## ggplot() +
##     geom_line(aes(scal.time, y, group = id.num), sitka) +
##     geom_line(aes(scal.time, predict(fit1), group = id.num), sitka, col = 'blue') +
##     geom_line(aes(scal.time, likelihood(cand[, 1], 1 / fit1$sigma^2), group = id.num), sitka, col = 'red') +
##     geom_line(aes(scal.time, likelihood(cand[, 2], 1 / fit1$sigma^2), group = id.num), sitka, col = 'green')


PriorLogFac <- function(n.terms, n.subject) {

    zeta.dim <- n.terms * (n.subject + 1)
    idx <- matrix(seq(1, zeta.dim), n.terms, n.subject + 1)
    beta.idx <- idx[1:2, 1]
    u.idx <- idx[-(1:2), 1]
    b.idx <- idx[1:2, -1, drop = FALSE]
    v.idx <- c(idx[-(1:2), -1])

    rm(list = setdiff(ls(), c("beta.idx", "u.idx", "b.idx", "v.idx",
                              "zeta.dim")))

    function(zeta, u.prec, v.prec, b.prec) {
        if (length(zeta) != zeta.dim) {
            stop("Length of zeta mismatch with n.terms and n.subject. ")
        }

        ## beta (standard normal prior)
        res <- NormKernel(zeta[beta.idx], 0, 1)

        ## u (truncated normal)
        res <- res + NormKernel(zeta[u.idx], 0, u.prec)

        ## b (truncated normal)
        for (i in seq(1, n.subject)) {
            res <- res + NormKernel(zeta[b.idx[, i]], 0, b.prec)
        }

        ## v (truncated normal)
        res <- res + NormKernel(zeta[v.idx], 0, v.prec)

        return(res)

    }
}

#' Construct a difference matrix
#'
#' \code{get_diff_mat} returns a difference matrix of an arbitrary order and size.
#'
#' This function returns a \eqn{k^{th}} order difference matrix \eqn{D}, such
#' that \eqn{Dx} gives a vector of \eqn{k^{th}} order differences of \eqn{x}.
#' The parameter \code{size} is usually the dimension of \eqn{x}.
#'
#' @param size the column size of the difference matrix. Required.
#'
#' @param k the order of difference. Required.
#'
#' @return A \code{size - k} by \code{size} matrix.

get_diff_mat <- function(size, k) {
    if (size <= k) {
        stop("Order of difference greater than column size.")
    }

    D <- diag(1, size)
    D[(col(D) + 1) == row(D)] <- -1
    D <- D[-1, ]
    if (k == 1) {
        return(D)
    } else {
        return(get_diff_mat(size - 1, k - 1) %*% D)
    }
}


#' Construct an equidistant B-splines design matrix
#'
#' \code{get_design_bs} returns an equidistant B-splines design matrix without
#' explicitly specifying the location of the knots.
#'
#' The knots locations are the minimum and maximum of \code{x}, and \code{K}
#' equidistant points between these two extrema. The B-splines are equidistant,
#' i.e. no multiple knots at both ends of \eqn{x}. This function uses
#' \code{splines::splineDesign}.
#'
#' @param x predictor vector. Required.
#' @param K number of inner knots, or a vector of all knots (interior,
#'     boundaries, and knots outside extrema). Required.
#' @param deg the degree of polynomial which the B-splines span. Required.
#' @param EPS tolerance error.
#'
#' @return a list with components `design' (\code{length(x)} by \code{K + deg +
#'     1} design matrix) and all `knots' (interior and boundaries knots, and
#'     knots outside extrema).
get_design_bs <- function(x, K, deg, EPS = 1e-6) {
    res <- list()

    ## get the knots
    if (length(K) == 1 && is.numeric(K)) {
        dist <- diff(range(x)) / (K + 1)
        knots <- seq(min(x) - (deg * dist) - EPS, max(x) + (deg * dist) + EPS,
                     len = K + 2 * (deg + 1))
        names(knots) <- NULL
        res$knots <- knots
    } else if (is.vector(K) && is.numeric(K)) {
        knots <- K
        res$knots <- knots
    } else {
        stop("Supplied knots must a numeric vector.")
    }
    names(knots) <- NULL

    res$design <- splines::splineDesign(knots, x, ord = deg + 1)
    res$knots <- knots
    return(res)
}


## Return the constraint matrix A, where Ax > 0.
## size : column size of the matrix
## shape : increasing, decreasing
get_constmat_bs <- function(size, shape) {
    D <- get_diff_mat(size, 1)
    if (shape == "increasing") {
        return(D)
    } else if (shape == "decreasing") {
        warning("decresing cases untested.")
        return(-1 * D)
    } else {
        stop("Unrecognised shape.")
    }
}

## Return a design matrix at quartile (or supplied) knots, and all knots.
## The knots are taken at the quartiles, and min and max of x.
## x: all the explanatory data
## K: number of inner knots, or a vector of all knots (including extrema)
## deg: degree of the spline polynomial
get_design_tpf <- function(x, K, deg) {

    res <- list()

    ## get the knots
    if (length(K) == 1 && is.numeric(K)) {
        knots <- quantile(unique(x), seq(0, 1, len = K + 2))[-c(1, K + 2)]
        names(knots) <- NULL
        res$knots <- c(min(x), knots, max(x))
    } else if (is.vector(K) && is.numeric(K)) {
        knots <- K[-c(1, length(K))]
        res$knots <- K
    } else {
        stop("Supplied knots must a numeric vector.")
    }
    names(knots) <- NULL

    ## get the design matrix
    splines <- outer(x, knots, `-`)
    ## if (deg == 1) {
    ##     splines <- splines * (splines > 0)
    ##     res$design <- cbind(1, x, splines, deparse.level = 0)
    ## } else if (deg == 2) {
    ##     splines <- splines^2 * (splines > 0)
    ##     res$design <- cbind(1, x, x^2, splines, deparse.level = 0)
    ## } else {
    ##     stop("Invalid degree.")
    ## }
    ## colnames(res$design) <- NULL

    ## this chunck can be generalised to higher degree polynomials
    if (deg > 0 && deg < 4) {
        splines <- splines^deg * (splines > 0)
        res$design <- cbind(1, poly(x, degree = deg, raw = TRUE, simple = TRUE),
                            splines, deparse.level = 0)
        colnames(res$design) <- NULL
    } else {
        stop("Invalid degree.")
    }

    res
}

## Return the constraint matrix A, where Ax > 0.
## knots: all knots, i.e. extrama of x and inner knots. If deg = 1, this can
##        the number of inner knots
## shape: increasing, decreasing
## deg: degree of the polynomial splines
get_constmat_tpf <- function(knots, shape, deg) {
    if (deg == 1) {
        if (length(knots) == 1) {
            knots <- rep(NA, knots + 2) # a dummy knots vector
        }
        if (shape == "increasing") {
            A <- diag(length(knots) - 1)
            A[lower.tri(A)] <- 1
            A <- cbind(0, A)
            return(A)
        } else if (shape == "decreasing") {
            warning("decresing cases untested.")
            A <- diag(length(knots) - 1)
            A[lower.tri(A)] <- 1
            A <- cbind(0, A)
            return(-1 * A)
        } else {
            stop("Unrecognised shape.")
        }
    } else if (deg == 2) {
        if (!is.vector(knots, "numeric") || length(knots) == 1) {
            stop("Knots vector must be supplied for quadratic spline.")
        } else if (shape == "increasing") {
            A <- 2 * outer(knots[-(1:2)], knots[-c(1, length(knots))], `-`)
            A[upper.tri(A)] <- 0
            A <- rbind(0, 0, A, deparse.level = 0)
            A <- cbind(0, 1, 2 * knots, A, deparse.level = 0)
            return(A)
        } else if (shape == "decreasing") {
            warning("decresing cases untested.")
            A <- 2 * outer(knots[-(1:2)], knots[-c(1, length(knots))], `-`)
            A[upper.tri(A)] <- 0
            A <- rbind(0, 0, A, deparse.level = 0)
            A <- cbind(0, 1, 2 * knots, A, deparse.level = 0)
            return(-1 * A)
        } else {
            stop("Unrecognised shape.")
        }
    } else {
        stop("Unrecognised spline degree.")
    }
}

## delete the samples with indices 1:size
truncate_fm <- function(model, size) {
    trc_idx <- -1 * seq_len(size)
    helper_rm <- function(x) {
        if (is.array(x)) {
            if (length(dim(x)) == 3) {
                x[, , trc_idx]
            } else if (length(dim(x)) == 2) {
                x[, trc_idx]
            }
        } else if (is.vector(x)) {
            x[trc_idx]
        }  else {
            NULL
        }
    }
    helper_mean <- function(x) {
        if (is.array(x)) {
            if (length(dim(x)) == 3) {
                rowMeans(x, dims = 2)
            } else if (length(dim(x)) == 2) {
                rowMeans(x)
            }
        } else if (is.vector(x)) {
            mean(x)
        }  else {
            NULL
        }
    }
    model$samples <- rapply(model$samples, helper_rm, how = "list")
    model$means <- rapply(model$samples, helper_mean, how = "list")
    model
}


## Get an ordinary least squares estimate with a given response and design
## matrix.
get_ols <- function(response, design) {
    tcrossprod(solve(crossprod(design)), design) %*% as.vector(response)
}


## Plot multiple lines
## Requires: model$basis$knots (including extrema), model$basis$type,
##   model$basis$degree, model$info$lvl_pop, model$data, model$mean
##
plot_spline <- function(model, limits = NULL, plot_which = NULL, fine = 200, mle = FALSE) {

    EPS <- 1e-6
    knots <- model$basis$knots
    names(knots) <- NULL
    type <- model$basis$type
    deg <- model$basis$degree

    ## Rename the input data to prevent future accidental name change.
    data <- model$data
    colnames(data) <- c("x", "y", "sub", "pop")

    is_multi_pop <- !is.null(model$info$lvl_pop)

    ## Use the range of the predictor if 'limits' is missing
    if (is.null(limits)) {
        limits <- range(data$x)
    } else if (!(is.vector(limits) && length(limits) == 2)) {
        stop("'limits' must be a vector of length 2.")
    }

    ## Plot all population if 'plot_which' is missing
    if (is.null(plot_which)) {
        plot_which <- model$info$lvl_pop
    } else if (!is.vector(lvl_pop)) {
        stop("'plot_which' must be a vector.")
    } else if (!all(plot_which %in% lvl_pop)) {
        stop("Invalid population levels in 'plot_which'.")
    }

    ## 'x' axis of the plot
    knots_within <- knots[knots > (min(limits) - EPS) &
                          knots < (max(limits) + EPS)]
    plot_x <- unique(c(seq(min(limits), max(limits), length.out = fine),
                       knots_within))
    plot_x <- plot_x[order(plot_x)]

    ## Model matrix
    if (type == "tpf") {
        model_mat <- get_design_tpf(plot_x, knots, deg)$design
    } else if (type == "bs") {
        model_mat <- splines::splineDesign(knots, plot_x, ord = deg + 1,
                                           outer.ok = TRUE)
    } else if (type == "bs_hier") {
        model_mat <- splines::splineDesign(knots, plot_x, ord = deg + 1,
                                           outer.ok = TRUE)
        model_mat <- model_mat %*% model$basis$trans_mat
    } else {
        stop("Unknown type of model.")
    }

    if (mle) {
        ## Extract the mle (posterior mode)
        coef_pop <- model$mle$population
        dev_sub <- model$mle$subjects
    } else {
        ## Extract the posterior means
        coef_pop <- model$means$population
        dev_sub <- model$means$subjects
    }

    ## Helper function to calculate the y axis of the plot
    get_plot_y <- function(coef) {
        model_mat %*% coef
    }

    ## y axis of the plot
    if (is_multi_pop) {
        coef_sub <- mapply(`+`, dev_sub[plot_which], coef_pop[plot_which],
                           SIMPLIFY = FALSE)
        plot_y_pop <- lapply(coef_pop[plot_which], get_plot_y)
        plot_y_sub <- lapply(coef_sub[plot_which], get_plot_y)
    } else {
        coef_sub <- dev_sub + coef_pop
        ## Align the data format of a single population model to a multiple
        ## population model. Create a dummy name for the population
        plot_which <- "dummy"
        data$pop <- "dummy"
        plot_y_pop <- list(dummy = get_plot_y(coef_pop))
        plot_y_sub <- list(dummy = get_plot_y(coef_sub))
    }

    ## Reformat dataframe for ggplot
    plotdat_pop <- list()
    plotdat_sub <- list()
    for (i in plot_which) {
        plotdat_pop[[i]] <- data.frame(x = plot_x, y = plot_y_pop[[i]])
        plotdat_sub[[i]] <- reshape2::melt(plot_y_sub[[i]],
                                           varnames = c("x", "sub"),
                                           as.is = TRUE, value.name = "y")
        plotdat_sub[[i]]$x <- plot_x
    }
    aes <- ggplot2::aes
    geom_point <- ggplot2::geom_point
    geom_line <- ggplot2::geom_line

    for (i in plot_which) {
        ## for each population
        ggobj <- ggplot2::ggplot(mapping = aes(x, y, col = sub)) +
            geom_point(data = data[data$pop == i, ]) +             # dataset
            geom_line(aes(group = sub), data = plotdat_sub[[i]]) + # subject-specific curves
            geom_line(aes(col = NULL), data = plotdat_pop[[i]])    # population curve
        ## + geom_line(aes(col = NULL), data = ols, col = 'red')
        print(ggobj + theme_bw() + theme(legend.position="none"))
    }

    ## This is a reassurance step (reorder group levels)
    ## plot.data$grps <- factor(plot.data$grps, levels = colnames(coefs))

    ## plot.data$x <- rep(x, times = NCOL(coefs))
    ## pop.idx <- plot.data$grps %in% "population"


    ## colnames(data) <- c("x", "y", "grps")
    ## base <- ggplot2::ggplot()
    ## gg_pop <- ggplot2::geom_line(ggplot2::aes(plot_x, plot_y, group = grps, col = grps),
    ##                           plot.data[!pop.idx, ])

    ##         ggplot2::geom_point(ggplot2::aes(x, y, col = grps), data) +
    ##         ggplot2::geom_line(ggplot2::aes(x, value), plot.data[pop.idx, ],
    ##                            col = "black") +
    ##         ggplot2::geom_point(ggplot2::aes(x = knots.within,
    ##                                          y = rep(0, length(knots.within))),
    ##                             col = "red", shape = 4)
}

## convert an array of precision matrices to covariance matrices
prec_to_cov <- function(prec) {
    size <- dim(prec)[3]
    for (i in 1:size) {
        prec[, , i] <- chol2inv(chol(prec[, , i]))
    }
    prec
}

## return a vector of parameter names
para_names <- function(fm) {
    population <- fm$samples$population
    subjects <- fm$samples$subjects
    precision <- fm$samples$precision

    n_theta <- NROW(population)
    n_subs <- dim(subjects)[2]
    n_delta <- dim(subjects)[1]
    dim_sub1 <- NCOL(precision$sub1)

    theta_names <- rep(NA, n_theta)
    for (i in 1:n_theta) {
        theta_names[i] <- paste0("theta[", i, "]")
    }
    delta_names <- matrix(NA, n_delta, n_subs)
    for (i in 1:n_subs) {
        for (j in 1:n_delta) {
            delta_names[j, i] <- paste0("delta[", j, ",", i, "]")
        }
    }
    cov_names <- matrix(NA, dim_sub1, dim_sub1)
    for (i in 1:dim_sub1) {
        for (j in 1:dim_sub1) {
            cov_names[j, i] <- paste0("cov_delta1[", j, ",", i, "]")
        }
    }
    sig2_names <- c("sig2_theta", "sig2_delta2", "sig2_eps")
    c(theta_names, delta_names, cov_names, sig2_names)
}

## return a vector of statistics calculated from "fun"
## eg. sweep_posterior(fm, sd), sweep_posterior(fm, mean)
sweep_posterior <- function(fm, fun) {
    population <- fm$samples$population
    subjects <- fm$samples$subjects
    precision <- fm$samples$precision

    n_subs <- dim(subjects)[2]
    n_delta <- dim(subjects)[1]
    dim_sub1 <- NCOL(precision$sub1)

    stat_pop <- apply(population, 1, fun)
    stat_sub <- matrix(NA, n_delta, n_subs, dimnames = list(NULL, levels(grp)))
    for (i in levels(grp)) {
        stat_sub[, i] <- apply(subjects[, i, ], 1, fun)
    }
    stat_cov_pop <- fun(1 / precision$pop)
    stat_cov_sub1 <- matrix(NA, dim_sub1, dim_sub1)
    cov_sub1 <- prec_to_cov(precision$sub1)
    for (i in 1:dim_sub1) {
        stat_cov_sub1[, i] <- apply(cov_sub1[, i, ], 1, fun)
    }
    stat_cov_sub2 <- fun(1 / precision$sub2)
    stat_cov_eps <- fun(1 /precision$eps)
    stat_all <- c(stat_pop, stat_sub, stat_cov_sub1, stat_cov_pop,
                  stat_cov_sub2, stat_cov_eps)
    names(stat_all) <- para_names(fm)
    stat_all
}

## split Markov chains (from Aki)
## sims: a 2D array of samples (# iter * # chains)
split_chains <- function(sims) {
    if (is.vector(sims)) {
        dim(sims) <- c(length(sims), 1)
    }
    niter <- dim(sims)[1]
    half <- niter / 2
    cbind(sims[1:floor(half), ], sims[ceiling(half + 1):niter, ])
}

#' Traditional Rhat convergence diagnostic
#'
#' Compute the Rhat convergence diagnostic for a single parameter
#' For split-Rhat, call this with split chains.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for Rhat.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
rhat_rfun <- function(sims) {
    if (is.vector(sims)) {
        dim(sims) <- c(length(sims), 1)
    }
    chains <- ncol(sims)
    n_samples <- nrow(sims)
    chain_mean <- numeric(chains)
    chain_var <- numeric(chains)
    for (i in seq_len(chains)) {
        chain_mean[i] <- mean(sims[, i])
        chain_var[i] <- var(sims[, i])
    }
    var_between <- n_samples * var(chain_mean)
    var_within <- mean(chain_var)
    sqrt((var_between / var_within + n_samples - 1) / n_samples)
}


## takes fms and compute rhat for every parameters
rhats <- function(...) {
    flats <- flatten_chains(...)
    hats <- rep(NA, dim(flats)[3])
    names(hats) <- dimnames(flats)[[3]]

    for (i in names(hats)) {
        hats[i] <- rhat_rfun(flats[, , i])
    }
    hats
}

## convert fm$samples into a matrix of n_samples * n_parameters
flatten_chain <- function(fm) {
    n_terms <- NROW(fm$samples$population)
    n_subs <- dim(fm$samples$subjects)[2]
    dim_sub1 <- dim(fm$samples$precision$sub1)[2]
    size <- NCOL(fm$samples$population)
    ## order: pop, sub, cov_sub2, sig2_pop, sig2_sub2, sig2_eps
    para <- para_names(fm)
    flat <- matrix(NA, size, length(para), dimnames = list(NULL, para))

    flat[, grep("^theta", para)] <- t(fm$samples$population)
    dim(fm$samples$subject) <- c(n_terms * n_subs, size)
    flat[, grep("^delta", para)] <- t(fm$samples$subject)
    dim(fm$samples$precision$sub1) <- c(dim_sub1^2, size)
    flat[, grep("^cov_delta1", para)] <- t(fm$samples$precision$sub1)
    flat[, grep("sig2_theta", para)] <- fm$samples$precision$pop
    flat[, grep("sig2_delta2", para)] <- fm$samples$precision$sub2
    flat[, grep("sig2_eps", para)] <- fm$samples$precision$eps

    flat
}

## flatten multiple fms$samples into an array of n_samples * n_chain * n_parameters
## for Rhat calculation
flatten_chains <- function(...) {
    fms <- list(...)
    n_chains <- length(fms)
    size <- NCOL(fms[[1]]$samples$population)
    ## order: pop, sub, cov_sub2, sig2_pop, sig2_sub2, sig2_eps
    para <- para_names(fms[[1]])
    flats <- array(NA, c(size, n_chains, length(para)),
                        dimnames = list(NULL, NULL, para))
    for (i in 1:n_chains) {
        flats[, i, ] <- flatten_chain(fms[[i]])
    }
    flats
}

## combine multiple fms into one fm
combine_fm <- function(...) {
    fms <- list(...)
    n_chains <- length(fms)
    n_terms <- NROW(fms[[1]]$samples$population)
    n_subs <- dim(fms[[1]]$samples$subjects)[2]
    dim_sub1 <- dim(fms[[1]]$samples$precision$sub1)[2]
    ind_size <- rep(NA, n_chains)
    for (i in 1:n_chains) {
        ind_size[i] <- NCOL(fms[[i]]$samples$population)
    }
    size <- sum(ind_size)
    start_idx <- c(0, cumsum(ind_size)[-n_chains]) + 1
    end_idx <- cumsum(ind_size)
    samples <- list(population = matrix(NA, n_terms, size),
                    subjects = array(NA, c(n_terms, n_subs, size),
                                     dimnames = list(NULL, levels(grp))),
                    precision = list(pop = rep(NA, size),
                                     sub1 = array(NA, c(dim_sub1, dim_sub1, size)),
                                     sub2 = rep(NA, size),
                                     eps = rep(NA, size)))
    for (i in 1:n_chains) {
        idx <- seq(start_idx[i], end_idx[i])
        samples$population[, idx] <- fms[[i]]$samples$population
        samples$subjects[, , idx] <- fms[[i]]$samples$subjects
        samples$precision$pop[idx] <- fms[[i]]$samples$precision$pop
        samples$precision$sub1[, , idx] <- fms[[i]]$samples$precision$sub1
        samples$precision$sub2[idx] <- fms[[i]]$samples$precision$sub2
        samples$precision$eps[idx] <- fms[[i]]$samples$precision$eps
    }
    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    list(means = means, samples = samples)
}

## return a summary statistics for posterior
summary_matrix <- function(...) {
    combined <- combine_fm(...)
    quan025 <- function(x) quantile(x, 0.025, names = FALSE)
    quan500 <- function(x) quantile(x, 0.5, names = FALSE)
    quan975 <- function(x) quantile(x, 0.975, names = FALSE)

    tibble::tibble(Parameter = para_names(combined),
                   Rhat = rhats(...),
                   n_eff = sweep_posterior(combined, mcmcse::ess),
                   mean = sweep_posterior(combined, mean),
                   sd = sweep_posterior(combined, sd),
                   "2.5%" = sweep_posterior(combined, quan025),
                   "50%" = sweep_posterior(combined, quan500),
                   "97.5%" = sweep_posterior(combined, quan975))
}



