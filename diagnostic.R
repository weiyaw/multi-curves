########################################################################
##########################                    ##########################
##########################    MANIPULATION    ##########################
##########################                    ##########################
########################################################################

## convert an array of precision matrices to covariance matrices
prec_to_cov <- function(prec) {
    size <- dim(prec)[3]
    for (i in 1:size) {
        prec[, , i] <- chol2inv(chol(prec[, , i]))
    }
    prec
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
    if (!is.null(fm$samples$lp)) {
        flat[, grep("lp__", para)] <- fm$samples$lp
    }
    if (!is.null(fm$samples$ll)) {
        flat[, grep("ll__", para)] <- fm$samples$ll
    }
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
                        dimnames = list(NULL, paste("Chain", 1:n_chains), para))
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
                                     dimnames = dimnames(fms[[1]]$samples$subjects)),
                    precision = list(pop = rep(NA, size),
                                     sub1 = array(NA, c(dim_sub1, dim_sub1, size)),
                                     sub2 = rep(NA, size),
                                     eps = rep(NA, size)),
                    lp = rep(NA, size),
                    ll = rep(NA, size))
    for (i in 1:n_chains) {
        idx <- seq(start_idx[i], end_idx[i])
        samples$population[, idx] <- fms[[i]]$samples$population
        samples$subjects[, , idx] <- fms[[i]]$samples$subjects
        samples$precision$pop[idx] <- fms[[i]]$samples$precision$pop
        samples$precision$sub1[, , idx] <- fms[[i]]$samples$precision$sub1
        samples$precision$sub2[idx] <- fms[[i]]$samples$precision$sub2
        samples$precision$eps[idx] <- fms[[i]]$samples$precision$eps
        if (!is.null(fms[[i]]$samples$lp)) {
            samples$lp[idx] <- fms[[i]]$samples$lp
        }
        if (!is.null(fms[[i]]$samples$ll)) {
            samples$ll[idx] <- fms[[i]]$samples$ll
        }
    }
    means <- list(population = rowMeans(samples$population),
                  subjects = rowMeans(samples$subjects, dims = 2))
    list(means = means, samples = samples)
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


######################################################################
##########################                  ##########################
##########################    DIAGNOSTIC    ##########################
##########################                  ##########################
######################################################################

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
    c(theta_names, delta_names, cov_names, sig2_names, "lp__")
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
    stat_sub <- matrix(NA, n_delta, n_subs, dimnames = dimnames(subjects))
    for (i in colnames(stat_sub)) {
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

## return a vector of statistics calculated from "fun"
##
## this function takes n_samples * n_chain * n_parameters, combine all the
## chains, calculate the statistics, and return a vector of the statistics for
## each of the parameters.
## eg. sweep_posterior_flats(flats, sd), sweep_posterior(flats, mean)
sweep_posterior_flats <- function(flats, fun) {
    stat_all <- rep(NA, dim(flats)[3])
    names(stat_all) <- dimnames(flats)[[3]]
    for (i in names(stat_all)) {
        stat_all[i] <- fun(c(flats[, , i]))
    }
    stat_all
}

fft_next_good_size <- function(N) {
    ## Find the optimal next size for the FFT so that
    ## a minimum number of zeros are padded.
    if (N <= 2)
        return(2)
    while (TRUE) {
        m = N
        while ((m %% 2) == 0) m = m / 2
        while ((m %% 3) == 0) m = m / 3
        while ((m %% 5) == 0) m = m / 5
        if (m <= 1)
            return(N)
        N = N + 1
    }
}

#' Autocovariance estimates
#'
#' Compute autocovariance estimates for every lag for the specified
#' input sequence using a fast Fourier transform approach. Estimate
#' for lag t is scaled by N-t.
#'
#' @param y A numeric vector forming a sequence of values.
#'
#' @return A numeric vector of autocovariances at every lag (scaled by N-lag).
autocovariance <- function(y) {
    N <- length(y)
    M <- fft_next_good_size(N)
    Mt2 <- 2 * M
    yc <- y - mean(y)
    yc <- c(yc, rep.int(0, Mt2 - N))
    transform <- fft(yc)
    ac <- fft(Conj(transform) * transform, inverse = TRUE)
    ## use "biased" estimate as recommended by Geyer (1992)
    ac <- Re(ac)[1:N] / (N * 2 * N-1)
    ac
}

#' Autocorrelation estimates
#'
#' Compute autocorrelation estimates for every lag for the specified
#' input sequence using a fast Fourier transform approach. Estimate
#' for lag t is scaled by N-t.
#'
#' @param y A numeric vector forming a sequence of values.
#'
#' @return A numeric vector of autocorrelations at every lag (scaled by N-lag).
autocorrelation <- function(y) {
    ac <- autocovariance(y)
    ac <- ac / ac[1]
}

#' Rank normalization
#'
#' Compute rank normalization for a numeric array. First replace each
#' value by its rank. Average rank for ties are used to conserve the
#' number of unique values of discrete quantities. Second, normalize
#' ranks via the inverse normal transformation.
#'
#' @param x A numeric array of values.
#'
#' @return A numeric array of rank normalized values with the same
#'     size as input.
z_scale <- function(x) {
    S <- length(x)
    r <- rank(x, ties.method = 'average')
    z <- qnorm((r - 1 / 2) / S)
    if (!is.null(dim(x))) {
        ## output should have the input dimension
        z <- array(z, dim = dim(x), dimnames = dimnames(x))
    }
    z
}

#' Rank uniformization
#'
#' Compute rank uniformization for a numeric array. First replace each
#' value by its rank. Average rank for ties are used to conserve the
#' number of unique values of discrete quantities. Second, uniformize
#' ranks to scale [1/(2S), 1-1/(2S)], where S is the the number of values.
#'
#' @param x A numeric array of values.
#'
#' @return A numeric array of rank uniformized values with the same
#'     size as input.
u_scale <- function(x) {
    S <- length(x)
    r <- rank(x, ties.method = 'average')
    u <- (r - 1 / 2) / S
    if (!is.null(dim(x))) {
        ## output should have the input dimension
        u <- array(u, dim = dim(x), dimnames = dimnames(x))
    }
    u
}

#' Rank values
#'
#' Compute ranks for a numeric array. First replace each
#' value by its rank. Average rank for ties are used to conserve the
#' number of unique values of discrete quantities. Second, normalize
#' ranks via the inverse normal transformation.
#'
#' @param x A numeric array of values.
#'
#' @return A numeric array of ranked values with the same
#'     size as input.
r_scale <- function(x) {
    r <- rank(x, ties.method = 'average')
    if (!is.null(dim(x))) {
        ## output should have the input dimension
        r <- array(r, dim = dim(x), dimnames = dimnames(x))
    }
    r
}

split_chains <- function(sims) {
    ## split Markov chains
    ## Args:
    ##   sims: a 2D array of samples (# iter * # chains)
    if (is.vector(sims)) {
        dim(sims) <- c(length(sims), 1)
    }
    niter <- dim(sims)[1]
    half <- niter / 2
    cbind(sims[1:floor(half), ], sims[ceiling(half + 1):niter, ])
}

is_constant <- function(x, tol = .Machine$double.eps) {
    abs(max(x) - min(x)) < tol
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

#' Effective sample size
#'
#' Compute the effective sample size estimate for a sample of several chains
#' for one parameter. For split-ESS, call this with split chains.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for the effective sample size.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
ess_rfun <- function(sims) {
    if (is.vector(sims)) {
        dim(sims) <- c(length(sims), 1)
    }
    chains <- ncol(sims)
    n_samples <- nrow(sims)

    acov <- lapply(seq_len(chains), function(i) autocovariance(sims[, i]))
    acov <- do.call(cbind, acov)
    chain_mean <- apply(sims, 2, mean)
    mean_var <- mean(acov[1, ]) * n_samples / (n_samples - 1)
    var_plus <- mean_var * (n_samples - 1) / n_samples
    if (chains > 1)
        var_plus <- var_plus + var(chain_mean)

    ## Geyer's initial positive sequence
    rho_hat_t <- rep.int(0, n_samples)
    t <- 0
    rho_hat_even <- 1
    rho_hat_t[t + 1] <- rho_hat_even
    rho_hat_odd <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
    rho_hat_t[t + 2] <- rho_hat_odd
    while (t < nrow(acov) - 5 && !is.nan(rho_hat_even + rho_hat_odd) &&
           (rho_hat_even + rho_hat_odd > 0)) {
               t <- t + 2
               rho_hat_even = 1 - (mean_var - mean(acov[t + 1, ])) / var_plus
               rho_hat_odd = 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
               if ((rho_hat_even + rho_hat_odd) >= 0) {
                   rho_hat_t[t + 1] <- rho_hat_even
                   rho_hat_t[t + 2] <- rho_hat_odd
               }
           }
    max_t <- t
    ## this is used in the improved estimate
    if (rho_hat_even>0)
      rho_hat_t[max_t + 1] <- rho_hat_even

    ## Geyer's initial monotone sequence
    t <- 0
    while (t <= max_t - 4) {
        t <- t + 2
        if (rho_hat_t[t + 1] + rho_hat_t[t + 2] >
            rho_hat_t[t - 1] + rho_hat_t[t]) {
            rho_hat_t[t + 1] = (rho_hat_t[t - 1] + rho_hat_t[t]) / 2;
            rho_hat_t[t + 2] = rho_hat_t[t + 1];
        }
    }
    ess <- chains * n_samples
    ## Geyer's truncated estimate
    ## tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t])
    ## Improved estimate reduces variance in antithetic case
    tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t]) + rho_hat_t[max_t+1]
    ## Safety check for negative values and with max ess equal to ess*log10(ess)
    tau_hat <- max(tau_hat, 1/log10(ess))
    ess <- ess / tau_hat
    ess
}

#' Rhat convergence diagnostic
#'
#' Compute Rhat convergence diagnostic as the maximum of rank normalized
#' split-Rhat and rank normalized folded-split-Rhat for one parameter.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for the effective sample size.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
rhat <- function(sims) {
  bulk_rhat <- rhat_rfun(z_scale(split_chains(sims)))
  sims_folded <- abs(sims - median(sims))
  tail_rhat <- rhat_rfun(z_scale(split_chains(sims_folded)))
  max(bulk_rhat, tail_rhat)
}

## compute split-Rhat for every parameters calculates rank normalized
## split-Rhat and rank normalized folded-split-Rhat
## takes in n_samples * n_chain * n_parameters
rhat_flats <- function(flats) {
    hats <- rep(NA, dim(flats)[3])
    names(hats) <- dimnames(flats)[[3]]

    for (i in names(hats)) {
        hats[i] <- rhat(flats[, , i])
    }
    hats
}

## takes fms, flatten them and run rhat_flats
rhat_fms <- function(...) {
    flats <- flatten_chains(...)
    rhat_flats(flats)
}

## compute standard ESS for every parameters
## takes in n_samples * n_chain * n_parameters
ess_flats <- function(flats) {
    ess <- rep(NA, dim(flats)[3])
    names(ess) <- dimnames(flats)[[3]]

    for (i in names(ess)) {
        ess[i] <- ess_rfun(flats[, , i])
    }
    ess
}

## takes fms, flatten them and run rss_flats
ess_fms <- function(...) {
    flats <- flatten_chains(...)
    ess_flats(flats)
}

## return a summary statistics for posterior for flats matrix
summary_matrix_flats <- function(flats) {

    ## need working
    ## need to write sweep_posterior_flat

    quan025 <- function(x) quantile(x, 0.025, names = FALSE)
    quan500 <- function(x) quantile(x, 0.5, names = FALSE)
    quan975 <- function(x) quantile(x, 0.975, names = FALSE)

    tibble::tibble(Parameter = dimnames(flats)[[3]],
                   Rhat = rhat_flats(flats),
                   n_eff = ess_flats(flats),
                   mean = sweep_posterior_flats(flats, mean),
                   sd = sweep_posterior_flats(flats, sd),
                   "2.5%" = sweep_posterior_flats(flats, quan025),
                   "50%" = sweep_posterior_flats(flats, quan500),
                   "97.5%" = sweep_posterior_flats(flats, quan975))
}

## return a summary statistics for posterior
summary_matrix <- function(...) {
    flats <- flatten_chains(...)
    summary_matrix_flats(flats)
}

## rank plot
mcmc_hist_r_scale <- function(x, nbreaks = 50, ...) {
    max <- prod(dim(x)[1:2])
    bayesplot::mcmc_hist(r_scale(x),
                         breaks = seq(0, max, by = max / nbreaks) + 0.5,
                         ...) +
        theme(axis.line.y = element_blank())
}

