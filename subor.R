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

new.rmvtgauss.lin <- function (n, mu, Sigma, PrecMat = FALSE,
                               factorised = FALSE, Amat, Avec, start,
                               burnin = 1, thin = 1) {
    if (length(n) == 1) {
        n <- as.integer(n)
    } else {
        n <- length(n)
    }
    if (is.na(n) || n < 0) {
        stop("argument 'n' invalid")
    }
    if (missing(mu)) {
        stop("argument 'mu' is missing with no default")
    }
    mu <- as.numeric(mu)
    d <- length(mu)
    if (missing(Sigma)) {
        stop("argument 'Sigma' is missing with no default")
    }
    Sigma <- as.matrix(Sigma)
    if (missing(Amat)) {
        stop("argument 'Amat' is missing with no default")
    }
    Amat <- as.matrix(Amat)
    if (missing(Avec)) {
        stop("argument 'Avec' is missing with no default")
    }
    Avec <- as.numeric(Avec)
    if (missing(start)) {
        start <- rep(0, d)
    }
    start <- as.numeric(start)
    if (length(burnin) > 1) {
        warning("'burnin' has length > 1; using first component")
        burnin <- as.numeric(burnin[1])
    }
    if (length(thin) > 1) {
        warning("'thin' has length > 1; using first component")
        thin <- as.numeric(thin[1])
    }
    if (!is.numeric(Sigma) || !is.numeric(Amat) || !is.numeric(burnin) ||
        !is.numeric(thin)) {
        stop("some arguments are not numeric")
    }
    if (any(is.na(mu)) || any(is.na(Sigma)) || any(is.na(Amat)) ||
        any(is.na(Avec)) || any(is.na(start)) || is.na(burnin) ||
        is.na(thin)) {
        stop("some arguments contain NAs")
    }
    if (!isTRUE(all.equal(dim(Sigma), c(d, d)))) {
        stop("arguments 'mu' and 'Sigma' are incompatible")
    }
    if (length(start) != d) {
        stop("arguments 'mu' and 'start' are incompatible")
    }
    if (NROW(Amat) != d) {
        stop("arguments 'mu' and 'Amat' are incompatible")
    }
    if (NCOL(Amat) != length(Avec)) {
        stop("arguments 'Amat' and 'Avec' are incompatible")
    }

    if (PrecMat) {
        Sigma <- solve(Sigma)
    }

    if (!factorised) {
        Right <- chol(Sigma)
    }

    Bmat <- Right %*% Amat

    Bvec <- Avec - crossprod(Amat, mu)

    ## eigv.Sigma <- eigen(Sigma)
    ## sqrt.Sigma <- with(eigv.Sigma, tcrossprod(vectors %*% diag(sqrt(values)),
    ##                                           vectors))

    new.start <- backsolve(Right, start - mu, transpose = TRUE)
    ## new.start <- solve(sqrt.Sigma, start - mu)

    res <- matrix(NA, nrow = d, ncol = n)
    for (i in 1:burnin) {
        tmpmat <- new.start * Bmat
        for (dd in 1:d) {
            tmpvec <- Bvec - colSums(tmpmat[-dd, , drop = FALSE])
            ind <- Bmat[dd, , drop = FALSE] > 0
            if (sum(ind) > 0) {
                left <- max(tmpvec[ind]/Bmat[dd, ind])
            } else {
                left <- -Inf
            }
            ind <- Bmat[dd, , drop = FALSE] < 0
            if (sum(ind) > 0) {
                right <- min(tmpvec[ind]/Bmat[dd, ind])
            } else {
                right <- Inf
            }
            if (left > right) {
                browser()
                stop("inconsistent constraints")
            }
            new.start[dd] <- rtgauss(1, left = left, right = right)
            tmpmat[dd, ] <- new.start[dd] * Bmat[dd, , drop = FALSE]
        }
    }
    res[, 1] <- new.start
    if (n > 1) {
        for (i in 2:n) {
            for (j in 1:thin) {
                tmpmat <- new.start * Bmat
                for (dd in 1:d) {
                  tmpvec <- Bvec - colSums(tmpmat[-dd, , drop = FALSE])
                  ind <- Bmat[dd, , drop = FALSE] > 0
                  if (sum(ind) > 0) {
                    left <- max(tmpvec[ind]/Bmat[dd, ind])
                  } else {
                    left <- -Inf
                  }
                  ind <- Bmat[dd, , drop = FALSE] < 0
                  if (sum(ind) > 0) {
                    right <- min(tmpvec[ind]/Bmat[dd, ind])
                  } else {
                    right <- Inf
                  }
                  if (left > right) {
                      stop("inconsistent constraints")
                  }
                  new.start[dd] <- rtgauss(1, left = left, right = right)
                  tmpmat[dd, ] <- new.start[dd] * Bmat[dd, , drop = FALSE]
                }
            }
            res[, i] <- new.start
        }
    }

    return(crossprod(Right, res) + mu)
    ## return(crossprod(sqrt.Sigma, res) + mu)
}


#' Plot the means of linear spline models.
#'
#' \code{PlotLinMean} plots linear splines response curves that correspond to
#' each posterior means.
#'
#' Plot the linear splines response curves that correspond to the posterior
#' means of the population and subject-specific posteriors, over the region
#' given in \code{limits}. The locations of knots should be supplied. A layer of
#' scatterplot can be added on top of it.
#'
#' @param post a list of matrices whose columns are the samples from the
#'     posterior. \code{post[[i]]} should give a matrix that contains the
#'     samples from the \eqn{i^th} subject posterior. Likewise, \code{post$pop}
#'     corresponds to the population posterior. Required.
#'
#' @param knots a vector of numeric that specifies the locations of the knots
#'     used in the linear splines model. The length of \code{knots} should be
#'     two less than the \code{NROW} of the matrices in \code{post}. Required.
#'
#' @param limits a vector of length 2 that specifies the x region on which the
#'     response curves should be plotted. Required.
#'
#' @param data a data frame that gives that data for the scatterplot. First,
#'     second, and third columns should correspond to the x-coordinates,
#'     y-coordinates, and a factor that specifies the subjects to which the
#'     datapoints correspond. Required
#'
#' @return A nice graph of subject-specific curves.
PlotLinMean <- function(post, knots, limits, data) {

    n.curves <- length(post)
    n.terms <- NROW(post$pop)

    ## Initialise the regression coefficients for each response curve.
    mean.set <- matrix(NA, n.terms, n.curves)
    colnames(mean.set) <- names(post)
    mean.set[, n.curves] <- rowMeans(post$pop)
    for (i in seq_len(n.curves - 1)) {
        mean.set[, i] <- mean.set[, n.curves] + rowMeans(post[[i]])
    }
    PlotLinSpline(mean.set, knots, limits, data)
}


#' Plot a linear spline model.
#'
#' \code{PlotLinSpline} plots linear splines response curves from the given
#' coefficients.
#'
#' Plot the linear splines response curves from the given coefficients over the
#' region given in \code{limits}. The locations of knots should be supplied. A
#' layer of scatterplot can be added on top of it.
#'
#' @param coefs A matrix whose columns are the coefficients of corresponding
#'     curves. In other words, the number of columns is equal to the number of
#'     curves. The column names of this matrix will be used as the names of
#'     corresponding curves. The "population" curve will be colored in
#'     black. Required.
#'
#' @param knots A vector of numeric that specifies the locations of the knots
#'     used in the linear splines model. The length of \code{knots} should be
#'     two less than the \code{NROW} of the matrices in \code{post}. Required.
#'
#' @param bases The type of bases used for fitting the model. Possible values
#'     are "tpf": truncated power functions; "bs": b-splines.
#'
#' @param limits A vector of length 2 that specifies the x region on which the
#'     response curves should be plotted. Required.
#'
#' @param data A data frame that gives that data for the scatterplot. First,
#'     second, and third columns should correspond to the x-coordinates,
#'     y-coordinates, and a factor that specifies the subjects to which the
#'     datapoints correspond. The names of the factor should agree with that
#'     provided in the \code{coefs} matrix.
#'
#' @return A nice graph of subject-specific curves.
PlotLinSpline <- function(coefs, knots, bases = "tpf", limits, data) {

    library(ggplot2)

    if (is.vector(coefs)) {
        coefs <- as.matrix(coefs)
    } else if (!is.matrix(coefs)) {
        stop("coefs must be a matrix")
    }

    if (!is.vector(knots) || !is.vector(limits)) {
        stop("knots and limits must be vectors")
    }

    rownames(coefs) <- NULL
    names(knots) <- NULL

    ## "x" terms
    x <- c(seq(min(limits), max(limits), length.out = 100), knots)
    x <- x[order(x)]

    if (bases == "tpf" && is.vector(knots) && is.numeric(knots)) {
        ## spline basis terms for truncated power functions
        basis <- outer(x, knots, `-`)
        basis <- basis * (basis > 0)
        model.mat <- cbind(1, x, basis, deparse.level = 0)
    } else if (bases == "bs" && length(knots) == 1 && is.numeric(knots)) {
        ## spline basis terms for b-splines
        basis <- outer(x, knots, `-`)
        basis <- basis * (basis > 0)
        model.mat <- cbind(1, x, basis, deparse.level = 0)
    } else {
        stop("Incompatible given knots and bases type.")
    }

    ## Calculate the response and melt them into a data frame
    y <- model.mat %*% coefs
    plot.data <- reshape2::melt(y)
    colnames(plot.data) <- c("points", "grps", "value")

    # This step is not necessary (reorder group levels)
    plot.data$grps <- factor(plot.data$grps, levels = colnames(coefs))

    plot.data$x <- rep(x, times = NCOL(coefs))
    pop.idx <- plot.data$grps %in% "population"


    if (missing(data)) {
        ggplot2::ggplot() +
            ggplot2::geom_line(aes(x, y, group = grps, col = grps),
                               plot.data[!pop.idx, ]) +
            ggplot2::geom_line(aes(x, value), plot.data[pop.idx, ],
                               col = "black")
    } else if (is.data.frame(data)) {
        ggplot2::ggplot() +
            ggplot2::geom_line(aes(x, value, group = grps, col = grps),
                               plot.data[!pop.idx, ]) +
            ggplot2::geom_point(aes(data[[1]], data[[2]], col = data[[3]])) +
            ggplot2::geom_line(aes(x, value), plot.data[pop.idx, ],
                               col = "black")

    }
}



#' Plot a quadratic spline model.
#'
#' \code{PlotQuadSpline} plots quadratic splines response curves from the given
#' coefficients.
#'
#' Plot the quadratic splines response curves from the given coefficients over the
#' region given in \code{limits}. The locations of knots should be supplied. A
#' layer of scatterplot can be added on top of it.
#'
#' @param coefs a matrix whose columns are the coefficients of corresponding
#'     curves. In other words, the number of columns is equal to the number of
#'     curves. The column names of this matrix will be used as the names of
#'     corresponding curves. The "population" curve will be colored in
#'     black. Required.
#'
#' @param knots a vector of numeric that specifies the locations of the knots
#'     used in the linear splines model. The length of \code{knots} should be
#'     two less than the \code{NROW} of the matrices in \code{post}. Required.
#'
#' @param limits a vector of length 2 that specifies the x region on which the
#'     response curves should be plotted. Required.
#'
#' @param data a data frame that gives that data for the scatterplot. First,
#'     second, and third columns should correspond to the x-coordinates,
#'     y-coordinates, and a factor that specifies the subjects to which the
#'     datapoints correspond. The names of the factor should agree with that
#'     provided in the \code{coefs} matrix.
#'
#' @return A nice graph of subject-specific curves.
PlotQuadSpline <- function(coefs, knots, limits, data) {

    if (!(require(ggplot2) && require(reshape2))) {
        stop("Couldn't load ggplot2 and reshape2")
    }

    if (is.vector(coefs)) {
        coefs <- as.matrix(coefs)
    } else if (!is.matrix(coefs)) {
        stop("coefs must be a matrix")
    }

    if (!is.vector(knots) || !is.vector(limits)) {
        stop("knots and limits must be vectors")
    }

    rownames(coefs) <- NULL
    names(knots) <- NULL

    ## "time" terms
    time <- c(seq(min(limits), max(limits), length.out = 200), knots)
    time <- time[order(time)]

    ## spline basis terms
    basis <- outer(time, knots, `-`)
    basis <- basis * (basis > 0)
    basis <- basis^2

    ## model matrix
    model.mat <- cbind(1, time, time^2, basis, deparse.level = 0)

    ## Calculate the response and melt them into a data frame
    y <- model.mat %*% coefs
    plot.data <- melt(y)
    colnames(plot.data) <- c("points", "grps", "value")

    # This is a reassurance step (reorder group levels)
    plot.data$grps <- factor(plot.data$grps, levels = colnames(coefs))

    plot.data$time <- rep(time, times = NCOL(coefs))
    pop.idx <- plot.data$grps %in% "population"

    if (missing(data)) {
        ggplot() +
            geom_line(aes(time, y, group = grps, col = grps),
                      plot.data[!pop.idx, ]) +
            geom_line(aes(time, value), plot.data[pop.idx, ], col = "black")
    } else if (is.data.frame(data)) {
        ggplot() +
            geom_line(aes(time, value, group = grps, col = grps),
                      plot.data[!pop.idx, ]) +
            geom_point(aes(data[[1]], data[[2]], col = data[[3]])) +
            geom_line(aes(time, value), plot.data[pop.idx, ], col = "black")

    }
}
