
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
#' \code{DiffMat} returns a difference matrix of an arbitrary order and size.
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

DiffMat <- function(size, k) {
    if (size <= k) {
        stop("Order of difference greater than column size.")
    }

    D <- diag(1, size)
    D[(col(D) + 1) == row(D)] <- -1
    D <- D[-1, ]
    if (k == 1) {
        return(D)
    } else {
        return(DiffMat(size - 1, k - 1) %*% D)
    }
}


#' Construct an equidistant B-splines design matrix
#'
#' \code{BSplineDesign} returns an equidistant B-splines design matrix without
#' explicitly specifying the location of the knots.
#'
#' The knots locations are the minimum and maximum of \code{x}, and \code{K}
#' equidistant points between these two extrema. The B-splines are equidistant,
#' i.e. no multiple knots at both ends of \eqn{x}. This function uses
#' \code{splines::splineDesign}.
#'
#' @param x predictor vector. Required.
#' @param K number of inner knots. Required.
#' @param deg the degree of polynomial which the B-splines span. Required.
#' @param EPS tolerance error.
#'
#' @return a list with components `design' (\code{length(x)} by \code{K + deg + 1}
#'     design matrix) and all `knots' (inner and outer).
BSplineDesign <- function(x, K, deg, EPS = 1e-6) {
    dist <- diff(range(x)) / (K + 1)
    knots <- seq(min(x) - (deg * dist) - EPS, max(x) + (deg * dist) + EPS,
                 len = K + 2 * (deg + 1))
    res <- list()
    res$design <- splines::splineDesign(knots, x, ord = deg + 1)
    res$knots <- knots
    return(res)
}


## Return the constraint matrix A, where Ax > 0.
## size : column size of the matrix
## shape : increasing, decreasing
BSplineConstMat <- function(size, shape) {
    if (shape == "increasing") {
        A <- diag(1, size)
        A[(col(A) + 1) == row(A)] <- -1
        A <- A[-1, ]
        return(A)
    } else if (shape == "decreasing") {
        warning("decresing cases untested.")
        A <- diag(1, size)
        A[(col(A) + 1) == row(A)] <- -1
        A <- A[-1, ]
        return(-1 * A)
    } else {
        stop("Unrecognised shape.")
    }
}

## Return a design matrix at quartile (or supplied) knots, and all knots.
## x: all the explanatory data
## K: number of inner knots, or a vector of all knots (including extrema)
## deg: degree of the spline polynomial
TpfDesign <- function(x, K, deg) {

    res <- list()

    if (length(K) == 1 && is.numeric(K)) {
        knots <- quantile(unique(x), seq(0, 1, len = K + 2))[-c(1, K + 2)]
        names(knots) <- NULL
        res$knots <- c(min(x), knots, max(x))
    } else if (is.vector(K) && min(K) >= min(x) && max(K) <= max(x)) {
        knots <- K[-c(1, length(K))]
        res$knots <- K
    } else {
        stop("Invalid knots. Knots should be within the extrema.")
    }
    names(knots) <- NULL

    splines <- outer(x, knots, `-`)
    if (deg > 0 && deg < 3) {
        splines <- splines^deg * (splines > 0)
        res$design <- cbind(1, poly(x, degree = deg, raw = T, simple = TRUE),
                            splines, deparse.level = 0)
        colnames(res$design) <- NULL
    } else {
        stop("Invalid degree.")
    }

    ## if (deg == 1) {
    ##     splines <- splines * (splines > 0)
    ##     res$design <- cbind(1, x, splines, deparse.level = 0)
    ## } else if (deg == 2) {
    ##     splines <- splines^2 * (splines > 0)
    ##     res$design <- cbind(1, x, x^2, splines, deparse.level = 0)
    ## } else {
    ##     stop("Invalid degree.")
    ## }

    return(res)
}

## Return the constraint matrix A, where Ax > 0.
## knots: all knots, i.e. extrama of x and inner knots. If deg = 1, this can
##        the number of inner knots
## shape: increasing, decreasing
## deg: degree of the polynomial splines
TpfConstMat <- function(knots, shape, deg) {
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
