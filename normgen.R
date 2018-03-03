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




new.rtmg <- function (n, M, r, initial, f = NULL, g = NULL, q = NULL, burn.in = 30)
{
    d = nrow(M)
    if (ncol(M) != d) {
        stop("Error: M must be a square matrix.")
    }
    if (length(initial) != d) {
        stop("Error: wrong length for initial value vector.")
    }

    M = (M + t(M))/2
    eigs = eigen(M, symmetric = T, only.values = T)$values

    if (any(eigs <= 0)) {
        stop("Error: M must be positive definite.")
    }
    R = chol(M)
    Mir = solve(M, r)
    Ri = solve(R)
    initial2 = as.vector(R %*% initial) - as.vector(r %*% Ri)

    if (is.null(f) & is.null(g)) {
        numlin = 0
        f2 = NULL
        g2 = NULL
    }
    else {
        if (is.matrix(f) & is.vector(g)) {
            numlin = nrow(f)
            if (length(g) != numlin | ncol(f) != d) {
                stop("Error: inconsistent linear constraints. f must be an m-by-d matrix and g an m-dimensional vector.")

            }
            if (any(f %*% initial + g <= 0)) {
                stop("Error: initial point violates linear constraints.")
            }
            f2 = f %*% Ri
            g2 = as.vector(f %*% Mir + g)
        }
        else {
            stop("Error: for linear constraints, f must be a matrix and g a vector.\n")
        }
    }
    if (is.null(q)) {
        numquad = 0
        quads = NULL
    }
    else {
        if (is.list(q)) {
            ll = lapply(q, length)
            nc = c(do.call("cbind", ll))
            if (any(nc != 3)) {
                stop("Error: each element in q must be a list of length 3.\n")
            }
            numquad = length(q)
            quads = matrix(nrow = numquad * (d + 2), ncol = d)
            for (i in 1:numquad) {
                qci = q[[i]]
                t = initial %*% qci[[1]] %*% initial + qci[[2]] %*%
                  initial + qci[[3]]
                if (t <= 0) {
                  stop("Error: initial point violates quadratic constraints. \n")
                }
                else {
                  A = qci[[1]]
                  B = qci[[2]]
                  C = qci[[3]]
                  quads[((i - 1) * (d + 2) + 1):((i - 1) * (d +
                    2) + d), ] = crossprod(Ri, A %*% Ri)
                  quads[((i - 1) * (d + 2) + d + 1), ] = 2 *
                    crossprod(Mir, A %*% Ri) + crossprod(B, Ri)
                  C = C + crossprod(Mir, A %*% Mir) + crossprod(B, Mir)
                  quads[i * (d + 2), ] = c(C, rep(0, d - 1))
                }
            }
        }
        else {
            stop("Error: for quadratic constraints, q must be a list.\n")
        }
    }
    seed = sample(1:1e+07, 1)
    samples = .Call("rtmg", n + burn.in, seed, initial2, numlin,
        f2, g2, numquad, quads, PACKAGE = "tmg")
    samples = samples[(burn.in + 1):(burn.in + n), ]
    samples = tcrossprod(samples, Ri) + matrix(rep(Mir, n), nrow = n,
        byrow = T)
    return(samples)
}
