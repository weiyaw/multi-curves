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

## Get a penalised least squares estimate with a given response, design matrix
## and penalty matrix.
get_pls <- function(response, design, penalty) {
    inv_term <- chol2inv(chol(crossprod(design) + crossprod(penalty)))
    tcrossprod(inv_term, design) %*% as.vector(response)
}

## get the maximum-a-posteriori or maximum likelihood estimate from the model
get_max <- function(fm, type = "map") {
    if (type == "map") {
        if (is.null(fm$samples$lp)) {
            stop("Log-posterior not available.")
        } else {
            idx <- which.max(fm$samples$lp)
        }
    } else if (type == "mle") {
        if (is.null(fm$samples$ll)) {
            stop("Log-likelihood not available.")
        } else {
            idx <- which.max(fm$samples$ll)
        }
    } else {
        stop("Unknown type in get_max.")
    }
    list(population = fm$samples$population[, idx],
         subjects = fm$samples$subjects[, , idx],
         precision = list(pop = fm$samples$precision$pop[idx],
                          sub1 = fm$samples$precision$sub1[, , idx],
                          sub2 = fm$samples$precision$sub2[idx],
                          eps = fm$samples$precision$eps[idx]))
}

################################################################
##########################            ##########################
##########################    PLOT    ##########################
##########################            ##########################
################################################################

## Plot the population and subject curves (suitable for models with multiple populations)
## 'plot_which' takes a vector of names (for models with multiple populations)
## type can be either of "mean", "mle" or "map"
## If shade is true, thin the population samples by 10 and plot the thinned samples.

## Requires: model$basis$knots (including extrema), model$basis$type,
##   model$basis$degree, model$info$lvl_pop, model$data, model$mean,
##   model$samples$population (if shade is true)

plot_spline <- function(model, limits = NULL, plot_which = NULL, plot_type = "mean",
                        shade = FALSE) {

    EPS <- 1e-6
    fine <- 200                         # how fine the plot_x should be?
    knots <- model$basis$knots
    names(knots) <- NULL
    type <- model$basis$type
    deg <- model$basis$degree

    ## Rename the input data to prevent future accidental name change.
    data <- model$data
    colnames(data) <- c("x", "y", "sub", "pop")

    lvl_pop <- levels(data$pop)

    ## Use the range of the predictor if 'limits' is missing
    if (is.null(limits)) {
        limits <- range(data$x)
    } else if (!(is.vector(limits) && length(limits) == 2)) {
        stop("'limits' must be a vector of length 2.")
    }

    ## Plot all population if 'plot_which' is missing
    if (is.null(plot_which)) {
        plot_which <- lvl_pop
    } else if (!is.vector(plot_which)) {
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

    if (plot_type == "mean") {
        ## extract the posterior means
        coef_pop <- model$means$population
        dev_sub <- model$means$subjects
    } else if (plot_type == "map") {
        ## extract the map
        model_map <- get_max(model, "map")
        coef_pop <- model_map$population
        dev_sub <- model_map$subjects
    } else if (plot_type == "mle") {
        ## extract the mle (maximum likelihood)
        model_mle <- get_max(model, "mle")
        coef_pop <- model_mle$population
        dev_sub <- model_mle$subjects
    } else {
        stop("Unknown plotting type.")
    }

    ## Helper function to calculate the y axis of the plot
    get_plot_y <- function(coef) {
        model_mat %*% coef
    }

    ## y axis of the plot
    if (length(lvl_pop) == 0) {
        coef_sub <- dev_sub + coef_pop
        ## Align the data format of a single population model to a multiple
        ## population model. Create a dummy name for the population
        plot_which <- "dummy"
        data$pop <- "dummy"
        plot_y_pop <- list(dummy = get_plot_y(coef_pop))
        plot_y_sub <- list(dummy = get_plot_y(coef_sub))
    } else {
        coef_sub <- mapply(`+`, dev_sub[plot_which], coef_pop[plot_which],
                           SIMPLIFY = FALSE)
        plot_y_pop <- lapply(coef_pop[plot_which], get_plot_y)
        plot_y_sub <- lapply(coef_sub[plot_which], get_plot_y)
    }

    `%>%` <- tidyr::`%>%`               # define piping
    if (shade) {
        if (plot_which != "dummy") {
            stop("shade only works for models with single population.")
        }
        ## thin the sample by 10
        thin_idx <- seq(10, NCOL(model$samples$population), 10)
        thin_y <- get_plot_y(model$samples$population[, thin_idx])
        colnames(thin_y) <- thin_idx
        prob <- c(0.25, 0.75)
        prob_outer <- c(0.05, 0.95)
        plotdat_thin <- apply(thin_y, 1, function(x) quantile(x, c(prob, prob_outer))) %>%
            {tibble::as_tibble(t(.))} %>%
            dplyr::mutate(x = plot_x)
    }
    ## Reformat dataframe for ggplot
    plotdat_pop <- list()
    plotdat_sub <- list()
    for (i in plot_which) {
        plotdat_pop[[i]] <- data_frame(x = plot_x, y = as.numeric(plot_y_pop[[i]]))
        plotdat_sub[[i]] <- tibble::as_tibble(plot_y_sub[[i]]) %>%
            dplyr::mutate(x = plot_x) %>%
            tidyr::gather("sub", "y", -"x")
         }
    aes_ <- ggplot2::aes_
    geom_point <- ggplot2::geom_point
    geom_line <- ggplot2::geom_line

    ggls <- list()
    for (i in plot_which) {
        if (shade) {
            ggls$rib90 <- geom_ribbon(aes_(~x, ymin = ~`5%`, ymax = ~`95%`),
                                        plotdat_thin, fill = "grey95")
            ggls$rib50 <- geom_ribbon(aes_(~x, ymin = ~`25%`, ymax = ~`75%`),
                                      plotdat_thin, fill = "grey85")
            ## ggobj <- ggobj +
            ##     geom_ribbon(aes_(~x, ymin = ~`5%`, ymax = ~`95%`), plotdat_thin,
            ##             fill = "grey95") +
            ## geom_ribbon(aes_(~x, ymin = ~`25%`, ymax = ~`75%`), plotdat_thin,
            ##             fill = "grey85")
        }
        ggls$data <- geom_point(aes_(~x, ~y, col = ~sub), data[data$pop == i, ])
        ggls$sub <- geom_line(aes_(~x, ~y, col = ~sub, group = ~sub), plotdat_sub[[i]])
        ggls$pop <- geom_line(aes_(~x, ~y), data = plotdat_pop[[i]])
        ## ggobj <- ggobj +
        ##     geom_point(aes_(~x, ~y, col = ~sub), data[data$pop == i, ]) +
        ##     geom_line(aes_(~x, ~y, col = ~sub, group = ~sub), plotdat_sub[[i]]) +
        ##     geom_line(aes_(~x, ~y), data = plotdat_pop[[i]])
        print(ggplot2::ggplot() + ggls + theme_bw() + theme(legend.position="none"))
    }
    ggls
}












