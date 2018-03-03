source("subor.R")

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
#' @param knots A vector of numeric that specifies the knots locations used in
#'     the linear splines model. The length of \code{knots} should be two less
#'     than the \code{NROW} of the matrices in \code{post}. For B-splines, ALL
#'     knots are required, including those outside the range. Required.
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
PlotLinSpline <- function(coefs, knots, limits, data, bases = "tpf") {

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
    } else if (bases == "bs" && is.numeric(knots)) {
        ## spline basis terms for b-splines
        model.mat <- splines::splineDesign(knots, x, ord = 2, outer.ok = TRUE)
    } else {
        stop("Incompatible given knots and bases type.")
    }
    browser()
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


## Plot spline

PlotSpline <- function(model, limits, data, fine = 200) {

    EPS <- 1e-6
    knots <- model$basis$knots
    knots.within <- knots[knots > (min(limits) - EPS) &
                          knots < (max(limits) + EPS)]
    type <- model$basis$type
    deg <- model$basis$degree
    coef.means <- model$means

    ## Calculate the coefficients for plotting
    coefs <- as.matrix(as.data.frame(coef.means)) + coef.means$population
    coefs[, NCOL(coefs)] <- coef.means$population
    colnames(coefs) <- names(coef.means)

    if (is.vector(coefs)) {
        coefs <- as.matrix(coefs)
    } else if (!is.matrix(coefs)) {
        stop("coefs must be a matrix")
    }

    if (!is.vector(knots) || !is.vector(limits)) {
        stop("knots and limits must be vectors")
    }

    if (!is.factor(data[[3]])) {
        data[[3]] <- factor(data[[3]], levels = colnames(coefs))
    }

    rownames(coefs) <- NULL
    names(knots) <- NULL

    ## "x" terms (explanatory variable)
    x <- c(seq(min(limits), max(limits), length.out = fine), knots.within)
    x <- x[order(x)]

    ## model matrix
    if (type == "tpf") {
        model.mat <- TpfDesign(x, knots, deg)$design
    } else if (type == "bs") {
        model.mat <- splines::splineDesign(knots, x, ord = deg + 1,
                                           outer.ok = TRUE)
    } else {
        stop("Unknown type of model.")
    }

    ## Calculate the response and melt them into a data frame
    my <- model.mat %*% coefs
    plot.data <- reshape2::melt(my)
    colnames(plot.data) <- c("points", "grps", "value")

    ## This is a reassurance step (reorder group levels)
    plot.data$grps <- factor(plot.data$grps, levels = colnames(coefs))

    plot.data$x <- rep(x, times = NCOL(coefs))
    pop.idx <- plot.data$grps %in% "population"

    if (missing(data)) {
        ggplot2::ggplot() +
            ggplot2::geom_line(aes(x, value, group = grps, col = grps),
                               plot.data[!pop.idx, ]) +
            ggplot2::geom_line(aes(x, value), plot.data[pop.idx, ],
                                   col = "black")
    } else if (is.data.frame(data)) {
        colnames(data) <- c("x", "y", "grps")
        ggplot2::ggplot() +
            ggplot2::geom_line(ggplot2::aes(x, value, group = grps, col = grps),
                               plot.data[!pop.idx, ]) +
            ggplot2::geom_point(ggplot2::aes(x, y, col = grps), data) +
            ggplot2::geom_line(ggplot2::aes(x, value), plot.data[pop.idx, ],
                                   col = "black") +
            ggplot2::geom_point(ggplot2::aes(x = knots.within,
                                             y = rep(0, length(knots.within))),
                                    col = "red", shape = 4)

    }
}

