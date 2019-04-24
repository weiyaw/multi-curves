source("./subor.R")
stop("don't use this. use subor.R instead.")

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
        model.mat <- get_design_tpf(x, knots, deg)$design
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

    ## melt the response into a data frame
    ## plot_data_sub <- reshape2::melt(plot_y_sub, varnames = c("idx", "sub"),
    ##                                 as.is = TRUE, value.name = "y")

    ## plot_data_pop <- reshape2::melt(plot_y_pop, varnames = c("idx", "sub"),
    ##                                 as.is = TRUE, value.name = "y")


    ## if (type == "tpf") {
    ## ols_y <- get_ols(data$y[data$pop == plot_which[1]],
    ##                  get_design_tpf(data$x[data$pop == plot_which[1]], knots, deg)$design)
    ## } else if (type == "bs") {
    ## ols_y <- get_ols(data$y[data$pop == plot_which[1]],
    ##                  get_design_bs(data$x[data$pop == plot_which[1]], knots, deg)$design)
    ## } else {
    ##     stop("Unknown type of model.")
    ## }
    ## ols <- data.frame(x = plot_x,
    ##                   y = get_plot_y(ols_y))
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


plot_a_spline <- function(coef, limits) {


}

