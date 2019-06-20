library(tidyverse)
rm(list = ls())
source("~/Dropbox/master/algo/subor.R")

# generate some parameters
theta_original <- c(10, 13, 15.5, 17, 24.5, 25)
n_terms <- length(theta_original)
n_subs <- 5
knt <- seq(0, 10, 2)

## using hierarchical model
G <- cbind(1, poly(1:n_terms, deg = 1, raw = FALSE, simple = TRUE))
D <- get_diff_mat(n_terms, 1 + 1)
H <- crossprod(D, solve(tcrossprod(D)))
GH <- cbind(G, H)
theta <- c(solve(GH) %*% theta_original)   # (17.5, 13.267, -0.5, -1, 6, -7)
Amat <- get_constmat_bs(n_terms, "increasing") %*% GH

## using standard model
## theta <- theta_original
## Amat <- get_constmat_bs(n_terms, "increasing")

set.seed(1, "Mersenne-Twister")
delta <- t(mvtnorm::rmvnorm(n_subs, mean = rep(0, n_terms), sigma = diag(n_terms) * 10))
## delta <- t(tnorm::rmvtnorm(n_subs, mean = rep(0, n_terms), cov = diag(n_terms) * 10,
##                          initial = rep(0, n_terms),
##                          F = Amat, g = Amat %*% theta))

colnames(delta) <- paste0("sub", 1:n_subs)


## check if the parameters satisfy the constraints
all(Amat %*% theta >= 0)
all(Amat %*% (theta + delta) >= 0)

## add some noise
set.seed(2, "Mersenne-Twister")
x <- c(0, seq(0.5, 9.5, len = 10) + rnorm(10, 0, 0.1), 10)
Bmat <- get_design_bs(x, length(knt) - 2, 1)$design %*% GH # hierarchical
## Bmat <- get_design_bs(x, length(knt) - 2, 1)$design        # standard
y <- Bmat %*% (theta + delta) + rnorm(nrow(Bmat) * ncol(delta), 0, 0.5)
simdata <- as.tibble(cbind(x, y)) %>% gather(sub, y, -x)
## write_csv(simdata, "~/Dropbox/master/algo/data/simdata.csv")

## visualise the data
plot_x <- seq(0, 10, len = 300)

# hierarchical
plot_y <- get_design_bs(plot_x, length(knt) - 2, 1)$design %*% GH %*% (theta + delta)
plot_y_pop <- get_design_bs(plot_x, length(knt) - 2, 1)$design %*% GH %*% theta
# standard
## plot_y <- get_design_bs(plot_x, length(knt) - 2, 1)$design %*% (theta + delta)
## plot_y_pop <- get_design_bs(plot_x, length(knt) - 2, 1)$design %*% theta

plotdata <- as.tibble(cbind(plot_x, plot_y)) %>% gather(sub, plot_y, -plot_x)
plotdata_pop <- data_frame(plot_x, plot_y = as.numeric(plot_y_pop))
## write_csv(plotdata, "~/Dropbox/master/algo/data/simdata-curve.csv")
## write_csv(plotdata_pop, "~/Dropbox/master/algo/data/simdata-popcurve.csv")
## write_csv(plotdata, "~/Dropbox/master/algo/data/simdata-ridge.csv")
write_csv(plotdata, "~/Dropbox/master/algo/data/simdata-lmm.csv")
ggplot() +
    geom_line(aes(plot_x, plot_y, col = sub), plotdata) +
    geom_point(aes(x, y, col = sub), simdata) +
    geom_line(aes(plot_x, plot_y), plotdata_pop) +
    labs(x = 'x', y = 'y') +
    ylim(0, 35) +
    theme_bw() + theme(legend.position="none")



