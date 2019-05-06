# generate some parameters
library(tidyverse)
rm(list = ls())
source("~/Dropbox/master/algo/subor.R")
knt <- seq(0, 10, 2)
## theta <- c(5, 6, 15.5, 17, 24.5, 25)
theta <- c(10, 13, 15.5, 17, 24.5, 25)
n_terms <- length(theta)
n_subs <- 5
D <- get_diff_mat(n_terms, 1 + 1)
G <- cbind(1, poly(1:n_terms, deg = 1, raw = FALSE, simple = TRUE))
H <- crossprod(D, solve(tcrossprod(D)))
GH <- cbind(G, H)
Amat <- get_constmat_bs(n_terms, "increasing")
set.seed(1)
bvs <- t(tnorm::rmvtnorm(n_subs, mean = rep(0, n_terms), cov = diag(n_terms) * 10,
                         initial = rep(0, n_terms),
                         F = Amat %*% GH, g = Amat %*% theta))
delta <- GH %*% bvs
## delta <- t(mvtnorm::rmvnorm(n_subs, mean = rep(0, n_terms), sigma = diag(n_terms) * 10))
colnames(delta) <- paste0("sub", 1:n_subs)


## check if the parameters satisfy the constraints
all(Amat %*% theta >= 0)
all(Amat %*% (theta + delta) >= 0)

## add some noise
set.seed(2)
x <- c(0, seq(0.5, 9.5, len = 10) + rnorm(10, 0, 0.1), 10)
Bmat <- get_design_bs(x, length(knt) - 2, 1)$design
y <- Bmat %*% (theta + delta) + rnorm(nrow(Bmat) * ncol(delta), 0, 0.5)
simdata <- as.tibble(cbind(x, y)) %>% gather(sub, y, -x)

## visualise the data
plot_x <- seq(0, 10, len = 100)
plot_y <- get_design_bs(plot_x, length(knt) - 2, 1)$design %*% (theta + delta)
plot_y_pop <- get_design_bs(plot_x, length(knt) - 2, 1)$design %*% theta
plotdata <- as.tibble(cbind(plot_x, plot_y)) %>% gather(sub, plot_y, -plot_x)
ggplot() +
    geom_line(aes(plot_x, plot_y, col = sub), plotdata) +
    geom_point(aes(x, y, col = sub), simdata) +
    geom_line(aes(plot_x, plot_y_pop)) +
    theme_bw() + theme(legend.position="none")

## design matrices
rm(list = setdiff(ls(), c("simdata", "knt")))
source("~/Dropbox/master/algo/subor.R")
deg <- 1
K <- length(knt)
n_bsf <- K + deg + 1                # number of b-spline basis functions
D <- get_diff_mat(n_bsf, deg + 1)       # difference matrix
des_ls_pop <- get_design_bs(simdata$x, K, deg)
G_pop <- cbind(1, poly(1:n_bsf, deg = deg, raw = FALSE, simple = TRUE))
H_pop <- crossprod(D, solve(tcrossprod(D)))
B_pop <- des_ls_pop$design %*% cbind(G_pop, H_pop)
K_mat <- cbind(matrix(0, K, deg + 1), diag(K))
rm(list = c("deg", "K", "n_bsf", "D"))

source("~/Dropbox/master/algo/main-ridge.R")
set.seed(9)
fm1 <- bayes_ridge_sub(simdata$y, simdata$sub, B_pop, K_mat, 1 + 1, 0, 5000)
fm1$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = 1,
                  trans_mat = cbind(G_pop, H_pop))
fm1$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
plot_spline(fm1)
plot(fm1$samples$precision$pop)
plot(fm1$samples$population[1, ])

set.seed(1)
source("~/Dropbox/master/algo/main-mc.R")
fm1mc <- mc_sub(simdata$y, simdata$sub, B_pop, K_mat, fm1$samples$prec)
fm1mc$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = 1,
                  trans_mat = cbind(G_pop, H_pop))
fm1mc$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
source("~/Dropbox/master/algo/subor.R")
plot_spline(fm1mc)
plot(fm1mc$samples$precision$pop)
plot(fm1mc$samples$population[1, ])

source("~/Dropbox/master/algo/main-ridge.R")
set.seed(10)
fm1v2 <- bayes_ridge_sub_v2(simdata$y, simdata$sub, B_pop, K_mat, 1 + 1, 0, 5000)
fm1v2$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = 1,
                  trans_mat = cbind(G_pop, H_pop))
fm1v2$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
source("~/Dropbox/master/algo/subor.R")
plot_spline(fm1v2)
plot(fm1v2$samples$precision$pop)
plot(fm1v2$samples$population[1, ])

source("~/Dropbox/master/algo/main-mc.R")
set.seed(11)
fm1v2mc <- mc_sub(simdata$y, simdata$sub, B_pop, K_mat, fm1v2$samples$prec)
fm1v2mc$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = 1,
                  trans_mat = cbind(G_pop, H_pop))
fm1v2mc$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
source("~/Dropbox/master/algo/subor.R")
plot_spline(fm1v2mc)
plot(fm1v2mc$samples$precision$pop)
plot(fm1v2mc$samples$population[1, ])

## with monotonicity constraint
A_mat <- get_constmat_bs(NCOL(B_pop), "increasing") %*% cbind(G_pop, H_pop)
lower <- rep(0, NROW(A_mat))

source("~/Dropbox/master/algo/main-ridge.R")
set.seed(45)
fm1c <- bayes_ridge_cons_sub(simdata$y, simdata$sub, B_pop, K_mat, 1 + 1, A_mat,
                             0, 15000, init = fm1cv2$means)
plot(fm1c$samples$population[1, ])
plot(fm1c$samples$precision$pop)
plot(1 / fm1c$samples$precision$sub1[1, 1, ])
plot(1 / fm1c$samples$precision$sub2)
fm1c$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = 1,
                     trans_mat = cbind(G_pop, H_pop))
fm1c$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
plot_spline(fm1c)

source("~/Dropbox/master/algo/main-ridge.R")
set.seed(42)
fm1cv2 <- bayes_ridge_cons_sub_v2(simdata$y, simdata$sub, B_pop, K_mat, 1 + 1, A_mat,
                                0, 5000)

plot(fm1cv2$samples$population[1, ])
plot(fm1cv2$samples$subjects[1, 1, ])
plot(fm1cv2$samples$precision$pop)
plot(1 / fm1cv2$samples$precision$sub1[1, 1, ])
plot(1 / fm1cv2$samples$precision$sub2)
plot(1 / fm1cv2$samples$precision$eps)
fm1cv2$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = 1,
                     trans_mat = cbind(G_pop, H_pop))
fm1cv2$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
plot_spline(fm1cv2)


source("~/Dropbox/master/algo/main-mc.R")
set.seed(40)
fm1cmc <- mc_cons_sub(simdata$y, simdata$sub, B_pop, K_mat, A_mat, lower,
                      fm1c$samples$prec)
plot(fm1cmc$samples$population[1, ])
fm1cmc$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = 1,
                     trans_mat = cbind(G_pop, H_pop))
fm1cmc$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
plot_spline(fm1cmc)



