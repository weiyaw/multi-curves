rm(list = ls())
setwd("~/Dropbox/master/algo")
source("main-ridge.R")
library(tidyverse)

sitka10 <- read_table2("~/Dropbox/master/algo/data/sitka.txt") %>%
    filter(id.num >= 1 & id.num <= 10) %>%
    mutate(id.num = as.factor(id.num)) %>%
    select(days, log.size, id.num)

growth10boys <- fda::growth[-3] %>%
    map_dfc(function(x) as.tibble(x, rownames = "age")) %>%
    gather(sub, height, -age, -age1) %>%
    filter(grepl("^boy0[1-9]|^boy10", sub)) %>%
    mutate(age = as.numeric(age), sub = as.factor(sub)) %>%
    select(age, height, sub)

data <- growth10boys
deg <- 2

## 'x' invariant under scale and location transformation
## 'y' differ by a scale after scale transformation
x <- data[[1]] / max(data[[1]])
y <- data[[2]] / max(data[[2]])

## number of interior knots
K <- 8

## convert the group variable into a factor
if (is.factor(data[[3]])) {
    grp <- droplevels(data[[3]])
} else {
    grp <- factor(data[[3]], levels = unique(data[[3]]))
}

source("~/Dropbox/master/algo/subor.R")

## design matrices for the population curve (b-spline with transformation)
n_bsf <- K + deg + 1                # number of b-spline basis functions
D <- get_diff_mat(n_bsf, deg + 1)       # difference matrix
des_ls_pop <- get_design_bs(x, K, deg)
G_pop <- cbind(1, poly(1:n_bsf, deg = deg, raw = FALSE, simple = TRUE))
H_pop <- crossprod(D, solve(tcrossprod(D)))
B_pop <- des_ls_pop$design %*% cbind(G_pop, H_pop)
K_mat <- cbind(matrix(0, K, deg + 1), diag(K))

source("main-ridge.R")
set.seed(1)
fm1 <- bayes_ridge(y, cbind(X_pop, Z_pop), K_mat, 0, 10000)

source("main-ridge.R")
set.seed(24)
system.time(fm2_1 <- bayes_ridge_sub(y, grp, B_pop, K_mat, deg + 1, 0, 1000))
system.time(fm2_2 <- bayes_ridge_sub(y, grp, B_pop, K_mat, deg + 1, 1000, 5000))
system.time(fm2_3 <- bayes_ridge_sub(y, grp, B_pop, K_mat, deg + 1, 1000, 5000))
system.time(fm2_4 <- bayes_ridge_sub(y, grp, B_pop, K_mat, deg + 1, 1000, 5000))
fm2_2$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = deg,
                    trans_mat = cbind(G_pop, H_pop))
fm2_2$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
source("subor.R")
plot_spline(fm2_2)

plot(fm2$samples$precision$pop)
acf(fm2$samples$precision$pop)
source("diagnostic.R")
ess_rfun(flatten_chains(fm2_1, fm2_2)[, , "sig2_delta2"])

source("main-ridge.R")
A_mat <- get_constmat_bs(NCOL(B_pop), "increasing") %*% cbind(G_pop, H_pop)
set.seed(2)
fm2c <- bayes_ridge_cons_sub(y, grp, B_pop, K_mat, deg + 1, A_mat, 0, 5000)
fm2c$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = deg,
                    trans_mat = cbind(G_pop, H_pop))
fm2c$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
plot_spline(fm2c)
plot(fm2c$samples$population[1, ])

source("main-mc.R")
lower <- rep(0, NROW(A_mat))
fm2cmc <- mc_cons_sub(y, grp, B_pop, K_mat, A_mat, lower, fm2_1$samples$prec, 1000)
fm2cmc <- mc_cons_sub(y, grp, B_pop, K_mat, A_mat, lower, fm2_1$samples$prec, 1000)
fm2cmc$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = deg,
                     trans_mat = cbind(G_pop, H_pop))
fm2cmc$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
plot_spline(fm2cmc)

source("main-mc.R")
set.seed(21)
system.time(fm2_1_re <- mc_sub(y, grp, cbind(X_pop, Z_pop), K_mat, fm2_1$samples$precision, 1000))
fm2_1_re$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = deg,
                       trans_mat = cbind(G_pop, H_pop))
fm2_1_re$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
source("graphs.R")
plot_spline(fm2_1_re)

effectiveSize(mcmc(1 / fm2_1_re$samples$precision$eps))
effectiveSize(mcmc(1 / fm2_1_re$samples$precision$pop))
effectiveSize(mcmc(1 / fm2_1_re$samples$precision$sub2))
effectiveSize(mcmc(t(fm2_1_re$samples$population)))
plot(mcmc(t(fm2_1_re$samples$subject[, 2, ])))

## design matrices for the population curve (b-spline)
n_bsf <- K + deg + 1                # number of b-spline basis functions
des_ls_pop <- get_design_bs(x, K, deg)
B_pop <- des_ls_pop$design

K_mat <- get_diff_mat(n_bsf, deg + 1)

source("main-ridge.R")
set.seed(1)
fm3 <- bayes_ridge(y, B_pop, K_mat, 0, 5000)
plot_x <- seq(min(x), max(x), len = 100)
plot_y <- get_design_bs(plot_x, des_ls_pop$knots, deg)$design %*% fm3$means$population
ggplot() + geom_point(aes(x, y, col = grp)) + geom_line(aes(plot_x, plot_y))

source("main-ridge.R")
set.seed(2)
fm4 <- bayes_ridge_sub(y, grp, B_pop, K_mat, deg + 1, 0, 500)
fm4$basis <- list(type = 'bs', knots = des_ls_pop$knots, degree = deg)
fm4$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
source("graphs.R")
plot_spline(truncate_spline(fm4, 400))

source("main-ridge.R")
A_mat <- get_constmat_bs(NCOL(B_pop), "increasing")
set.seed(2)
fm4c <- bayes_ridge_cons_sub(y, grp, B_pop, K_mat, deg + 1, A_mat, 0, 10000)
fm4c$basis <- list(type = 'bs', knots = des_ls_pop$knots, degree = deg)
fm4c$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
plot_spline(fm4c)
plot(fm4c$samples$population[2, ], ylim = c(-2.15, -2.3))
plot(fm4c$samples$precision$pop)
plot(fm4c$samples$precision$sub1[1, 1, ])
plot(fm4c$samples$precision$sub2)



source("main-mc.R")
set.seed(21)
system.time(fm4_re <- mc_sub(y, grp, B_pop, K_mat, fm4$samples$precision, 500))
fm4_re$basis <- list(type = 'bs', knots = des_ls_pop$knots, degree = deg)
fm4_re$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
source("graphs.R")
plot_spline(fm4_re)

fm_stan <- readRDS("../stan/growth/bs-ridge.rds")
shinystan::launch_shinystan(fm_stan)

## design matrices for the population curve (tpf)
deg <- 2
des_ls_pop <- get_design_tpf(x, K, deg)
B_pop <- des_ls_pop$design
K_mat <- cbind(matrix(0, K, deg + 1), diag(K))

source("main-ridge.R")
set.seed(1)
fm5 <- bayes_ridge(y, B_pop, K_mat, 0, 10000)
plot_x <- seq(min(x), max(x), len = 100)
plot_y <- get_design_tpf(plot_x, des_ls_pop$knots, deg)$design %*% fm5$means$population
ggplot() + geom_point(aes(x, y, col = grp)) + geom_line(aes(plot_x, plot_y))


source("main-ridge.R")
set.seed(2)
system.time(fm6 <- bayes_ridge_sub(y, grp, B_pop, K_mat, deg + 1, 0, 5000))
fm6$basis <- list(type = 'tpf', knots = des_ls_pop$knots, degree = deg)
fm6$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
source("graphs.R")
plot_spline(fm6)

source("main-ridge.R")
A_mat <- get_constmat_tpf(des_ls_pop$knots, "increasing", deg)
set.seed(2)
fm6c <- bayes_ridge_cons_sub(y, grp, B_pop, K_mat, deg + 1, A_mat, 0, 5000)
plot(fm6c$samples$population[1, ])
plot(fm6c$samples$population[2, ])
plot(fm6c$samples$subjects[1, 1, ])
plot(1 / fm6c$samples$precision$pop, type = 'l')
plot(1 / fm6c$samples$precision$sub2)
fm6c$basis <- list(type = 'tpf', knots = des_ls_pop$knots, degree = deg)
fm6c$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
plot_spline(fm6c)


source("main-mc.R")
set.seed(21)
system.time(fm6_re <- mc_sub(y, grp, B_pop, K_mat, fm6$samples$precision, 1000))
fm6_re$basis <- list(type = 'tpf', knots = des_ls_pop$knots, degree = deg)
fm6_re$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
source("graphs.R")
plot_spline(fm6_re)
plot(fm6_re$samples$subjects[1, 9, ])
