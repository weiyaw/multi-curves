rm(list = ls())
setwd("~/Dropbox/master/algo")
source("main-ridge.R")
growth <- reshape2::melt(fda::growth[-3])
growth <- with(growth, data.frame(x = Var1, y = value, grp.sub = Var2, grp.pop = L1))
growth10 <- subset(growth,
                   grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                  "girl01", "girl02", "girl03", "girl04", "girl05"),
                   drop = TRUE)
growth10$grp.sub <- droplevels(growth10$grp.sub)

growth10boys <- subset(growth,
                       grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                      "boy06", "boy07", "boy08", "boy09", "boy10"),
                   drop = TRUE)
growth10boys$grp.sub <- droplevels(growth10boys$grp.sub)

growth20 <- subset(growth,
                   grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                  "boy06", "boy07", "boy08", "boy09", "boy10",
                                  "girl01", "girl02", "girl03", "girl04", "girl05",
                                  "girl06", "girl07", "girl08", "girl09", "girl10"),
                   drop = TRUE)
growth20$grp.sub <- droplevels(growth20$grp.sub)

growth30 <- subset(growth,
                   grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                  "boy06", "boy07", "boy08", "boy09", "boy10",
                                  "boy11", "boy12", "boy13", "boy14", "boy15",
                                  "girl01", "girl02", "girl03", "girl04", "girl05",
                                  "girl06", "girl07", "girl08", "girl09", "girl10",
                                  "girl11", "girl12", "girl13", "girl14", "girl15"),
                   drop = TRUE)
growth30$grp.sub <- droplevels(growth30$grp.sub)

growthboys <- subset(growth, grp.pop == "hgtm", drop = TRUE)

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
G_pop <- cbind(1, poly(1:n_bsf, deg = deg, raw = T, simple = T))
H_pop <- crossprod(D, solve(tcrossprod(D)))
B_pop <- des_ls_pop$design %*% cbind(G_pop, H_pop)
K_mat <- cbind(matrix(0, K, deg + 1), diag(K))

source("main-ridge.R")
set.seed(1)
fm1 <- bayes_ridge(y, cbind(X_pop, Z_pop), K_mat, 0, 10000)

source("main-ridge.R")
set.seed(24)
system.time(fm2_1 <- bayes_ridge_sub(y, grp, B_pop, K_mat, deg + 1, 1000, 10000))
system.time(fm2_2 <- bayes_ridge_sub(y, grp, B_pop, K_mat, deg + 1, 1000, 10000))
fm2_2$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = deg,
                    trans_mat = cbind(G_pop, H_pop))
fm2_2$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
source("graphs.R")
plot_spline(fm2_2)

plot(mcmc(1 / fm2_2$samples$precision$sub2))
plot(mcmc(1 / fm2_1$samples$precision$sub2))
gelman.diag(mcmc.list(mcmc(1 / fm2_1$samples$precision$sub2), mcmc(1 / fm2_2$samples$precision$sub2)))
gelman.plot(mcmc.list(mcmc(1 / fm2_1$samples$precision$sub2), mcmc(1 / fm2_2$samples$precision$sub2)))
effectiveSize(mcmc(1 / fm2_2$samples$precision$sub2))
effectiveSize(mcmc(fm2_2$samples$population[1, ]))

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

plot(fm4$samples$precision$eps)
plot(fm4$samples$precision$pop)
plot(fm4$samples$precision$sub2)
plot(fm4$samples$population[1, 900:1000])
plot(fm4$samples$population[2, 900:1000])

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
plot(fm6$samples$population[1, ])
plot(1 / fm6$samples$precision$pop, type = 'l')
sd(1 / fm6$samples$precision$sub2)


source("main-mc.R")
set.seed(21)
system.time(fm6_re <- mc_sub(y, grp, B_pop, K_mat, fm6$samples$precision, 1000))
fm6_re$basis <- list(type = 'tpf', knots = des_ls_pop$knots, degree = deg)
fm6_re$data <- data.frame(x = x, y = y, grp_sub = grp, grp_pop = NA)
source("graphs.R")
plot_spline(fm6_re)
plot(fm6_re$samples$subjects[1, 9, ])
