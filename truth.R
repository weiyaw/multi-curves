## generate some parameters
rm(list = ls())
source("~/Dropbox/master/algo/subor.R")
knt <- seq(0, 10, 2)
theta <- c(0, 1, 1.5, 3, 4.5, 5)
n_terms <- length(theta)
n_subs <- 5
Amat <- get_constmat_bs(n_terms, "increasing")
set.seed(1)
delta <- t(tnorm::rmvtnorm(n_subs, mean = rep(0, n_terms), cov = diag(n_terms) * 10,
                           initial = rep(0, n_terms),
                           F = Amat, g = Amat %*% theta))
colnames(delta) <- paste0("sub", 1:n_subs)


## check if the parameters satisfy the constraints
all(Amat %*% theta >= 0)
all(Amat %*% (theta + delta) >= 0)

## add some noise
set.seed(2)
x <- c(0, seq(0.5, 9.5, len = 10) + rnorm(10, 0, 0.1), 10)
Bmat <- get_design_bs(x, length(knt) - 2, 1)$design
y <- Bmat %*% (theta + delta) + rnorm(nrow(Bmat) * ncol(delta), 0, 0.1)
simdata <- as.tibble(cbind(x, y)) %>% gather(sub, y, -x)

## visualise the data
plot_x <- seq(0, 10, len = 100)
plot_y <- get_design_bs(plot_x, length(knt) - 2, 1)$design %*% (theta + delta)
plot_y_pop <- get_design_bs(plot_x, length(knt) - 2, 1)$design %*% theta
plotdata <- as.tibble(cbind(plot_x, plot_y)) %>% gather(sub, plot_y, -plot_x)
ggplot() +
    geom_line(aes(plot_x, plot_y, col = sub), plotdata) +
    geom_point(aes(x, y, col = sub), simdata) +
    geom_line(aes(plot_x, plot_y_pop))

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
fm1 <- bayes_ridge_sub(simdata$y, simdata$sub, B_pop, K_mat, 1 + 1, 0, 5000)

fm1$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = 1,
                  trans_mat = cbind(G_pop, H_pop))
fm1$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
source("~/Dropbox/master/algo/subor.R")
plot_spline(fm1)

source("~/Dropbox/master/algo/main-mc.R")
A_mat <- get_constmat_bs(NCOL(B_pop), "increasing") %*% cbind(G_pop, H_pop)
lower <- rep(0, NROW(A_mat))
fm1cmc <- mc_cons_sub(simdata$y, simdata$sub, B_pop, K_mat, A_mat, lower,
                      fm1$samples$prec, 1000)
fm1cmc$basis <- list(type = 'bs_hier', knots = des_ls_pop$knots, degree = 1,
                     trans_mat = cbind(G_pop, H_pop))
fm1cmc$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
plot_spline(fm1cmc)



