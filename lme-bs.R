rm(list = ls())
library(nlme)
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


data <- growth10boys
deg <- 2

## 'x' invariant under scale and location transformation
## 'y' differ by a scale after scale transformation
x <- data[[1]] / max(data[[1]])
y <- data[[2]] / max(data[[2]])

## number of inner knots
K_pop <- 8
K_sub <- 8

## convert the group variable into a factor
if (is.factor(data[[3]])) {
    grp <- droplevels(data[[3]])
} else {
    grp <- factor(data[[3]], levels = unique(data[[3]]))
}

## dummy variable of ones
ones <- rep(1, length(grp))
m <- length(unique(grp))

source("~/Dropbox/master/algo/subor.R")
n_bsf <- K_pop + deg + 1                # number of b-spline basis functions
D <- get_diff_mat(n_bsf, deg + 1)       # difference matrix

## design matrices for the population curve
des_ls_pop <- get_design_bs(x, K_pop, deg)
B_pop <- des_ls_pop$design
G_pop <- cbind(1, poly(1:n_bsf, deg = deg, raw = F, simple = T))
H_pop <- crossprod(D, solve(tcrossprod(D)))
X_pop <- B_pop %*% G_pop
Z_pop <- B_pop %*% H_pop

fm <- lme(fixed = y ~ X_pop - 1,
          random = list(ones = pdIdent(~ Z_pop - 1)))

## ## design matrix corresponding to the range to plot (population model)
## B_plot_pop <- get_design_bs(plot_x, des_ls_pop$knots, deg = deg)$design

## plot_y_pop <- B_plot_pop %*% G_pop %*% as.numeric(fixef(fm)) +
##     B_plot_pop %*% H_pop %*% as.numeric(ranef(fm))

## library(ggplot2)
## ori_dat <- data.frame(x = x, y = y, sub = grp)
## plot_dat_pop <- data.frame(x = plot_x, y = plot_y_pop)

## ggplot() + geom_point(aes(x, y, col = sub), ori_dat) +
##     geom_line(aes(x, y), plot_dat_pop)

## design matrices for subject-specific curves
des_ls_sub <- get_design_bs(x, K_sub, deg)
B_sub <- des_ls_sub$design
G_sub <- G_pop
H_sub <- H_pop
X_sub <- B_sub %*% G_sub
Z_sub <- B_sub %*% H_sub


## covariance structures of random effects
pop_pd <- pdIdent(~ Z_pop - 1)
sub_pd <- pdBlocked(list(pdSymm(~ X_sub - 1), pdIdent(~ Z_sub - 1)))

fm <- lme(fixed = y ~ X_pop - 1,
          random = list(ones = pop_pd, grp = sub_pd))


## PLOT LME RESULTS
## the range of x to plot
plot_x <- c(seq(min(x), max(x), length = 200),
            des_ls_pop$knots[(deg + 1):(deg + 1 + K_pop)])
plot_x <- sort(unique(plot_x))


## design matrix corresponding to the range to plot (ss model)
B_plot_pop <- get_design_bs(plot_x, des_ls_pop$knots, deg = deg)$design
B_plot_sub <- get_design_bs(plot_x, des_ls_sub$knots, deg = deg)$design
C_mat <- cbind(B_plot_pop %*% G_pop, B_plot_pop %*% H_pop,
               B_plot_sub %*% G_sub, B_plot_sub %*% H_sub)
plot_y_pop <- C_mat %*% c(as.numeric(coef(fm, level = 1)), rep(0, deg + 1 + K_sub))
plot_y_sub <- C_mat %*% t(as.matrix(coef(fm, level = 2)))
## plot_y_pop <- C_mat %*% c(as.numeric(fixef(fm)), rep(0, deg + 1 + K_sub))
## plot_y_sub <- C_mat %*% t(as.matrix(coef(fm)))
colnames(plot_y_sub) <- levels(grp)

## plotting from stan_fit object
## be careful of plot_x and plot_y
fm_stan <- readRDS('~/Dropbox/master/stan/growth/bs-quad.rds')
ch_no <- 2                              # which chain to choose
betas <- get_posterior_mean(fm_stan)[1:3, ch_no]
us <- get_posterior_mean(fm_stan)[4:(4+7), ch_no]
bs <- matrix(get_posterior_mean(fm_stan)[(4+8):(4+8+m*3-1), ch_no], deg + 1)
vs <- matrix(get_posterior_mean(fm_stan)[(4+8+m*3):(4+8+m*3+(m*8)-1), ch_no], K_sub)
plot_y_pop <- C_mat %*% c(betas, us, rep(0, n_bsf))
plot_y_sub <- C_mat %*% rbind(matrix(c(betas, us), n_bsf, m), rbind(bs, vs))
colnames(plot_y_sub) <- levels(grp)


library(ggplot2)
ori_dat <- data.frame(x = x, y = y, sub = grp)
plot_dat_pop <- data.frame(x = plot_x, y = plot_y_pop)
plot_dat_sub <- reshape2::melt(plot_y_sub, varnames = c("x", "sub"), as.is = TRUE, value.name = "y")
plot_dat_sub$x <- plot_x

ggplot(mapping = aes(x, y, col = sub)) +
    geom_point(data = ori_dat) +
    geom_line(aes(group = sub), data = plot_dat_sub) +
    geom_line(aes(col = NULL), data = plot_dat_pop, lwd = 1) +
    xlab("Age (scaled)") + ylab("Height (scaled)") +
    ## labs(title = paste("pop:", K_pop, "; sub:", K_sub, "; deg:", deg,
    ##                    "; n:", length(levels(grp)))) +
    theme_bw() +
    theme(legend.position="none")







## marginal cov
library(rstan)

Bmat <- cbind(X_pop, Z_pop)
K <- diag(c(rep(0, deg + 1), rep(1, NCOL(B_pop) - deg - 1)))
## Bmat <- B_pop
## K <- D
stan_data <- list(n_terms = NCOL(Bmat), n_samples = length(y), n_penalty = NROW(K),
                  Bmat = Bmat, y = y, K = K)
## stan_data <- list(n_terms = NCOL(B_pop), n_samples = length(y), n_penalty = NROW(D),
##                   Bmat = B_pop, y = y, K = D)

fm1 <- stan(file = "~/Dropbox/master/algo/marginal_cov.stan",
           data = stan_data, iter = 3000, warmup = 2000, chains = 1,
           cores = 4, algo = "NUTS",
           control = list(adapt_delta = 0.98, max_treedepth = 16), seed = 1)


## fm2 <- vb(c_file, data = stan_data)

shinystan::launch_shinystan(fm1)

source("main-bs.R")
set.seed(1)
var_samples <- marginal_cov(5000, 0, data = list(Bmat = B_pop, y = y, K = D))
plot(var_samples$eps)

log_like <- function(pop, eps) {
    n_terms <- NCOL(D)
    n_samples <- length(y)
    xBmat <- crossprod(Bmat)
    xK <- crossprod(K)
    xy <- crossprod(y)
    Bmatxy <- crossprod(B_pop, y)

    L <- xK / pop + xBmat / eps

    if (all(pop > 0, eps > 0)) {
        -0.5 * (n_terms * log(pop) + n_samples * log(eps)) -
            0.5 / eps * (crossprod(Bmatxy, solve(L, Bmatxy)) - xy)
    } else {
        stop()
        -Inf
    }
}
mapply(log_like, rep(0.001, 1, 100), MoreArgs=list(eps = 0.01))
