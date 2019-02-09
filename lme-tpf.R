rm(list = ls())
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

x <- data[[1]] ## / max(data[[1]])
y <- data[[2]] ## / max(data[[2]])

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

## design matrices for the population curve
source("~/Dropbox/master/algo/subor.R")
des_ls_pop <- get_design_tpf(x, K_pop, deg)
raw_pop <- des_ls_pop$design[, 1:(deg + 1)]
G_pop <- solve(crossprod(raw_pop)) %*%
    crossprod(raw_pop, cbind(1, poly(x, degree = deg, simple = TRUE)))
X_pop <- raw_pop %*% G_pop
Z_pop <- des_ls_pop$design[, (deg + 2):(deg + 1 + K_pop)]

## design matrices for subject-specific curves
des_ls_sub <- get_design_tpf(x, K_sub, deg)
raw_sub <- des_ls_sub$design[, 1:(deg + 1)]
G_sub <- G_pop
X_sub <- raw_sub %*% G_sub
sub_cols <- seq(deg + 2, deg + 1 + K_sub, 1) # which columns of Z_pop is used for Z_sub
Z_sub <- des_ls_sub$design[, sub_cols]


library(nlme)
## covariance structures of random effects
pop_pd <- pdDiag(~ Z_pop)
## sub_pd <- pdBlocked(list(pdSymm(~ X_sub - 1), pdIdent(~ Z_sub - 1)))
sub_pd <- pdBlocked(list(pdDiag(~ X_sub - 1), pdDiag(~ Z_sub - 1)))

## fm <- lme(fixed = y ~ X_pop - 1,
##           random = list(ones = pop_pd, grp = sub_pd),
##           control = list(maxIter = 250, msMaxIter = 250, niterEM = 250,
##                          msVerbose = F, eval.max = 200))
## fm_struct <- fm$modelStruct$reStruct
## fm_cov <- lapply(as.matrix(fm_struct), function(x) {x * fm$sigma^2})
## sig2_eps <- fm$sigma^2
## sig2_u <- tail(diag(fm_cov$ones, names = FALSE), 1)
## sig2_eps / sig2_u


fm <- lme(fixed = y ~ X_pop - 1, random = list(ones = pop_pd))


## PLOT LME RESULTS
## the range of x to plot
plot_x <- c(seq(min(x), max(x), length = 200), des_ls_pop$knots, des_ls_sub$knots)
plot_x <- sort(unique(plot_x))

## design matrix corresponding to the range to plot
des_plot_pop <- get_design_tpf(plot_x, des_ls_pop$knots, deg = deg)
des_plot_sub <- get_design_tpf(plot_x, des_ls_sub$knots, deg = deg)
## C_mat <- cbind(des_plot_pop$design, des_plot_sub$design)
C_mat <- cbind(des_plot_pop$design[, 1:(deg + 1)] %*% G_pop,
               des_plot_pop$design[, (deg + 2):(deg + 1 + K_pop)],
               des_plot_sub$design[, 1:(deg + 1)] %*% G_sub,
               des_plot_sub$design[, sub_cols])

plot_y_pop <- C_mat %*% c(as.numeric(coef(fm, level = 1)),
                          rep(0, deg + 1 + length(sub_cols)))
plot_y_sub <- C_mat %*% t(as.matrix(coef(fm, level = 2)))
## plot_y_pop <- C_mat %*% c(as.numeric(nlme::fixef(fm)), rep(0, deg + 1 + K_sub))
## plot_y_sub <- C_mat %*% t(as.matrix(coef(fm)))
colnames(plot_y_sub) <- levels(grp)
## plot_y_pop <- des_plot_pop$design %*% as.numeric(coef(fm, level = 1))

## plotting from stan_fit object
## be careful of plot_x and plot_y
## fm_stan <- readRDS('~/Dropbox/master/stan/growth/no-penalty-full.rds')
## ch_no <- 2                              # which chain to choose
## betas <- get_posterior_mean(fm_stan)[1:3, ch_no]
## us <- get_posterior_mean(fm_stan)[4:(4+7), ch_no]
## bs <- matrix(get_posterior_mean(fm_stan)[(4+8):(4+8+m*3-1), ch_no], deg + 1)
## vs <- matrix(get_posterior_mean(fm_stan)[(4+8+m*3):(4+8+m*3+(m*8)-1), ch_no], K_sub)
## plot_y_pop <- C_mat %*% c(betas, us, rep(0, deg + 1 + K_sub))
## plot_y_sub <- C_mat %*% rbind(matrix(c(betas, us), deg + 1 + K_sub, m), rbind(bs, vs))
## colnames(plot_y_sub) <- levels(grp)

## construct suitable data frames for ggplot
ori_dat <- data.frame(x = x, y = y, sub = grp)
plot_dat_pop <- data.frame(x = plot_x, y = plot_y_pop)
plot_dat_sub <- reshape2::melt(plot_y_sub, varnames = c("x", "sub"), as.is = TRUE, value.name = "y")
plot_dat_sub$x <- plot_x
library(ggplot2)

## plot 10 sub at a time
sub_names <- unique(plot_dat_sub$sub)
num_pics <- length(sub_names)
for (i in seq(1, ceiling(length(sub_names) / num_pics))) {
    small_plot_dat <- subset(plot_dat_sub,
                             sub %in% sub_names[seq((i - 1) * num_pics + 1, i * num_pics)],
                             drop = TRUE)

    small_ori_dat <- subset(ori_dat,
                            sub %in% sub_names[seq((i - 1) * num_pics + 1, i * num_pics)],
                            drop = TRUE)

    ggobj <- ggplot(mapping = aes(x, y, col = sub)) +
        geom_point(data = small_ori_dat) +
        geom_line(aes(group = sub), data = small_plot_dat) +
        geom_line(aes(col = NULL), data = plot_dat_pop, lwd = 1) +
        ## geom_point(aes(x = des_ls_sub$knots, y = rep(0, length(des_ls_sub$knots)),
        ##                col = NULL), col = 'blue', pch = 2) +
        labs(title = paste("pop:", K_pop, "; sub:", K_sub, "; deg:", deg,
                           "; n:", length(levels(grp)))) +
        theme(legend.position="none")
    print(ggobj)
}


## manual EBLUPS
library(Matrix)
n <- length(y)
m <- length(levels(grp))
fm_struct <- fm$modelStruct$reStruct
fm_cov <- lapply(as.matrix(fm_struct), function(x) {x * fm$sigma^2})

X <- des_ls_pop$design[, 1:(deg + 1)]
Z_list <- list()
for (i in levels(grp)) {
    Z_list[[i]] <- des_ls_sub$design[grp == i, ]
}
Z <- cbind(des_ls_pop$design[, (deg + 2):(deg + 1 + K_pop)], bdiag(Z_list))
R <- Diagonal(n, fm$sigma^2)
G <- bdiag(list(fm_cov[[1]] * 100000, bdiag(rep(fm_cov[2], m))))

V <- tcrossprod(Z %*% G, Z) + R
Vi <- solve(V)
beta <- solve(crossprod(X, Vi %*% X)) %*% crossprod(X, Vi %*% y)
u <- tcrossprod(G, Z) %*% Vi %*% (y - X %*% beta)

## plots for manual EBLUPS
coef_pop <- c(as.numeric(beta), u[1:K_pop], rep(0, deg + 1 + K_sub))
coef_sub <- rbind(matrix(c(as.numeric(beta), u[1:K_pop]), deg + 1 + K_pop, m),
                  matrix(u[-(1:K_pop)], deg + 1 + K_sub))
plot_y_pop <- C_mat %*% coef_pop
plot_y_sub <- C_mat %*% coef_sub
colnames(plot_y_sub) <- levels(grp)

## calculate smoothing parameter
sig2_eps <- fm$sigma^2
sig2_u <- tail(diag(fm_cov$grp, names = FALSE), 1)
sig2_eps / sig2_u

## pop only eps : 0.001052391, u: 0.2678267 = 0.003929373
## random intercept eps: 0.0004935247, u: 0.0006089544 = 0.8104461
## full eps: 6.601439e-06, u: 4.500937 = 1.466681e-06



