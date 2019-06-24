# generate some parameters
library(tidyverse)
library(gridExtra)
library(bayesplot)
library(doMC)                           # for parallel loops
registerDoMC(4)                         # use 4 threads
source("~/Dropbox/master/algo/subor.R")
knt <- seq(0, 10, 2)
simdata <- read_csv("~/Dropbox/master/algo/data/simdata.csv")
plotdata <- read_csv("~/Dropbox/master/algo/data/simdata-curve.csv")
plotdata_pop <- read_csv("~/Dropbox/master/algo/data/simdata-popcurve.csv")
plotdata_ridge <- read_csv("~/Dropbox/master/algo/data/simdata-ridge.csv")
plotdata_lmm <- read_csv("~/Dropbox/master/algo/data/simdata-lmm.csv")

## compara gen-data from either ridge and LMM
ridge_curve <- ggplot() +
    geom_line(aes(plot_x, plot_y, col = sub), plotdata_ridge) +
    geom_line(aes(plot_x, plot_y), plotdata_pop) +
    labs(x = 'x', y = 'y') +
    ylim(0, 35) +
    theme_bw() + theme(legend.position="none")

lmm_curve <- ggplot() +
    geom_line(aes(plot_x, plot_y, col = sub), plotdata_lmm) +
    geom_line(aes(plot_x, plot_y), plotdata_pop) +
    labs(x = 'x', y = 'y') +
    ylim(0, 35) +
    theme_bw() + theme(legend.position="none")

ridge_vs_lmm <- grid.arrange(ridge_curve, lmm_curve, ncol = 2)

## ggsave("~/Dropbox/master/thesis/images/ridge-vs-lmm.pdf", ridge_vs_lmm,
##        width = 10, height = 6)


## generate the truth (dotted)
true_curve <- list(sub = geom_line(aes(plot_x, plot_y, col = sub), plotdata, lty = 2,
                                   alpha = 0.9),
                   pop = geom_line(aes(plot_x, plot_y), plotdata_pop, lty = 2),
                   theme_bw(),
                   theme(legend.position="none"),
                   labs(x = 'x', y = 'y'))

## design matrices (B-spline LMM model)
rm(list = setdiff(ls(), c("simdata", "knt", "true_curve")))
source("~/Dropbox/master/algo/subor.R")
deg <- 1
K <- length(knt) - 2                    # inner knots
n_bsf <- K + deg + 1                    # number of b-spline basis functions
D <- get_diff_mat(n_bsf, deg + 1)       # difference matrix
type <- "bs"                            # "bs", "bs-ridge" or "tpf"

## design matrices, Gmat and Hmat
des_info <- get_design_bs(simdata$x, K, deg)                         # bs
Gmat <- cbind(-1/sqrt(n_bsf), poly(1:n_bsf, deg = deg, raw = FALSE)) # bs unraw
Hmat <- crossprod(D, solve(tcrossprod(D)))                           # bs
Bmat <- des_info$design %*% cbind(Gmat, Hmat)                        # bs
Kmat <- cbind(matrix(0, K, deg + 1), diag(K))                        # bs
rm(list = c("K", "n_bsf", "D"))

source("~/Dropbox/master/algo/main-ridge.R")
set.seed(100, kind = "L'Ecuyer-CMRG")
fm1_ls <- foreach(i = 1:4) %dopar% {
    init <- list(pop = get_pls(simdata$y, Bmat, Kmat) + rnorm(NCOL(Bmat), sd = 100))
    fm <- bayes_ridge_sub(simdata$y, simdata$sub, Bmat, Kmat, deg + 1, 1000, 2000,
                          init = init)
    fm$basis <- list(type = 'bs_hier', knots = des_info$knots, degree = deg, # bs
                     trans_mat = cbind(Gmat, Hmat))
    fm$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
    fm
}
RNGkind("Mersenne-Twister")

source("~/Dropbox/master/algo/subor.R")
g1ls <- list()
for (fm in fm1_ls) {
    g1ls <- c(g1ls, plot_spline(fm, shade = TRUE))
}
g1curve_all <- ggplot() + g1ls[grep("pop|sub|data", names(g1ls))] + theme_bw() +
    theme(legend.position="none")

## ggsave("~/Dropbox/master/thesis/images/truth-gibbs-all.pdf", g1curve_all,
##        width = 5, height = 6)

source("~/Dropbox/master/algo/diagnostic.R")
## visualise the truth curve and add it to g1curve
fm1 <- do.call(combine_fm, fm1_ls)
g1curve_true <- ggplot() + true_curve + plot_spline(fm1, shade = FALSE)

## regression curve
## ggsave("~/Dropbox/master/thesis/images/truth-gibbs-true.pdf", g1curve_true,
##        width = 5, height = 6)

## diagnostic
library(xtable)
source("~/Dropbox/master/algo/diagnostic.R")
flat1 <- do.call(flatten_chains, fm1_ls)
long1 <- summary_matrix_flats(flat1)
short1 <- long1 %>% filter(Rhat > 1.01 | n_eff < 500)
short1 <- long1 %>% filter(grepl("theta\\[1\\]|delta\\[1,1\\]", Parameter))
## print(xtable(short1), include.rownames=FALSE, tabular.environment = "tabular")

## diagnostic plots
g1combo_theta1 <- mcmc_combo(flat1, pars = "theta[1]", c("dens_overlay", "trace"))
g1combo_delta11 <- mcmc_combo(flat1, pars = "delta[1,1]", c("dens_overlay", "trace"))
g1combo_sum11 <- as.tibble(flat1[, , "theta[1]"] + flat1[, , "delta[1,1]"]) %>%
    gather(Chain, "theta[1] + delta[1,1]") %>%
    separate(Chain, c(NA, "Chain"), " ") %>%
    mcmc_combo(c("dens_overlay", "trace"))

## ggsave("~/Dropbox/master/thesis/images/truth-gibbs-combo-theta1.pdf", g1combo_theta1,
##        width = 10, height = 3)
## ggsave("~/Dropbox/master/thesis/images/truth-gibbs-combo-delta11.pdf", g1combo_delta11,
##        width = 10, height = 3)
## ggsave("~/Dropbox/master/thesis/images/truth-gibbs-combo-sum11.pdf", g1combo_sum11,
##        width = 10, height = 3)


## VERSION 2 OF GIBBS
source("~/Dropbox/master/algo/main-ridge.R")
set.seed(101, kind = "L'Ecuyer-CMRG")
fm1v2_ls <- foreach(i = 1:4) %dopar% {
    init <- list(pop = get_pls(simdata$y, Bmat, Kmat) + rnorm(NCOL(Bmat), sd = 100))
    fm <- bayes_ridge_sub_v2(simdata$y, simdata$sub, Bmat, Kmat, deg + 1, 1000, 2000,
                       init = init)
    if (type == "tpf") {
        fm$basis <- list(type = 'tpf', knots = des_info$knots, degree = deg) # tpf
    } else if (type == "bs") {
        fm$basis <- list(type = 'bs_hier', knots = des_info$knots, degree = deg, # bs
                         trans_mat = cbind(Gmat, Hmat))
    } else if (type == "bs-ridge") {
        fm$basis <- list(type = 'bs', knots = des_info$knots, degree = deg) # bs-ridge
    }
    fm$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
    fm
}
RNGkind("Mersenne-Twister")

source("~/Dropbox/master/algo/subor.R")
g1v2ls <- list()
for (fm in fm1v2_ls) {
    g1v2ls <- c(g1v2ls, plot_spline(fm, shade = TRUE, silent = TRUE))
}
g1v2curve_all <- ggplot() + g1v2ls[grep("pop|sub|data", names(g1v2ls))] + theme_bw() +
    theme(legend.position="none")

## ggsave("~/Dropbox/master/thesis/images/truth-gibbsv2-all.pdf", g1v2curve_all,
##        width = 5, height = 6)

source("~/Dropbox/master/algo/diagnostic.R")
## visualise the truth curve and add it to g1v2curve
fm1v2 <- do.call(combine_fm, fm1v2_ls)
g1v2curve_true <- ggplot() + true_curve + plot_spline(fm1v2, shade = FALSE, silent = TRUE)

## regression curve
## ggsave("~/Dropbox/master/thesis/images/truth-gibbsv2-true.pdf", g1v2curve_true,
##        width = 5, height = 6)


## diagnostic
source("~/Dropbox/master/algo/diagnostic.R")
flat1v2 <- do.call(flatten_chains, fm1v2_ls)
long1v2 <- summary_matrix_flats(flat1v2)
## print(xtable(long1v2), include.rownames=FALSE, tabular.environment = "longtable")
short1v2 <- long1v2 %>% filter(Rhat > 1.01 | n_eff < 500)
short1v2 <- long1v2 %>% filter(grepl("theta\\[1\\]|delta\\[1,1\\]", Parameter))
## print(xtable(short1v2), include.rownames=FALSE, tabular.environment = "tabular")

## diagnostic plots
g1v2combo_theta1 <- mcmc_combo(flat1v2, pars = "theta[1]", c("dens_overlay", "trace"))
g1v2combo_delta11 <- mcmc_combo(flat1v2, pars = "delta[1,1]", c("dens_overlay", "trace"))
g1v2combo_sum11 <- as.tibble(flat1v2[, , "theta[1]"] + flat1v2[, , "delta[1,1]"]) %>%
    gather(Chain, "theta[1] + delta[1,1]") %>%
    separate(Chain, c(NA, "Chain"), " ") %>%
    mcmc_combo(c("dens_overlay", "trace"))

## ggsave("~/Dropbox/master/thesis/images/truth-gibbsv2-combo-theta1.pdf",
##        g1v2combo_theta1, width = 10, height = 3)
## ggsave("~/Dropbox/master/thesis/images/truth-gibbsv2-combo-delta11.pdf",
##        g1v2combo_delta11, width = 10, height = 3)
## ggsave("~/Dropbox/master/thesis/images/truth-gibbsv2-combo-sum11.pdf",
##        g1v2combo_sum11, width = 10, height = 3)


## TRY DIFFERENT MODELS, change the model parameters at the beginning and rerun v2
## bs-ridge
## ggsave("~/Dropbox/master/thesis/images/truth-bsridge-combo-theta1.pdf",
##        g1v2combo_theta1, width = 10, height = 3.5)
## ggsave("~/Dropbox/master/thesis/images/truth-bsridge-combo-delta11.pdf",
##        g1v2combo_delta11, width = 10, height = 3.5)
## ggsave("~/Dropbox/master/thesis/images/truth-bsridge-combo-sum11.pdf",
##        g1v2combo_sum11, width = 10, height = 3.5)

## tpf
## ggsave("~/Dropbox/master/thesis/images/truth-tpf-combo-theta1.pdf",
##        g1v2combo_theta1, width = 10, height = 3.5)
## ggsave("~/Dropbox/master/thesis/images/truth-tpf-combo-delta11.pdf",
##        g1v2combo_delta11, width = 10, height = 3.5)
## ggsave("~/Dropbox/master/thesis/images/truth-tpf-combo-sum11.pdf",
##        g1v2combo_sum11, width = 10, height = 3.5)











## with monotonicity constraint
## with monotonicity constraint
## with monotonicity constraint
## with monotonicity constraint

## design matrices (tpf of B-spline LMM model)
rm(list = setdiff(ls(), c("simdata", "knt", "true_curve")))
source("~/Dropbox/master/algo/subor.R")
deg <- 1
K <- length(knt) - 2                    # inner knots
n_bsf <- K + deg + 1                    # number of b-spline basis functions
D <- get_diff_mat(n_bsf, deg + 1)       # difference matrix
type <- "tpf"                            # "bs" or "tpf"

# generate some parameters
library(tidyverse)
library(gridExtra)
library(bayesplot)
library(doMC)                           # for parallel loops
registerDoMC(4)                         # use 4 threads
source("~/Dropbox/master/algo/subor.R")
knt <- seq(0, 10, 2)
simdata <- read_csv("~/Dropbox/master/algo/data/simdata.csv")
plotdata <- read_csv("~/Dropbox/master/algo/data/simdata-curve.csv")
plotdata_pop <- read_csv("~/Dropbox/master/algo/data/simdata-popcurve.csv")

## generate the truth (dotted)
true_curve <- list(sub = geom_line(aes(plot_x, plot_y, col = sub), plotdata, lty = 2,
                                   alpha = 0.9),
                   pop = geom_line(aes(plot_x, plot_y), plotdata_pop, lty = 2),
                   theme_bw(),
                   theme(legend.position="none"),
                   labs(x = 'x', y = 'y'))

## design matrices, Gmat and Hmat
if (type == "tpf") {
    des_info <- get_design_tpf(simdata$x, K, deg) # tpf
    Bmat <- des_info$design                       # tpf
    Kmat <- cbind(matrix(0, K, deg + 1), diag(K)) # tpf
    Amat <- get_constmat_tpf(des_info$knots, "increasing", deg)
    lower <- rep(0, NROW(Amat))
} else if (type == "bs") {
    des_info <- get_design_bs(simdata$x, K, deg)                         # bs
    Gmat <- cbind(-1/sqrt(n_bsf), poly(1:n_bsf, deg = deg, raw = FALSE)) # bs unraw
    Hmat <- crossprod(D, solve(tcrossprod(D)))                           # bs
    Bmat <- des_info$design %*% cbind(Gmat, Hmat)                        # bs
    Kmat <- cbind(matrix(0, K, deg + 1), diag(K))                        # bs
    Amat <- get_constmat_bs(NCOL(Bmat), "increasing") %*% cbind(Gmat, Hmat)
    lower <- rep(0, NROW(Amat))
}
rm(list = c("K", "n_bsf", "D"))

library(nlme)
pop <- rep(1, length(simdata$sub))
Xmat <- unname(Bmat[ , 1:(deg + 1)])
Zmat <- unname(Bmat[, -(1:(deg + 1))])
fit <- lme(y ~ Xmat - 1, simdata, list(pop = pdIdent(~Zmat - 1),
                                       sub = pdBlocked(list(pdSymm(~Xmat - 1),
                                                            pdIdent(~Zmat - 1)))))
prec <- lapply(as.matrix(fit$modelStruct$reStruct), function(x) x * fit$sigma^2)
prec$pop <- 1 / diag(prec$pop)[[1]]
prec$sub1 <- unname(solve(prec$sub[1:(deg + 1), 1:(deg + 1)]))
prec$sub2 <- 1 / diag(prec$sub)[[deg + 2]]
prec$sub <- NULL

source("~/Dropbox/master/algo/main-ridge.R")
set.seed(103, kind = "L'Ecuyer-CMRG")
fm1cv2_ls <- foreach(i = 1:4) %dopar% {
    init <- list(pop = c(tnorm::rmvtnorm(1, mean = get_pls(simdata$y, Bmat, Kmat),
                                         initial = c(1, 1, rep(0, 4)),
                                         F = Amat, g = -1 * lower)))
    fm <- bayes_ridge_cons_sub_v2(simdata$y, simdata$sub, Bmat, Kmat, deg + 1,
                                  Amat, 1000, 2000, init, prec = prec)
    if (type == "tpf") {
        fm$basis <- list(type = 'tpf', knots = des_info$knots, degree = deg) # tpf
    } else if (type == "bs") {
        fm$basis <- list(type = 'bs_hier', knots = des_info$knots, degree = deg, # bs
                         trans_mat = cbind(Gmat, Hmat))
    }
    fm$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
    fm
}
RNGkind("Mersenne-Twister")

source("~/Dropbox/master/algo/subor.R")
g1cv2ls <- list()
for (fm in fm1cv2_ls) {
    g1cv2ls <- c(g1cv2ls, plot_spline(fm, shade = TRUE, silent = TRUE))
}
g1cv2curve_all <- ggplot() + g1cv2ls[grep("pop|sub|data", names(g1cv2ls))] + theme_bw() +
    theme(legend.position="none")

## ggsave("~/Dropbox/master/thesis/images/truth-gibbscv2-all.pdf", g1cv2curve_all,
##        width = 5, height = 6)

source("~/Dropbox/master/algo/diagnostic.R")
## visualise the truth curve and add it to g1cv2curve
fm1cv2 <- do.call(combine_fm, fm1cv2_ls)
g1cv2curve_true <- ggplot() + true_curve +
    plot_spline(fm1cv2, shade = FALSE, silent = TRUE)

## regression curve
## ggsave("~/Dropbox/master/thesis/images/truth-gibbscv2-true.pdf", g1cv2curve_true,
##        width = 5, height = 6)


## diagnostic
source("~/Dropbox/master/algo/diagnostic.R")
flat1cv2 <- do.call(flatten_chains, fm1cv2_ls)
long1cv2 <- summary_matrix_flats(flat1cv2)
## print(xtable(long1cv2), include.rownames=FALSE, tabular.environment = "longtable")
short1cv2 <- long1cv2 %>% filter(Rhat > 1.01 | n_eff < 500)
short1cv2 <- long1cv2 %>% filter(grepl("theta\\[1\\]|delta\\[1,1\\]", Parameter))
## print(xtable(short1cv2), include.rownames=FALSE, tabular.environment = "tabular")

## diagnostic plots
g1cv2combo_theta1 <- mcmc_combo(flat1cv2, pars = "theta[1]", c("dens_overlay", "trace"))
g1cv2combo_delta11 <- mcmc_combo(flat1cv2, pars = "delta[1,1]", c("dens_overlay", "trace"))



source("~/Dropbox/master/algo/diagnostic.R")
flat1cv2 <- do.call(flatten_chains, list(fm1cv2))[, , 1:43, drop = FALSE]
tail(summary_matrix_flats(flat1cv2), n = 10)
mcmc_combo(flat1cv2[, , 1:43], pars = "theta[1]", c("dens", "trace"))



