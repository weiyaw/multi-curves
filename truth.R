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

true_curve <- ggplot() +
    geom_line(aes(plot_x, plot_y, col = sub), plotdata, lty = 2, alpha = 0.9) +
    geom_line(aes(plot_x, plot_y), plotdata_pop, lty = 2) +
    theme_bw() + theme(legend.position="none") +
    labs(x = 'x', y = 'y')

## design matrices (B-spline hierarchical model)
rm(list = setdiff(ls(), c("simdata", "knt", "true_curve")))
source("~/Dropbox/master/algo/subor.R")
deg <- 1
K <- length(knt) - 2                    # inner knots
n_bsf <- K + deg + 1                    # number of b-spline basis functions
D <- get_diff_mat(n_bsf, deg + 1)       # difference matrix
type <- "tpf"                            # "bs", "bs-ridge" or "tpf"
raw <- FALSE                             # raw poly or not

## design matrices, Gmat and Hmat
if (type == "tpf") {
    des_info <- get_design_tpf(simdata$x, K, deg) # tpf
    Bmat <- des_info$design                       # tpf
    Kmat <- cbind(matrix(0, K, deg + 1), diag(K)) # tpf
} else if (type == "bs" & raw) {
    des_info <- get_design_bs(simdata$x, K, deg)           # bs
    Gmat <- cbind(1, poly(1:n_bsf, deg = deg, raw = TRUE)) # bs raw
    Hmat <- crossprod(D, solve(tcrossprod(D)))             # bs
    Bmat <- des_info$design %*% cbind(Gmat, Hmat)          # bs
    Kmat <- cbind(matrix(0, K, deg + 1), diag(K))          # bs
} else if (type == "bs" & !raw) {
    des_info <- get_design_bs(simdata$x, K, deg)                         # bs
    Gmat <- cbind(-1/sqrt(n_bsf), poly(1:n_bsf, deg = deg, raw = FALSE)) # bs unraw
    Hmat <- crossprod(D, solve(tcrossprod(D)))                           # bs
    Bmat <- des_info$design %*% cbind(Gmat, Hmat)                        # bs
    Kmat <- cbind(matrix(0, K, deg + 1), diag(K))                        # bs
} else if (type == "bs-ridge") {
    des_info <- get_design_bs(simdata$x, K, deg)  # bs-ridge
    Bmat <- des_info$design                       # bs-ridge
    Kmat <- D                                     # bs-ridge
}
rm(list = c("K", "n_bsf", "D"))

source("~/Dropbox/master/algo/main-ridge.R")
set.seed(100, kind = "L'Ecuyer-CMRG")
fm1_ls <- foreach(i = 1:4) %dopar% {
    init <- list(pop = get_pls(simdata$y, Bmat, Kmat) + rnorm(NCOL(Bmat), sd = 100))
    fm <- bayes_ridge_sub(simdata$y, simdata$sub, Bmat, Kmat, deg + 1, 1000, 2000,
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
g1ls <- list()
for (fm in fm1_ls) {
    g1ls <- c(g1ls, plot_spline(fm, shade = TRUE))
}
g1curve_all <- ggplot() + g1ls[grep("pop|sub|data", names(g1ls))] + theme_bw() +
    theme(legend.position="none")

## ggsave("~/Dropbox/master/thesis/images/truth-gibbs-all.pdf", g1curve_all,
##        width = 5, height = 7)

source("~/Dropbox/master/algo/diagnostic.R")
## visualise the truth curve and add it to g1curve
fm1 <- do.call(combine_fm, fm1_ls)
g1curve_true <- true_curve + plot_spline(fm1, shade = FALSE)

## regression curve
## ggsave("~/Dropbox/master/thesis/images/truth-gibbs-true.pdf", g1curve_true,
##        width = 5, height = 7)

## diagnostic
library(xtable)
source("~/Dropbox/master/algo/diagnostic.R")
flat1 <- do.call(flatten_chains, fm1_ls)
long1 <- summary_matrix_flats(flat1)
short1 <- long1 %>% filter(Rhat > 1.01 | n_eff < 500)
short1 <- long1 %>% filter(grepl("theta\\[1\\]|delta\\[1,1\\]", Parameter))
print(xtable(short1), include.rownames=FALSE, tabular.environment = "tabular")

## diagnostic plots
g1combo_theta1 <- mcmc_combo(flat1, pars = "theta[1]", c("dens_overlay", "trace"))
g1combo_delta11 <- mcmc_combo(flat1, pars = "delta[1,1]", c("dens_overlay", "trace"))
g1combo_sum11 <- as.tibble(flat1[, , "theta[1]"] + flat1[, , "delta[1,1]"]) %>%
    gather(Chain, "theta[1] + delta[1,1]") %>%
    separate(Chain, c(NA, "Chain"), " ") %>%
    mcmc_combo(c("dens_overlay", "trace"))

## ggsave("~/Dropbox/master/thesis/images/truth-gibbs-combo-theta1.pdf", g1combo_theta1,
##        width = 10, height = 3.5)
## ggsave("~/Dropbox/master/thesis/images/truth-gibbs-combo-delta11.pdf", g1combo_delta11,
##        width = 10, height = 3.5)
## ggsave("~/Dropbox/master/thesis/images/truth-gibbs-combo-sum11.pdf", g1combo_sum11,
##        width = 10, height = 3.5)


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
    g1v2ls <- c(g1v2ls, plot_spline(fm, shade = TRUE))
}
g1v2curve_all <- ggplot() + g1v2ls[grep("pop|sub|data", names(g1v2ls))] + theme_bw() +
    theme(legend.position="none")

## ggsave("~/Dropbox/master/thesis/images/truth-gibbsv2-all.pdf", g1v2curve_all,
##        width = 5, height = 7)
## ggsave("~/Dropbox/master/thesis/images/truth-bsridge-all.pdf", g1v2curve_all,
##        width = 5, height = 7)
## ggsave("~/Dropbox/master/thesis/images/truth-tpf-all.pdf", g1v2curve_all,
##        width = 5, height = 7)


source("~/Dropbox/master/algo/diagnostic.R")
## visualise the truth curve and add it to g1v2curve
fm1v2 <- do.call(combine_fm, fm1v2_ls)
g1v2curve_true <- true_curve + plot_spline(fm1v2, shade = FALSE)

## regression curve
## ggsave("~/Dropbox/master/thesis/images/truth-gibbsv2-true.pdf", g1v2curve_true,
##        width = 5, height = 7)
## ggsave("~/Dropbox/master/thesis/images/truth-bsridge-true.pdf", g1v2curve_true,
##        width = 5, height = 7)
## ggsave("~/Dropbox/master/thesis/images/truth-tpf-true.pdf", g1v2curve_true,
##        width = 5, height = 7)


## diagnostic
source("~/Dropbox/master/algo/diagnostic.R")
flat1v2 <- do.call(flatten_chains, fm1v2_ls)
long1v2 <- summary_matrix_flats(flat1v2)
## print(xtable(short_sum), include.rownames=FALSE, tabular.environment = "longtable")
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
##        g1v2combo_theta1, width = 10, height = 3.5)
## ggsave("~/Dropbox/master/thesis/images/truth-gibbsv2-combo-delta11.pdf",
##        g1v2combo_delta11, width = 10, height = 3.5)
## ggsave("~/Dropbox/master/thesis/images/truth-gibbsv2-combo-sum11.pdf",
##        g1v2combo_sum11, width = 10, height = 3.5)


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
A_mat <- get_constmat_bs(NCOL(Bmat), "increasing") %*% cbind(Gmat, Hmat)
lower <- rep(0, NROW(A_mat))

source("~/Dropbox/master/algo/main-ridge.R")
set.seed(45)
fm1c <- bayes_ridge_cons_sub(simdata$y, simdata$sub, Bmat, Kmat, 1 + 1, A_mat,
                             1000, 2000)

fm1c$basis <- list(type = 'bs_hier', knots = des_info$knots, degree = 1,
                     trans_mat = cbind(Gmat, Hmat))
fm1c$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
plot_spline(fm1c)

source("~/Dropbox/master/algo/main-ridge.R")
set.seed(42)
fm1cv2 <- bayes_ridge_cons_sub_v2(simdata$y, simdata$sub, Bmat, Kmat, 1 + 1, A_mat,
                                0, 5000)

fm1cv2$basis <- list(type = 'bs_hier', knots = des_info$knots, degree = 1,
                     trans_mat = cbind(Gmat, Hmat))
fm1cv2$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
source("~/Dropbox/master/algo/subor.R")
plot_spline(fm1cv2)




