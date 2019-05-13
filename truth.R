# generate some parameters
library(tidyverse)
library(gridExtra)
library(bayesplot)
library(doMC)                           # for parallel loops
registerDoMC(4)                         # use 4 threads
source("~/Dropbox/master/algo/subor.R")
knt <- seq(0, 10, 2)
simdata <- read_csv("~/Dropbox/master/algo/data/simdata.csv")

## design matrices (B-spline hierarchical model)
rm(list = setdiff(ls(), c("simdata", "knt")))
source("~/Dropbox/master/algo/subor.R")
deg <- 1
K <- length(knt) - 2                    # inner knots
n_bsf <- K + deg + 1                    # number of b-spline basis functions
D <- get_diff_mat(n_bsf, deg + 1)       # difference matrix
type <- "bs"                            # "bs", "bs-ridge" or "tpf"
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
    bayes_ridge_sub(simdata$y, simdata$sub, Bmat, Kmat, deg + 1, 1000, 2000,
                    init = init)
}
RNGkind("Mersenne-Twister")

source("~/Dropbox/master/algo/diagnostic.R")
source("~/Dropbox/master/algo/subor.R")
fm <- do.call(combine_fm, fm1_ls)
if (type == "tpf") {
    fm$basis <- list(type = 'tpf', knots = des_info$knots, degree = deg) # tpf
} else if (type == "bs") {
    fm$basis <- list(type = 'bs_hier', knots = des_info$knots, degree = deg, # bs
                     trans_mat = cbind(Gmat, Hmat))
} else if (type == "bs-ridge") {
    fm$basis <- list(type = 'bs', knots = des_info$knots, degree = deg) # bs-ridge
}
fm$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
g1curve <- plot_spline(fm, shade = TRUE) + theme_bw() + theme(legend.position="none")

## regression curve
ggsave("~/Dropbox/master/thesis/images/truth-curve-gibbs.pdf", g1curve,
       width = 7, height = 5)

## diagnostic
source("~/Dropbox/master/algo/diagnostic.R")
flat1 <- do.call(flatten_chains, fm1_ls)
long1 <- summary_matrix_flats(flat1)
## print(xtable(short_sum), include.rownames=FALSE, tabular.environment = "longtable")
short1 <- long1 %>% filter(Rhat > 1.01 | n_eff < 500)
g1diag <- grid.arrange(mcmc_trace(flat1, "theta[1]"),
                       mcmc_dens_overlay(flat1, "theta[1]"),
                       nrow = 1)

## diagnostic plots
ggsave("~/Dropbox/master/thesis/images/truth-diag-gibbs.pdf", g1diag,
       width = 14, height = 5)

source("~/Dropbox/master/algo/main-ridge.R")
set.seed(101, kind = "L'Ecuyer-CMRG")
fm1v2_ls <- foreach(i = 1:4) %dopar% {
    init <- list(pop = get_pls(simdata$y, Bmat, Kmat) + rnorm(NCOL(Bmat), sd = 100))
    bayes_ridge_sub_v2(simdata$y, simdata$sub, Bmat, Kmat, deg + 1, 1000, 2000,
                       init = init)
}
RNGkind("Mersenne-Twister")

source("~/Dropbox/master/algo/diagnostic.R")
source("~/Dropbox/master/algo/subor.R")
fm <- do.call(combine_fm, fm1v2_ls)
## for (fm in fm1v2_ls) {
    if (type == "tpf") {
        fm$basis <- list(type = 'tpf', knots = des_info$knots, degree = deg) # tpf
    } else if (type == "bs") {
        fm$basis <- list(type = 'bs_hier', knots = des_info$knots, degree = deg, # bs
                         trans_mat = cbind(Gmat, Hmat))
    } else if (type == "bs-ridge") {
        fm$basis <- list(type = 'bs', knots = des_info$knots, degree = deg) # bs-ridge
    }
    fm$data <- simdata %>% mutate(grp_sub = sub, grp_pop = NA, sub = NULL)
    plot_spline(fm, shade = TRUE)
## }

## diagnostic
source("~/Dropbox/master/algo/diagnostic.R")
flat1v2 <- do.call(flatten_chains, fm1v2_ls)
long1v2 <- summary_matrix_flats(flat1v2)
## print(xtable(short_sum), include.rownames=FALSE, tabular.environment = "longtable")
short1v2 <- long1v2 %>% filter(Rhat > 1.01 | n_eff < 500)
grid.arrange(mcmc_trace(flat1v2, "theta[1]"),
             mcmc_dens_overlay(flat1v2, "theta[1]"),
             nrow = 1)






## with monotonicity constraint
A_mat <- get_constmat_bs(NCOL(Bmat), "increasing") %*% cbind(Gmat, Hmat)
lower <- rep(0, NROW(A_mat))

source("~/Dropbox/master/algo/main-ridge.R")
set.seed(45)
fm1c <- bayes_ridge_cons_sub(simdata$y, simdata$sub, Bmat, Kmat, 1 + 1, A_mat,
                             0, 5000)

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




