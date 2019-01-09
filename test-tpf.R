
rm(list = ls())
setwd("~/Dropbox/master/algo/")
source("main-tpf.R")
source("graphs.R")

sitka <- read.table("data/sitka.txt", header = T)
sitka <- with(sitka, data.frame(x = days / 674,
                                y = log.size,
                                grp.sub = factor(id.num),
                                grp.pop = ozone))

sitka5 <- read.table("data/sitka5.txt", header = T)
sitka5 <- with(sitka5, data.frame(x = days / 674,
                                  y = log.size,
                                  grp.sub = factor(id.num),
                                  grp.pop = ozone))

sitka10 <- read.table("data/sitka10.txt", header = T)
sitka10 <- with(sitka10, data.frame(x = days / 674,
                                    y = log.size,
                                    grp.sub = factor(id.num),
                                    grp.pop = ozone))

growth <- reshape2::melt(fda::growth[-3])
growth <- with(growth, data.frame(x = Var1, y = value, grp.sub = Var2, grp.pop = L1))
growth10 <- subset(growth,
                   grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                  "girl01", "girl02", "girl03", "girl04", "girl05"),
                   drop = TRUE)
growth10$grp.sub <- droplevels(growth10$grp.sub)

growth20 <- subset(growth,
                   grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                  "girl01", "girl02", "girl03", "girl04", "girl05",
                                  "boy06", "boy07", "boy08", "boy09", "boy10",
                                  "girl06", "girl07", "girl08", "girl09", "girl10"),
                   drop = TRUE)
growth20$grp.sub <- droplevels(growth20$grp.sub)


source("main-tpf.R")
## linear spline
fm1_1 <- lme_tpf(sitka10, K = 5, deg = 1)
source("graphs.R")
plot_spline(fm1_1)
fm2_1 <- lme_tpf(growth10, K = 8, deg = 1)
plot_spline(fm2_1)

## quadratic spline (algorithm diverges)
## only works if the random effect of the quadratic polynomial term is removed
fm1_2 <- lme_tpf(sitka10, K = 5, deg = 2)
plot_spline(fm1_2)
fm2_2 <- lme_tpf(growth10, K = 8, deg = 2)
plot_spline(fm2_2)


data <- growth10
K_pop <- 8
K_sub <- 7
deg <- 2

x <- data[[1]] / max(data[[1]])
y <- data[[2]] / max(data[[2]])

## convert the group variable into a factor
if (is.factor(data[[3]])) {
    grp <- droplevels(data[[3]])
} else {
    grp <- factor(data[[3]], levels = unique(data[[3]]))
}

## dummy variable of ones
ones <- rep(1, length(grp))

## design matrices for the population curve
des_ls_pop <- get_design_tpf(x, K_pop, deg)
X_pop <- des_ls_pop$design[, 1:(deg + 1)]
Z_pop <- des_ls_pop$design[, (deg + 2):(deg + 1 + K_pop)]

## design matrices for subject-specific curves
des_ls_sub <- get_design_tpf(x, K_sub, deg)
X_sub <- des_ls_sub$design[, 1:(deg + 1)]
Z_sub <- des_ls_sub$design[, (deg + 2):(deg + 1 + K_sub)]

## covariance structures of random effects
pop_pd <- nlme::pdIdent(~ Z_pop - 1)
sub_pd <- nlme::pdBlocked(list(nlme::pdSymm(~ X_sub - 1), nlme::pdIdent(~ Z_sub - 1)))

fm <- nlme::lme(fixed = y ~ X_pop - 1, random = list(ones = pop_pd, grp = sub_pd),
                control = list(maxIter = 50, msMaxIter = 150, niterEM = 150))

## the range of x to plot
plot_x <- c(seq(min(data$x), max(data$x), length = 200), des_ls_pop$knots, des_ls_sub$knots)
plot_x <- sort(unique(plot_x))

## design matrix corresponding to the range to plot
des_plot_pop <- get_design_tpf(plot_x, des_ls_pop$knots, deg = deg)
des_plot_sub <- get_design_tpf(plot_x, des_ls_sub$knots, deg = deg)
C_mat <- cbind(des_plot_pop$design, des_plot_sub$design)
plot_y_pop <- C_mat %*% c(as.numeric(coef(fm, level = 1)), rep(0, deg + 1 + K_sub))
plot_y_sub <- C_mat %*% t(as.matrix(coef(fm, level = 2)))
colnames(plot_y_sub) <- levels(grp)

## construct suitable data frames for ggplot
ori_dat <- with(data, data.frame(x = x, y = y, sub = grp.sub))
plot_dat_pop <- data.frame(x = plot_x, y = plot_y_pop)
plot_dat_sub <- reshape2::melt(plot_y_sub, varnames = c("x", "sub"), as.is = TRUE, value.name = "y")
plot_dat_sub$x <- plot_x
library(ggplot2)
ggplot(mapping = aes(x, y, col = sub)) +
    geom_point(data = ori_dat) +
    geom_line(aes(group = sub), data = plot_dat_sub) +
    geom_line(aes(col = NULL), data = plot_dat_pop)


## Fit a lme model (quadratic)
## sub.pd.q <- pdBlocked(list(pdSymm(~ time.q + I(time.q^2)), pdIdent(~ Z.q - 1)))
## sub.pd.q <- pdSymm(~ time.q + I(time.q^2))
## sub.pd.q <- pdBlocked(list(pdSymm(~ time.q), pdIdent(~ Z.q - 1)))
## quad.fm <- lme(fixed = y ~ time.q + I(time.q^2),
##                random = list(pop.level = pop.pd.q, sub.level = sub.pd.q))



## TEST MIXED MODEL LINEAR SPLINE (SubjectsTpf)
source("main-tpf.R")
## system.time(fm2 <- get_tpf_old(sitka, 5, 2, size = 100, burn = 0))
saveRDS(fm2, "tpf-long.rds")
source("graphs.R")
fm1 <- readRDS("simulations/single/tpf-lin.rds")
fm2 <- readRDS("simulations/single/tpf-quad.rds")
PlotSpline(fm1, range(sitka$x), sitka)
PlotSpline(fm2, range(sitka$x), sitka)

fm3 <- readRDS("simulations/bspline-lin.rds")
fm4 <- readRDS("simulations/bspline-quad.rds")
PlotSpline(fm3, range(sitka$x), sitka)
PlotSpline(fm4, range(sitka$x), sitka)

## TEST MIXED MODEL LINEAR SPLINE WITH MULTIPLE POPULATION (SubjectsTpfMul)
## LOAD DATA sitka10 and sitka
rm(list = ls())


source("main-tpf.R")
source("graphs.R")
## single population (sitka10)
system.time(fm1 <- SubjectsTpf(sitka10, 5, deg = 2, shape = "increasing", size = 1000))
plot_spline(fm1)

## single population (sitka)
system.time(fm2 <- SubjectsTpf(sitka, 5, deg = 2, shape = "increasing", size = 10000, burn = 0))
plot_spline(fm2)
plot(fm2$samples$population[1, ])
plot(fm2$samples$population[2, ])
plot(fm2$samples$population[3, ])
plot(fm2$samples$population[4, ])
plot(fm2$samples$population[5, ])
plot(fm2$samples$population[6, ])
plot(fm2$samples$population[7, ])
plot(fm2$samples$population[8, ])

plot(fm2$samples$subjects[1, 1, ])
plot(fm2$samples$subjects[2, 1, ])
plot(fm2$samples$subjects[3, 1, ])
plot(fm2$samples$subjects[4, 1, ])
plot(fm2$samples$subjects[5, 1, ])
plot(fm2$samples$subjects[6, 1, ])
plot(fm2$samples$subjects[7, 1, ])
plot(fm2$samples$subjects[8, 1, ])


## multiple population (sitka)
system.time(fm3 <- SubjectsTpfMul(sitka, 5, deg = 2, shape = "increasing", size = 2000, burn = 0))

plot(fm1$samples$population$`1`[1, ])
plot(fm1$samples$population$`1`[2, ])
plot(fm1$samples$population$`1`[3, ])
plot(fm1$samples$population$`1`[4, ])
plot(fm1$samples$population$`1`[5, ])
plot(fm1$samples$population$`1`[6, ])
plot(fm1$samples$population$`1`[7, ])
plot(fm1$samples$population$`1`[8, ])

plot(fm1$samples$subjects$`1`[1, 1, ])
plot(fm1$samples$subjects$`1`[2, 1, ])
plot(fm1$samples$subjects$`1`[3, 1, ])
plot(fm1$samples$subjects$`1`[4, 1, ])
plot(fm1$samples$subjects$`1`[5, 1, ])
plot(fm1$samples$subjects$`1`[6, 1, ])
plot(fm1$samples$subjects$`1`[7, 1, ])
plot(fm1$samples$subjects$`1`[8, 1, ])

fm2 <- readRDS("simulations/multi/multi-1k-bugless.rds")
## multiple population (1000 burn, independent start, garbage)
fm3 <- readRDS("simulations/multi/multi-10k.rds")
source("graphs.R")
plot_spline(fm1)
plot_spline(fm2)
plot_spline(fm3)
plot_spline(fm4)

## multiple population (1000 burn, previous start)
source("main-tpf.R")
system.time(fm4 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 1000, burn = 0))

plot(fm4$samples$population$`1`[1, ])
plot(fm4$samples$population$`1`[2, ])
plot(fm4$samples$population$`1`[3, ])
plot(fm4$samples$population$`1`[4, ])
plot(fm4$samples$population$`1`[5, ])
plot(fm4$samples$population$`1`[6, ])
plot(fm4$samples$population$`1`[7, ])
plot(fm4$samples$population$`1`[8, ])

plot(fm4$samples$subjects$`1`[1, 1, ])
plot(fm4$samples$subjects$`1`[2, 1, ])
plot(fm4$samples$subjects$`1`[3, 1, ])
plot(fm4$samples$subjects$`1`[4, 1, ])
plot(fm4$samples$subjects$`1`[5, 1, ])
plot(fm4$samples$subjects$`1`[6, 1, ])
plot(fm4$samples$subjects$`1`[7, 1, ])
plot(fm4$samples$subjects$`1`[8, 1, ])

source("main-tpf.R")
system.time(fm4 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 1000, burn = 0))

fm5 <- readRDS("simulations/multi/multi-sitka.rds")
system.time(fm5 <- SubjectsTpfMul(sitka10, 5, deg = 2, shape = "increasing", size = 10000, burn = 0))

source("graphs.R")
fm7t <- truncate(fm7, 2000)
plot_spline(fm7t)

fm2t <- truncate(fm2, 300)
plot_spline(fm2t)


#### BERKELY GROWTH DATASET ####
fm6 <- readRDS("simulations/multi/multi-growth.rds")
fm7 <- readRDS("simulations/multi/multi-growth-uncon.rds")
growth <- reshape2::melt(fda::growth[-3])
growth <- with(growth, data.frame(x = Var1 / max(Var1), y = value / max(value), grp.sub = Var2, grp.pop = L1))
growth10 <- subset(growth,
                   grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                  "girl01", "girl02", "girl03", "girl04", "girl05"),
                   drop = TRUE)
growth10$grp.sub <- droplevels(growth10$grp.sub)

growth20 <- subset(growth,
                   grp.sub %in% c("boy01", "boy02", "boy03", "boy04", "boy05",
                                  "girl01", "girl02", "girl03", "girl04", "girl05",
                                  "boy06", "boy07", "boy08", "boy09", "boy10",
                                  "girl06", "girl07", "girl08", "girl09", "girl10"),
                   drop = TRUE)
growth20$grp.sub <- droplevels(growth20$grp.sub)




set.seed(1)
source("main-tpf.R")
system.time(fm8 <- sub_tpf(growth10, 8, deg = 2, penalty = FALSE, shape = "increasing", size = 10000, burn = 0, verbose = T))

plot(fm8$samples$population[1, ])
plot(fm8$samples$population[2, ])
plot(fm8$samples$population[3, ])
plot(fm8$samples$population[4, ])
plot(fm8$samples$population[5, ])
plot(fm8$samples$population[6, ])
plot(fm8$samples$population[7, ])
plot(fm8$samples$population[8, ])
plot(fm8$samples$population[9, ])
plot(fm8$samples$population[10, ])
plot(fm8$samples$population[11, ])

plot(fm8$samples$subjects[1, 1, ])
plot(fm8$samples$subjects[2, 1, ])
plot(fm8$samples$subjects[3, 1, ])
plot(fm8$samples$subjects[4, 1, ])
plot(fm8$samples$subjects[5, 1, ])
plot(fm8$samples$subjects[6, 1, ])
plot(fm8$samples$subjects[7, 1, ])
plot(fm8$samples$subjects[8, 1, ])
plot(fm8$samples$subjects[9, 1, ])
plot(fm8$samples$subjects[10, 1, ])
plot(fm8$samples$subjects[11, 1, ])

plot(1 / fm8$samples$precision$pop)
plot(1 / fm8$samples$precision$sub)
plot(fm8$samples$precision$eps)

source("graphs.R")
plot_spline(fm8, mle = T)
plot_spline(truncate_spline(fm8, 5000))


set.seed(2)
source("main-tpf.R")
system.time(fm9 <- sub_tpf(growth10, 8, deg = 2, penalty = FALSE, shape = "increasing", size = 10000, burn = 0, verbose = T))

plot(fm9$samples$population[1, ])
plot(fm9$samples$population[2, ])
plot(fm9$samples$population[3, ])
plot(fm9$samples$population[4, ])
plot(fm9$samples$population[5, ])
plot(fm9$samples$population[6, ])
plot(fm9$samples$population[7, ])
plot(fm9$samples$population[8, ])
plot(fm9$samples$population[9, ])
plot(fm9$samples$population[10, ])
plot(fm9$samples$population[11, ])

hist(fm9$samples$subjects[1, 10, ])
hist(fm9$samples$subjects[2, 10, ])
hist(fm9$samples$subjects[3, 10, ])
hist(fm9$samples$subjects[4, 10, ])
hist(fm9$samples$subjects[5, 10, ])
hist(fm9$samples$subjects[6, 10, ])
hist(fm9$samples$subjects[7, 10, ])
hist(fm9$samples$subjects[8, 10, ])
hist(fm9$samples$subjects[9, 10, ])
hist(fm9$samples$subjects[10, 10, ])
hist(fm9$samples$subjects[11, 10, ])

plot(1 / fm9$samples$precision$pop)
plot(1 / fm9$samples$precision$sub)
plot(fm9$samples$precision$eps)

plot_spline(fm9)
plot_spline(truncate_spline(fm9, 10000))

set.seed(1)
source("main-tpf.R")
system.time(fm10 <- SubjectsTpf(growth, 8, deg = 2, shape = "increasing", size = 20000, burn = 0, verbose = T))

plot(fm10$samples$population[1, ])
plot(fm10$samples$population[2, ])
plot(fm10$samples$population[3, ])
plot(fm10$samples$population[4, ])
plot(fm10$samples$population[5, ])
plot(fm10$samples$population[6, ])
plot(fm10$samples$population[7, ])
plot(fm10$samples$population[8, ])


plot(fm10$samples$subjects[1, 1, ])
plot(fm10$samples$subjects[2, 1, ])
plot(fm10$samples$subjects[3, 1, ])
plot(fm10$samples$subjects[4, 1, ])
plot(fm10$samples$subjects[5, 1, ])
plot(fm10$samples$subjects[6, 1, ])
plot(fm10$samples$subjects[7, 1, ])
plot(fm10$samples$subjects[8, 1, ])

plot(1 / fm10$samples$precision$pop)
plot(1 / fm10$samples$precision$sub)
plot(fm10$samples$precision$eps)

plot_spline(fm10)
plot_spline(truncate_spline(fm10, 10000))

system.time(fm11 <- SubjectsTpfMul(growth, 8, deg = 2, shape = "increasing", size = 10000, burn = 0))
plot(fm11$samples$population$hgtf[1, ])
plot(fm11$samples$population$hgtf[2, ])
plot(fm11$samples$population$hgtf[3, ])
plot(fm11$samples$population$hgtf[4, ])
plot(fm11$samples$population$hgtf[5, ])
plot(fm11$samples$population$hgtf[6, ])
plot(fm11$samples$population$hgtf[7, ])
plot(fm11$samples$population$hgtf[8, ])


plot(fm11$samples$subjects$hgtf[1, 1, ])
plot(fm11$samples$subjects$hgtf[2, 1, ])
plot(fm11$samples$subjects$hgtf[3, 1, ])
plot(fm11$samples$subjects$hgtf[4, 1, ])
plot(fm11$samples$subjects$hgtf[5, 1, ])
plot(fm11$samples$subjects$hgtf[6, 1, ])
plot(fm11$samples$subjects$hgtf[7, 1, ])
plot(fm11$samples$subjects$hgtf[8, 1, ])


plot_spline(fm11)
fm11t <- truncate_spline(fm11, 1500)
plot(1 / fm10t$samples$precision$pop)
plot(fm10t$samples$subject$hgtm[3, "boy01", ])

plot_spline(fm11t)
get_array_invs <- function(ary) {
    array(apply(ary, 3, solve), dim(ary))
}

## models without random effects on spline coefficients
source("pop-tpf.R")
set.seed(1)
fmpop <- pop_tpf(growth10, 8, deg = 2, random = 11, shape = "increasing", size = 100, burn = 0, verbose = T)
plot_spline(fmpop)


## FAKE data
source("main-tpf.R")
design_ls <- get_design_tpf(growth10$x, 5, 2)

fake_data <- data.frame()
