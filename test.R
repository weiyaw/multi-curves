sitka <- read.table("data/tsitka.txt", header=T)
lines(log.size ~ days, sitka, col = id.num)

x0 <- c(3,1)
d <- 3
sig <- matrix(c(3,2,2,3), 2)
eValues <- eigen(solve(sig))$values
eVectors <- eigen(solve(sig))$vectors
len <- d / sqrt(eValues)
plot(c(x0[1] + max(len), x0[1] - max(len)) , c(x0[2] + max(len), x0[2] - max(len)))
ellipse(x0, sig, d, col = 'red')
ending1 <- x0 + eVectors[,1] * len[1]
ending2 <- x0 + eVectors[,2] * len[2]
arrows(x0[1], x0[2], ending1[1], ending1[2])
arrows(x0[1], x0[2], ending2[1], ending2[2])

x0 <- as.vector(t(eVectors) %*% c(3, 1))
d <- 3
sig <- diag(eValues)
eValuesP <- eigen(solve(sig))$values
eVectorsP <- eigen(solve(sig))$vectors
len <- d / sqrt(eValuesP)
plot(c(x0[1] + max(len), x0[1] - max(len)) , c(x0[2] + max(len), x0[2] - max(len)))
ellipse(x0, sig, d, col = 'red')
ending1 <- x0 + eVectorsP[,1] * len[1]
ending2 <- x0 + eVectorsP[,2] * len[2]
arrows(x0[1], x0[2], ending1[1], ending1[2])
arrows(x0[1], x0[2], ending2[1], ending2[2])

