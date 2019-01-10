# Test method on eeukijk kroonenberg data
rm(list = ls())

library(devtools)
install_github("maryclare/LANOVA")
library(LANOVA)
library(RColorBrewer)
library(leapp)

data(fusarium)
Y <- fusarium[["Y"]]

# Obtain tuning parameter estimates
lanova.start <- proc.time()
ett <- estTuneTest(Y, lowpen = FALSE)

# Reject null of no interactions
ett$test.stat > qnorm(0.95)
# P-value: want to compute Pr(z > test.stat)
exp(pnorm(ett$test.stat, log = TRUE, lower.tail = FALSE))

fit <- estM(Y, ett)
lanova.stop <- proc.time()
lanova.time <- lanova.stop - lanova.start


# Prep some useful objects for later
n <- dim(Y)[1]
p <- dim(Y)[2]
q <- dim(Y)[3]

H.1 <- diag(n) - rep(1, n)%*%t(rep(1, n))/n
H.2 <- diag(p) - rep(1, p)%*%t(rep(1, p))/p
H.3 <- diag(q) - rep(1, q)%*%t(rep(1, q))/q



C.f <- fit$C.f
M.f <-fit$M.f

ipod.s.start <- proc.time()
inds.n <- array(dim = dim(Y))
for (i in 1:n) {
  inds.n[i, , ] <- i
}
inds.p <- array(dim = dim(Y))
for (i in 1:p) {
  inds.p[, i, ] <- i
}
inds.q <- array(dim = dim(Y))
for (i in 1:q) {
  inds.q[, , i] <- i
}


X <- model.matrix(~factor(c(inds.n))+factor(c(inds.p))+factor(c(inds.q))+factor(c(inds.n)):factor(c(inds.p))+
                    factor(c(inds.n)):factor(c(inds.q))+factor(c(inds.n)):factor(c(inds.q)))
H <- X%*%solve(crossprod(X))%*%t(X)

ipod.s <- IPOD(X = X, Y = c(Y), H = H, method = "soft")
C.ipod.s <- array(ipod.s$gamma, dim = dim(Y))
Theta.ipod.s <- array(lm(I(c(Y) - c(C.ipod.s))~X-1)$fitted.values, dim = dim(Y)) + C.ipod.s
ipod.s.stop <- proc.time()

ipod.s.time <- ipod.s.stop - ipod.s.start

ipod.h.start <- proc.time()
inds.n <- array(dim = dim(Y))
for (i in 1:n) {
  inds.n[i, , ] <- i
}
inds.p <- array(dim = dim(Y))
for (i in 1:p) {
  inds.p[, i, ] <- i
}
inds.q <- array(dim = dim(Y))
for (i in 1:q) {
  inds.q[, , i] <- i
}


X <- model.matrix(~factor(c(inds.n))+factor(c(inds.p))+factor(c(inds.q))+factor(c(inds.n)):factor(c(inds.p))+
                    factor(c(inds.n)):factor(c(inds.q))+factor(c(inds.n)):factor(c(inds.q)))
H <- X%*%solve(crossprod(X))%*%t(X)

ipod.h <- IPOD(X = X, Y = c(Y), H = H, method = "hard")
C.ipod.h <- array(ipod.h$gamma, dim = dim(Y))
Theta.ipod.h <- array(lm(I(c(Y) - c(C.ipod.h))~X-1)$fitted.values, dim = dim(Y)) + C.ipod.h
ipod.h.stop <- proc.time()
ipod.h.time <- ipod.h.stop - ipod.h.start

(ipod.s.time/lanova.time)[3]
(ipod.h.time/lanova.time)[3]

# Attempt marginal maximum likelihood
lam <- ett$lambda.c
sig.sq <- ett$sigma.2.z
y.tilde <- (diag(nrow(H)) - H)%*%c(Y)
lambdas <- numeric(1000)
sig.sqs <- numeric(length(lambdas))
i <- 1
diff <- Inf
while (i == 1 | diff > 10^(-7)) {
  samples <- bayes.lasso(y = y.tilde, X = diag(nrow(H)), lambda = lam,
                             sigma.sq = sig.sq, num.samp = 100000, delta = 1)
  rss <- apply(samples$samples.beta, 1, function(x) {sum((y.tilde - x))^2})
  sig.sq <- mean(rss)/(length(y.tilde) - ncol(X))
  lam <- length(y.tilde)/(mean(rowSums(abs(samples$samples.beta))))
  sig.sqs[i] <- sig.sq
  lambdas[i] <- lam
  par(mfrow = c(1, 2))
  plot(lambdas[1:i])
  abline(h = ett$lambda.c, lty = 2)
  plot(sig.sqs[1:i])
  abline(h = ett$sigma.2.z, lty = 2)
  if (i > 1) {
    diff <- mean(c(abs(lambdas[i] - lambdas[i - 1]), abs(sig.sqs[i] - sig.sqs[i - 1])))
    cat("Diff=", diff, "\n")
  }
  i <- i + 1
  
}

# Soft thresholding sets them *all* to zero
plot(c(C.f), c(C.ipod.s))

plot(c(C.f), c(C.ipod.h))
abline(a = 0, b = 1)
sum(C.f != 0)
mean(C.f != 0)

sum(C.ipod.h != 0)
mean(C.ipod.h != 0)

table("LANOVA" = c(C.f == 0), "IPOD H" = c(C.ipod.h == 0))

cols <- brewer.pal(11, "RdBu")
cols[6] <- "white"
breaks <- c(-1*c(10^(-7), quantile(c(abs(C.f)[C.f != 0], abs(C.ipod.h)[C.ipod.h != 0]), c(1:5/5)))[6:1],
            c(10^(-7), quantile(c(abs(C.f)[C.f != 0], abs(C.ipod.h)[C.ipod.h != 0]), c(1:5/5))))
cols.C.f <- apply(C.f, c(1, 2, 3), function(x) {
  ifelse(x <= breaks[2], cols[1],
         ifelse(x <= breaks[3], cols[2],
                ifelse(x <= breaks[4], cols[3],
                       ifelse(x <= breaks[5], cols[4],
                              ifelse(x <= breaks[6], cols[5],
                                     ifelse(x <= breaks[7], cols[6],
                                            ifelse(x <= breaks[8], cols[7],
                                                   ifelse(x <= breaks[9], cols[8],
                                                          ifelse(x <= breaks[10], cols[9],
                                                                 ifelse(x <= breaks[11], cols[10], cols[11]))))))))))
  
})
cols.C.ipod.h <- apply(C.ipod.h, c(1, 2, 3), function(x) {
  ifelse(x <= breaks[2], cols[1],
         ifelse(x <= breaks[3], cols[2],
                ifelse(x <= breaks[4], cols[3],
                       ifelse(x <= breaks[5], cols[4],
                              ifelse(x <= breaks[6], cols[5],
                                     ifelse(x <= breaks[7], cols[6],
                                            ifelse(x <= breaks[8], cols[7],
                                                   ifelse(x <= breaks[9], cols[8],
                                                          ifelse(x <= breaks[10], cols[9],
                                                                 ifelse(x <= breaks[11], cols[10], cols[11]))))))))))
  
})

years <- 1990:1993

pdf("~/Dropbox/EB_Fusion/LassoANOVA/Draft/CSDASubmission_September/fusarium.pdf", family = "Times",
    width = 6, height = 5)
par(mfrow = c(2, 4))
par(cex.axis = 1)
par(cex.lab = 1)
par(cex.main = 1)
par(mar = c(3, 3, 1, 0.1))
par(mgp = c(1.75, 0.75, 0))
par(oma = c(0.5, 0.5, 1, 0))
cex.top <- 0.9
for (i in 1:n) {
  plot(seq(1, p, length.out = 10), seq(1, q, length.out = 10), type = "n", xlab = "", ylab = "",
       main = years[i], axes = FALSE)
  mtext("LANOVA Estimates", 3, outer = TRUE, cex = cex.top)
  # axis(1, 0:(p + 1), c("", c("7F", "12H", "13H", "14H", "15H", "16G", "17G"), ""), las = 2)
  # mtext("Blight Strains", 1, line = 2.5, cex = 0.75)
  axis(1, 0:(p + 1), c("", rep("", p), ""), las = 2)
  if (i == 1) {
    axis(2, 0:(q + 1), c("", c(paste("N", 1:5, sep = ""),
                               paste("F", 11:15, sep = ""),
                               paste("G", 16:20, sep = ""),
                               paste("H", 21:25, sep = "")), ""), las = 1)
    mtext("Wheat Varieties", 2, line = 2.5, cex = 0.75)
  } else {axis(2, c(0, q + 1), rep("", 2))}
  for (j in 1:p) {
  points(rep(j, q), 1:q, col = cols.C.f[i, j, ],
         pch = 15, cex = 2)
  }
  abline(h = 5.5, lty = 3)
  abline(h = 10.5, lty = 3)
  abline(h = 15.5, lty = 3)
  abline(v = 1.5, lty = 3)
  abline(v = 5.5, lty = 3)
}
mtext("IPOD Hard Thresholding Estimates", 3, outer = TRUE, line = -18,  cex = cex.top)
for (i in 1:n) {
  plot(seq(1, p, length.out = 10), seq(1, q, length.out = 10), type = "n", xlab = "", ylab = "",
       main = years[i], axes = FALSE)
  axis(1, 0:(p + 1), c("", c("7F", "12H", "13H", "14H", "15H", "16G", "17G"), ""), las = 2)
  mtext("Blight Strains", 1, line = 2.5, cex = 0.75)
  if (i == 1) {
    axis(2, 0:(q + 1), c("", c(paste("N", 1:5, sep = ""),
                               paste("F", 11:15, sep = ""),
                               paste("G", 16:20, sep = ""),
                               paste("H", 21:25, sep = "")), ""), las = 1)
    mtext("Wheat Varieties", 2, line = 2.5, cex = 0.75)
  } else {axis(2, c(0, q + 1), rep("", 2))}
  for (j in 1:p) {
    points(rep(j, q), 1:q, col = cols.C.ipod.h[i, j, ],
           pch = 15, cex = 2)
  }
  abline(h = 5.5, lty = 3)
  abline(h = 10.5, lty = 3)
  abline(h = 15.5, lty = 3)
  abline(v = 1.5, lty = 3)
  abline(v = 5.5, lty = 3)
}
dev.off()

par(mfrow = c(1, 1))
plot(c((H.3%x%H.2%x%H.1)%*%c(Y)), c((H.3%x%H.2%x%H.1)%*%c(M.f)), xlab = "Residuals", ylab = "Estimates of C")
abline(a = 0, b = 1, lty = 2)
plot(c((H.3%x%H.2%x%H.1)%*%c(Y)), c(C.f), xlab = "Residuals", ylab = "Estimates of C")

hist(c((H.3%x%H.2%x%H.1)%*%c(Y)))
hist(c((H.3%x%H.2%x%H.1)%*%c(M.f)), add = TRUE, col = "red")

blights <- c("7F", "12H", "13H", "14H", "15H", "16G", "17G")
wheats <- c(paste("N", 1:5, sep = ""),
            paste("F", 11:15, sep = ""),
            paste("G", 16:20, sep = ""),
            paste("H", 21:25, sep = ""))

zero.f <- cbind(years[which(C.ipod.h != 0 & C.f == 0, arr.ind = TRUE)[, 1]],
      blights[which(C.ipod.h != 0 & C.f == 0, arr.ind = TRUE)[, 2]],
      wheats[which(C.ipod.h != 0 & C.f == 0, arr.ind = TRUE)[, 3]])

zero.h <- cbind(years[which(C.ipod.h == 0 & C.f != 0, arr.ind = TRUE)[, 1]],
      blights[which(C.ipod.h == 0 & C.f != 0, arr.ind = TRUE)[, 2]],
      wheats[which(C.ipod.h == 0 & C.f != 0, arr.ind = TRUE)[, 3]])

C.f[which(years == 1992), , which(wheats %in% c("H21", "H23"))]
C.ipod.h[which(years == 1992), , which(wheats %in% c("H21", "H23"))]
wheats
