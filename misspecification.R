rm(list = ls())

library(devtools)
install_github("maryclare/LANOVA")
library(LANOVA)
library(leapp)
library(parallel)

set.seed(1)

n <- 25 # 25
p <- 25 # 25

H.1 <- diag(n) - rep(1, n)%*%t(rep(1, n))/n
H.2 <- diag(p) - rep(1, p)%*%t(rep(1, p))/p

s2e <- 1
k <- 2
s2cs <- c(1/k, 1, (2*k)/k)*s2e

# prs <- seq(0.1, 0.9, by = 0.1)
prs <- seq(0, 1, by = 0.1)

sims <- 500 # Should be 500

# WARNING: This takes a very long time to run if 500 simulations are performed
sim.prs <- mclapply(prs, function(pr) {

  mse.mle <- array(NA, dim = c(length(s2cs), sims))
  mse.sur <- array(NA, dim = c(length(s2cs), sims))
  mse.uni <- array(NA, dim = c(length(s2cs), sims))
  mse.add <- array(NA, dim = c(length(s2cs), sims))
  mse.lan <- array(NA, dim = c(length(s2cs), sims))
  mse.low.1 <- array(NA, dim = c(length(s2cs), sims))
  mse.low.5 <- array(NA, dim = c(length(s2cs), sims))
  mse.ipod.s <- array(NA, dim = c(length(s2cs), sims))
  mse.ipod.h <- array(NA, dim = c(length(s2cs), sims))
  
  allz.lan <- array(NA, dim = c(length(s2cs), sims))
  sig.sq.zs <- array(NA, dim = c(length(s2cs), sims))
  sig.sq.ma <- array(NA, dim = c(length(s2cs), sims))
  lambda.cs <- array(NA, dim = c(length(s2cs), sims))
  lanova.times <- ipod.h.times <- ipod.s.times <- array(NA, dim = c(length(s2cs), sims))

  for (s2c in s2cs) {
    cat("s2c=", s2c, "\n")

    for (i in 1:sims) {
      cat("i=", i, "\n")
      m <- 0
      a<- rep(0, n)
      b<- rep(0, p)
      # C<-matrix( rnorm(n*p, mean = 0, sd = sqrt(s2c/(pr)))*rbinom(n*p,1,pr), nrow = n, ncol = p)
      C<-matrix( rnorm(n*p, mean = 0, sd = sqrt(s2c))*rbinom(n*p,1,pr), nrow = n, ncol = p)
      
      Theta<- m + outer(a,b,"+") + C

      Y<-Theta + matrix(rnorm(n*p, mean = 0, sd = sqrt(s2e)),n,p)

      lanova.start <- proc.time()
      ett <- estTuneTest(Y)
      
      # Modify for estimation
      if (ett$sigma.2.z <= 10^(-4)) {
        Theta.f <- Y
      } else {
        
        fit <- estM(Y, ett)
        
        Theta.f <- fit$mu.f + outer(fit$a.f, fit$b.f, "+") + fit$C.f
      }
      lanova.end <- proc.time()
    
      Z <- crossprod(t(H.1), crossprod(t(Y), H.2))
        
      # Comparison to IPOD
      ipod.s.start <- proc.time()

      X <- cbind(rep(1, prod(dim(Y))), (rep(1, p)%x%diag(n))[, -1], (diag(p)%x%rep(1, n))[, -1])
      H <- crossprod(t(X), tcrossprod(solve(crossprod(X)), X))

      ipod.s <- IPOD(X = X, Y = c(Y), H = H, method = "soft")
      c.ipod.s <- ipod.s$gamma
      Theta.ipod.s <- matrix(lm(I(c(Y) - c.ipod.s)~X-1)$fitted.values + c.ipod.s, nrow = n, ncol = p)

      ipod.s.end <- proc.time()

      ipod.h.start <- proc.time()

      ipod.h <- IPOD(X = X, Y = c(Y), H = H, method = "hard")
      c.ipod.h <- ipod.h$gamma
      Theta.ipod.h <- matrix(lm(I(c(Y) - c.ipod.h)~X-1)$fitted.values + c.ipod.s, nrow = n, ncol = p)

      ipod.h.end <- proc.time()

      ipod.s.time <- (ipod.s.end - ipod.s.start)
      ipod.h.time <- (ipod.h.end - ipod.h.start)

      lanova.time <- (lanova.end - lanova.start)
      
      Theta.add <- Y - Z
      Z.svd <- svd(Z)
  
      Theta.low.1 <- Theta.add +  Z.svd$u%*%diag(ifelse(1:length(Z.svd$d) <= 1, Z.svd$d, 0))%*%t(Z.svd$v)
      Theta.low.5 <- Theta.add +  Z.svd$u%*%diag(ifelse(1:length(Z.svd$d) <= 5, Z.svd$d, 0))%*%t(Z.svd$v)

      sig.sq.mad <- (median(c(abs(Z)))/0.6745)^2
      
      lambda.c.uni <- sqrt(sig.sq.mad*2*log(n*p))
      
      Theta.uni <- (Y - Z) + softt(Z, lambda.c.uni*1*2)
      
      lams.sure <- seq(0, lambda.c.uni, length.out = 10000)
      sure.risk <- unlist(lapply(lams.sure, function(x) {
        mean(c(sig.sq.mad + pmin(Z^2, x^2) - 2*sig.sq.mad*(Z^2 <= x^2)))
      }))
      lam.sure.min <- lams.sure[sure.risk == min(sure.risk)][1]
      
      Theta.sur <- (Y - Z) + softt(Z, lam.sure.min*1*2)
      
      mse.mle[which(s2c == s2cs), i] <- mean(as.vector(Y - Theta)^2)
      mse.lan[which(s2c == s2cs), i] <- mean(as.vector(Theta.f - Theta)^2)
      mse.uni[which(s2c == s2cs), i] <- mean(as.vector(Theta.uni - Theta)^2)
      mse.sur[which(s2c == s2cs), i] <- mean(as.vector(Theta.sur - Theta)^2)
      mse.low.1[which(s2c == s2cs), i] <- mean(as.vector(Theta.low.1 - Theta)^2)
      mse.low.5[which(s2c == s2cs), i] <- mean(as.vector(Theta.low.5 - Theta)^2)
      mse.ipod.h[which(s2c == s2cs), i] <- mean(as.vector(Theta.ipod.h - Theta)^2)
      mse.ipod.s[which(s2c == s2cs), i] <- mean(as.vector(Theta.ipod.s - Theta)^2)
      mse.add[which(s2c == s2cs), i] <- mean(as.vector(Theta.add - Theta)^2)
      
      
      ipod.s.times[which(s2c == s2cs), i] <- ipod.s.time[3]
      ipod.h.times[which(s2c == s2cs), i] <- ipod.h.time[3]
      sig.sq.ma[which(s2c == s2cs), i] <- sig.sq.mad
      sig.sq.zs[which(s2c == s2cs), i] <- ett$sigma.2.z
      lambda.cs[which(s2c == s2cs), i] <- ett$lambda.c
      lanova.times[which(s2c == s2cs), i] <- lanova.time[3]
    }
  }

  return(list("mse.mle" = mse.mle, 
              "mse.sur" = mse.sur,
              "mse.uni" = mse.uni, 
              "mse.lan" = mse.lan, "mse.low.1" = mse.low.1, "mse.low.5" = mse.low.5,
              "mse.ipod.h" = mse.ipod.h, "mse.ipod.s" = mse.ipod.s, 
              "mse.add" = mse.add, 
              "sig.sq.zs" = sig.sq.zs, 
              "sig.sq.ma" = sig.sq.ma, "lambda.cs" = lambda.cs,
              "ipod.s.time" = ipod.s.times,
              "ipod.h.time" = ipod.h.times,
              "lanova.time" = lanova.times))

}, mc.cores = length(prs))

save.image(paste("misspecification_", k, ".RData", sep = ""))
