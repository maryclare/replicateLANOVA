rm(list = ls())

library(devtools)
install_github("maryclare/LANOVA")
library(LANOVA)
library(RColorBrewer)
library(xtable)
library(leapp)
library(parallel)
library(wavelets)
library(gnorm)

# Get colors for later
cols <- brewer.pal(11, "RdBu")
cols[6] <- cols[3]
cols <- cols[c(2, 6, length(cols) - 1)]
set.seed(1)

n <- 25 # 25
p <- 25 # 25

H.1 <- diag(n) - rep(1, n)%*%t(rep(1, n))/n
H.2 <- diag(p) - rep(1, p)%*%t(rep(1, p))/p

s2cs <- c(1/2, 1, 1.5)
s2e <- 1

sims <- 500


mse.mle <- array(NA, dim = c(length(s2cs), sims))
mse.add <- array(NA, dim = c(length(s2cs), sims))
mse.lan <- array(NA, dim = c(length(s2cs), sims))
mse.low.1 <- array(NA, dim = c(length(s2cs), sims))
mse.low.5 <- array(NA, dim = c(length(s2cs), sims))
mse.ipod.s <- array(NA, dim = c(length(s2cs), sims))
mse.ipod.h <- array(NA, dim = c(length(s2cs), sims))
allz.lan <- array(NA, dim = c(length(s2cs), sims))
sig.sq.zs <- array(NA, dim = c(length(s2cs), sims))
sig.sq.zs.mad <- array(NA, dim = c(length(s2cs), sims))
lambda.cs <- array(NA, dim = c(length(s2cs), sims))

for (s2c in s2cs) {
  
  for (i in 1:sims) {  
    m <- 0 
    a<- rep(0, n)
    b<- rep(0, p)
    C<-sqrt(s2c)*matrix( rgnorm(n*p, mu = 0, alpha = sqrt(gamma(1)/gamma(3)), beta = 1), nrow = n, ncol = p)
    
    Theta<- m + outer(a,b,"+") + C
    
    Y<-Theta + matrix(rnorm(n*p, mean = 0, sd = sqrt(s2e)),n,p)
    
    ett <- estTuneTest(Y)
    # 
    # # Modify for estimation
    # if (ett$sigma.2.z <= 10^(-4)) {
    #   Theta.f <- Y
    # } else {
    #   
    #   fit <- estM(Y, ett)
    #   
    #   Theta.f <- fit$mu.f + outer(fit$a.f, fit$b.f, "+") + fit$C.f
    # }
    R <- H.1%*%Y%*%H.2
    
    sig.sq.zs[which(s2c == s2cs), i] <- ett$sigma.2.z
    sig.sq.zs.mad[which(s2c == s2cs), i] <- (median(c(abs(R)))/0.6745)^2
    # lambda.cs[which(s2c == s2cs), i] <- ett$lambda.c
  }
}

plot((sig.sq.zs[1, ]), (sig.sq.zs.mad[1, ]))

plot((sig.sq.zs[1, ] - 1)^2, (sig.sq.zs.mad[1, ] - 1)^2)
abline(a = 0, b = 1)
points((sig.sq.zs[2, ] - 1)^2, (sig.sq.zs.mad[2, ] - 1)^2, col = "red")
points((sig.sq.zs[3, ] - 1)^2, (sig.sq.zs.mad[3, ] - 1)^2, col = "blue")

apply((sig.sq.zs.mad - 1)^2, 1, mean)
apply((sig.sq.zs - 1)^2, 1, mean)

apply((sig.sq.zs.mad - 1)^2 > (sig.sq.zs - 1)^2, 1, mean)
      