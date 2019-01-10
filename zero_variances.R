rm(list = ls())

library(devtools)
install_github("maryclare/LANOVA")
install_github("maryclare/gnorm")
library(LANOVA)
library(RColorBrewer)
library(gnorm)
# Get colors for later
set.seed(1)

n <- 25 # 25
p <- 25 # 25

s2cs <- c(1/2, 1, 1.5)
s2e <- 1

sims <- 10000

sig.sq.zs <- array(NA, dim = c(length(s2cs), sims))
sig.sq.cs <- array(NA, dim = c(length(s2cs), sims))

for (s2c in s2cs) {
    for (i in 1:sims) {
      m <- 0
      a<- rep(0, n)
      b<- rep(0, p)
      C<-matrix( sqrt(s2c)*rgnorm(n*p, 0, 1, 1)/sqrt(2), nrow = n, ncol = p)
      
      Theta<- m + outer(a,b,"+") + C
      
      Y<-Theta + matrix(rnorm(n*p, mean = 0, sd = sqrt(s2e)),n,p)
      
      y.bar.2 <- mean(Y^2)
      y.bar.4 <- mean(Y^4)
      
      es4c <- y.bar.4/3 - y.bar.2^2
      es2c <- sqrt(ifelse(es4c > 0, es4c, 0))
      es2z <- y.bar.2 - es2c
      
      sig.sq.zs[which(s2c == s2cs), i] <- es2z
      sig.sq.cs[which(s2c == s2cs), i] <- es2c
    }
}

apply(sig.sq.zs, 1, function(x) {mean(x == 0)})
apply(sig.sq.cs, 1, function(x) {mean(x == 0)})

