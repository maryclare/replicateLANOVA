### Compare Information
rm(list = ls())

library(Deriv)
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
drule[["erfc"]] <- alist(x=-2*exp(-x^2)/sqrt(pi))
nl.lik <- function(z, s2e, s2c) {
  lam <- sqrt(2/s2c)
  (lam/4)*exp(s2e*lam^2/2)*(exp(-z*lam)*erfc(-z/sqrt(2*s2e) + sqrt(s2e)*lam/sqrt(2)) + 
                              exp(z*lam)*erfc(z/sqrt(2*s2e) + sqrt(s2e)*lam/sqrt(2)))
}
dl.s2e <- function(z, s2e, s2c) {
  # Deriv(log(nl.lik(z, s2e, s2c)), c("s2e"))
  .e2 <- sqrt(2/s2c)
  .e4 <- sqrt(2 * s2e)
  .e5 <- sqrt(s2e)
  .e7 <- .e2 * .e5/1.4142135623731
  .e8 <- z * .e2
  .e9 <- z/.e4
  .e11 <- 0.353553390593274 * (.e2/.e5)
  .e13 <- exp(-.e8)
  .e14 <- exp(.e8)
  .e15 <- .e7 - .e9
  .e16 <- .e7 + .e9
  .e17 <- z/(2 * (s2e * .e4))
  0.25 * (((erfc(.e15) * .e13 + erfc(.e16) * .e14)/s2c - (2 * 
                                                            ((.e11 - .e17) * exp(-.e16^2) * .e14) + 2 * ((.e11 + 
                                                                                                            .e17) * exp(-.e15^2) * .e13))/sqrt(pi)) * exp(s2e/s2c) * 
            .e2/nl.lik(z, s2e, s2c))
}
dl.s2c <- function(z, s2e, s2c) {
  # Deriv(log(nl.lik(z, s2e, s2c)), c("s2c"))
  .e2 <- sqrt(2/s2c)
  .e3 <- sqrt(s2e)
  .e7 <- .e2 * .e3/1.4142135623731
  .e8 <- z/sqrt(2 * s2e)
  .e9 <- z * .e2
  .e10 <- .e7 - .e8
  .e11 <- .e7 + .e8
  .e13 <- erfc(.e10)
  .e14 <- erfc(.e11)
  .e15 <- exp(-.e9)
  .e16 <- exp(.e9)
  .e17 <- .e13 * .e15
  .e18 <- .e14 * .e16
  .e19 <- sqrt(pi)
  -(0.25 * (((.e16 * (z * .e14 - 1.41421356237309 * (exp(-.e11^2) * 
                                                       .e3/.e19)) + s2e * (.e17 + .e18) * .e2 - (1.41421356237309 * 
                                                                                                   (exp(-.e10^2) * .e3/.e19) + z * .e13) * .e15) * .e2 + 
               .e17 + .e18) * exp(s2e/s2c)/(s2c^2 * nl.lik(z, s2e, s2c) * 
                                              .e2)))
  }
iinv.11 <- function(z, s2e, s2c) {
  dl.s2e(z, s2e, s2c)^2*nl.lik(z, s2e, s2c)
}
iinv.12 <- function(z, s2e, s2c) {
  dl.s2e(z, s2e, s2c)*dl.s2c(z, s2e, s2c)*nl.lik(z, s2e, s2c)
}
iinv.22 <- function(z, s2e, s2c) {
  dl.s2c(z, s2e, s2c)^2*nl.lik(z, s2e, s2c)
}
# Compute variance of s2c moment estimator
convert.s2c <- function(x, y) {
  ((y/3 - x^2))^{1/2}
}  
convgrad.s2c <- Deriv(convert.s2c, c("x", "y"))

convert.s2e <- function(x, y) {
  x - ((y/3 - x^2))^{1/2}
}  
convgrad.s2e <- Deriv(convert.s2e, c("x", "y"))

momvar.s2c.non <- function(s2e, s2c) {
  a <- (5*s2c^2 + 4*s2c*s2e + 2*s2e^2)
  b <- c <- 6*(14*s2c^3 + 13*s2c^2*s2e + 6*s2c*s2e^2 + 2*s2e^3)
  d <- 12*(207*s2c^4 + 204*s2c^3*s2e + 99*s2c^2*s2e^2 + 32*s2c*s2e^3 + 8*s2e^4)
  vmat <- rbind(c(a, b), c(c, d))
  estgrad.s2c <- convgrad.s2c(s2c + s2e, 6*s2c^2 + 6*s2c*s2e + 3*s2e^2)
  return(t(estgrad.s2c)%*%vmat%*%estgrad.s2c)
}

momvar.s2e.non <- function(s2e, s2c) {
  a <- (5*s2c^2 + 4*s2c*s2e + 2*s2e^2)
  b <- c <- 6*(14*s2c^3 + 13*s2c^2*s2e + 6*s2c*s2e^2 + 2*s2e^3)
  d <- 12*(207*s2c^4 + 204*s2c^3*s2e + 99*s2c^2*s2e^2 + 32*s2c*s2e^3 + 8*s2e^4)
  vmat <- rbind(c(a, b), c(c, d))
  estgrad.s2e <- convgrad.s2e(s2c + s2e, 6*s2c^2 + 6*s2c*s2e + 3*s2e^2)
  return(t(estgrad.s2e)%*%vmat%*%estgrad.s2e)
}


vals.s2e <- exp(seq(log(0.01), log(1), length.out = 25))
vals.s2c <- exp(seq(log(0.01), log(1), length.out = length(vals.s2e)))
momvar.s2c <- momvar.s2e <- matrix(NA, nrow = length(vals.s2c), ncol = length(vals.s2e))
info.s2c <- info.s2e <- matrix(NA, nrow = length(vals.s2c), ncol = length(vals.s2e))
for (i in 1:length(vals.s2c)) {
  s2c.i <- vals.s2c[i]
  lc.i <- sqrt(2/s2c.i)
  for (j in 1:length(vals.s2e)) {
    
    # Compute variance of s2c and s2e estimates using vdv theorem 4.1
    a <- (5*s2c.i^2 + 4*s2c.i*vals.s2e[j] + 2*vals.s2e[j]^2)
    b <- c <- 6*(14*s2c.i^3 + 13*s2c.i^2*vals.s2e[j] + 6*s2c.i*vals.s2e[j]^2 + 2*vals.s2e[j]^3)
    d <- 12*(207*s2c.i^4 + 204*s2c.i^3*vals.s2e[j] + 99*s2c.i^2*vals.s2e[j]^2 + 32*s2c.i*vals.s2e[j]^3 + 8*vals.s2e[j]^4)
    vmat <- rbind(c(a, b), c(c, d))
    
    e.theta.p <- rbind(c(1, 1), c(12*s2c.i + 6*vals.s2e[j], 6*(vals.s2e[j] + s2c.i)))
    parvar <- solve(e.theta.p)%*%vmat%*%t(solve(e.theta.p))
    
    momvar.s2c[i, j] <- diag(parvar)[1] # momvar.s2c.non(s2e = vals.s2e[j], s2c = s2c.i), gives same solution
    momvar.s2e[i, j] <- diag(parvar)[2] # momvar.s2e.non(s2e = vals.s2e[j], s2c = s2c.i), gives same solution
    a <- integrate(iinv.11, lower = -500/lc.i, upper = 500/lc.i, 
                   s2e = vals.s2e[j], s2c = s2c.i)$value
    b <- c <- integrate(iinv.12, lower = -500/lc.i, upper = 500/lc.i, 
                        s2e = vals.s2e[j], s2c = s2c.i)$value
    d <- integrate(iinv.22, lower = -500/lc.i, upper = 500/lc.i, 
                   s2e = vals.s2e[j], s2c = s2c.i)$value
    ii <- (solve(rbind(c(a, b), c(c, d))))
    info.s2c[i, j] <- ii[2, 2]
    info.s2e[i, j] <- ii[1, 1]
  }
}

pdf("comptomle.pdf", family = "Times")
par(mfrow = c(1, 1))
# par(oma = c(0, 0, 0, 0))
par(mar = c(5, 4, 4, 2) + 1)
par(cex.axis = 1.5)
par(cex.lab = 1.5)
par(cex.main = 1.5)
contour(vals.s2c, vals.s2e, info.s2c/momvar.s2c, xlab = expression(sigma[c]^2), ylab = expression(sigma[z]^2),
        main = expression(paste("Asymptotic Relative Efficiency of ", tilde(sigma)[c]^2, " versus ", hat(sigma)[c]^2, "", sep = "")),
        labcex = 1.5)
dev.off()
