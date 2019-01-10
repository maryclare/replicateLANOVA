library(statmod)
library(Matrix)

sample.beta <- function(X = X, y = y, sigma.sq, tau, Xty = NULL) {
  
  p <- ncol(X)
  n <- nrow(X)
  
  diag.X <- isDiagonal(X)
  
  if (!diag.X) {
    D <- diag(tau^2)
    Phi <- X/sqrt(sigma.sq)
    alpha <- y/sqrt(sigma.sq)
  
    u <- tau*rnorm(p)
    delta <- rnorm(n)
    v <- crossprod(t(Phi), u) + delta
    w <- crossprod(solve(tcrossprod(tcrossprod(Phi, D), Phi) + diag(n)), alpha - v) # This step could be made better
    beta <- u + crossprod(D, crossprod(Phi, w))
  } else {
    XtX.diag <- diag(X)^2
    var <- 1/(XtX.diag/sigma.sq + 1/tau^2)
    mean <- var*(Xty/sigma.sq)
    beta <- mean + sqrt(var)*rnorm(length(var))
  }
  
  return(beta) 
  
}

sample.tau <- function(beta, lambda) {
  
  p <- length(beta)
  mu <- sqrt(lambda^2/beta^2)
  shape <- lambda^2
  
  tau.inv <- sqrt(rinvgauss(p, mean = mu, shape = rep(shape, p)))
  
  return(1/tau.inv)
  
}


bayes.lasso <- function(X = X, y = y, sigma.sq, lambda, num.samp, delta = 0) {
  
  p <- ncol(X)
  Xty <- crossprod(X, y)
  
  beta <- crossprod(solve(crossprod(X) + delta*diag(p)), crossprod(X, y))
  
  samples.beta <- samples.tau <- matrix(nrow = num.samp, ncol = p)
  lik <- numeric(num.samp)
  
  
  for (i in 1:num.samp) {
    cat("i = ", i, "\n")
    tau <- sample.tau(beta = beta, lambda = lambda)
    beta <- sample.beta(X = X, y = y, tau = tau, sigma.sq = sigma.sq, Xty = Xty)
    samples.beta[i, ] <- beta
    samples.tau[i, ] <- tau
    lik[i] <- exp(sum(dnorm(y - crossprod(t(X), beta), mean = 0, sd = sqrt(sigma.sq), log = TRUE)) + 
                    p*log(lambda/2) - (lambda*sum(abs(beta))))
  }
  
  return(list("samples.beta" = samples.beta,
              "samples.tau" = samples.tau,
              "lik" = lik))
  
}

bayes.ridge <- function(X = X, y = y, sigma.sq, tau.sq, num.samp, delta = 0) {
  
  p <- ncol(X)
  
  beta <- crossprod(solve(crossprod(X) + delta*diag(p)), crossprod(X, y))
  
  samples.beta <- matrix(nrow = num.samp, ncol = p)
  lik <- numeric(num.samp)
  
  
  for (i in 1:num.samp) {
    beta <- sample.beta(X = X, y = y, tau = sqrt(tau.sq)*rep(1, p), sigma.sq = sigma.sq)
    samples.beta[i, ] <- beta
    lik[i] <- exp(sum(dnorm(y - crossprod(t(X), beta), mean = 0, sd = sqrt(sigma.sq), log = TRUE)) + 
                    sum(dnorm(beta, mean = 0, sd = sqrt(tau.sq), log = TRUE)))
  }
  
  return(list("samples.beta" = samples.beta,
              "lik" = lik))
  
}

bayes.unifo <- function(X = X, y = y, sigma.sq, num.samp, delta = 0) {
  
  p <- ncol(X)
  
  beta <- crossprod(solve(crossprod(X) + delta*diag(p)), crossprod(X, y))
  
  samples.beta <- matrix(nrow = num.samp, ncol = p)
  lik <- numeric(num.samp)
  
  Phi <- X/sqrt(sigma.sq)
  V <- solve(crossprod(Phi))
  alpha <- y/sqrt(sigma.sq)
  mu <- crossprod(V, crossprod(Phi, alpha))
  V.half <- chol(V)
  
  for (i in 1:num.samp) {
    beta <- crossprod(V.half, rnorm(p)) + mu
    samples.beta[i, ] <- beta
    lik[i] <- exp(sum(dnorm(y - crossprod(t(X), beta), mean = 0, sd = sqrt(sigma.sq), log = TRUE)))
  }
  
  return(list("samples.beta" = samples.beta,
              "lik" = lik))
  
}
