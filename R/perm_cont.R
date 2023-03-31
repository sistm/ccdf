#' Permutation procedure when Z is continuous
#'
#'
#'@param Y a numeric vector of length \code{n} containing the
#'preprocessed expressions from \code{n} samples (or cells).
#'
#'@param X a numeric or factor vector of length \code{n}
#'containing the variable to be tested (the condition to be tested). 
#' 
#'@param Z a numeric vector of length \code{n}
#'containing the covariate. Multiple variables are not allowed.
#'
#'@export
#'
#'@import stats
#' 
#'@return \code{X_star} a vector of permuted \code{X}.
#'
#'@examples
#' 
#'if(interactive()){
#'X <- rbinom(n=100, size = 1, prob = 0.5)
#'Z <- rnorm(100,0,1)
#'Y <- ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,0.5,1))
#'res <- perm_cont(Y,X,Z)}



perm_cont <- function(Y,X,Z){
  prob <- matrix(0,length(Z),length(Z))
  modmat <- model.matrix(~Z)
  reg_coefs <- solve(crossprod(modmat)) %*% t(modmat) %*% X
  X_star <- rep(NA,length(X))
  n <- length(Z)
  
  fit <- modmat %*% reg_coefs
  for (z in 1:length(Z)){
      pred_z <- fit[z,1]
      prob[z, -z] <- (pred_z/(pred_z - fit[-z]))^2
  }
  
  for (i in 1:length(X)){
    X_star[i] <- sample(X[-i], size=1, prob = prob[i,-i])
  }

  return(X_star)
}
