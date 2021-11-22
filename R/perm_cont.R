#' Permutation procedure when Z is continuous
#'
#'
#'@param Y a numeric vector of size \code{n} containing the
#'preprocessed expressions from \code{n} samples (or cells).
#'
#'@param X a numeric or factor vector of size \code{n}
#'containing the variable to be tested (the condition to be tested). 
#' 
#'@param Z a numeric vector of size \code{n}
#'containing the covariate. Multiple variables are not allowed.
#'
#' @param sample_group a vector of length \code{n} indicating whether the samples
#'should be grouped (e.g. paired samples or longitudinal data). Coerced
#'to be a \code{factor}. Default is \code{NULL} in which case no grouping is
#'performed.
#'
#'@export
#'
#'@import stats
#'
#' @importFrom lme4 lmer
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



perm_cont <- function(Y,X,Z, sample_group=NULL, method=c("linear regression", "mixed model")){
  prob <- matrix(NA,length(Z),length(Z))
  new_prob <- matrix(NA,length(Z),length(Z))
  if(length(method)>1) method <- method[1]
  if (method=="mixed model" & is.null(sample_group)) {
      warning("Warning")
      break
  }
  else if(method=="mixed model" & !is.null(sample_group))
      reg <- lmer(X ~ Z + (1|sample_group))
  else reg <- lm(X~Z)
  X_star <- rep(NA,length(X))

  for (z in 1:length(Z)){
    for (l in 1:length(Z)){
      prob[z,l] <- (abs(predict(reg,newdata=data.frame(Z=Z[z]))-predict(reg,newdata=data.frame(Z=Z[l])))/
          (sum(abs(predict(reg,newdata=data.frame(Z=Z[z]))-predict(reg,newdata=data.frame(Z=Z))))))[1]
    }
    new_prob[z,-z] <- (1/prob[z,-z])/sum(1/prob[z,-z])
    new_prob[z,z] <- 0
  }

  for (i in 1:length(X)){
    X_star[i] <- sample(X[-i], size=1, prob = prob[i,-i])
  }

  return(X_star)
}
