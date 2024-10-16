#' Asymptotic test for conditional independance
#' 
#' Test the conditional independence of Y and X given Z.
#'
#'@param Y a numeric vector of length \code{n} to test for conditional independence 
#'with \code{X} adjusted on \code{Z}
#'
#'@param X a data frame of size \code{n x p} of numeric or factor vector(s) 
#'containing the variable(s) to be tested for conditional independence 
#'against \code{X} adjusted on \code{Z}.
#' 
#'@param Z a data frame of size \code{n x q} of numeric or factor vector(s) 
#'containing the covariate(s) to condition the independence 
#'test upon.
#'
#'@param space_y a logical flag indicating whether the y thresholds are spaced. 
#'When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and 
#'the maximum of the observations is used. Default is \code{FALSE}.
#'
#'@param number_y an integer value indicating the number of y thresholds (and therefore
#'the number of regressions) to perform the test. Default is \code{n_Y_all}.
#'
#'@importFrom survey pchisqsum
#'
#' @export
#' 
#'@return A data frame with the following elements:
#'\itemize{
#'   \item \code{raw_pval} contains the raw p-values for a given gene.
#'   \item \code{Stat} contains the test statistic for a given gene.
#' }
#' 
#' @examples
#' 
#'X <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
#'Y <- ((X==1)*rnorm(n = 100,0,1)) + ((X==0)*rnorm(n = 100,0.5,1))
#'res_asymp <- cit_asymp(Y,data.frame(X=X))
#'
#'
#' Z <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
#' X <- as.numeric(Z)-1  + rnorm(n=100, sd=1)
#' r <- 1000
#' Y <- replicate(r, as.numeric(Z)-1)
#' YY <- (Y==1)*rnorm(n = 100*r,0,1) + (Y==0)*rnorm(n = 100*r,0.5,1)
#'pvals_sim <- pbapply::pbsapply(1:1000, function(i){
#' res_asymp <- cit_asymp(YY[,i],data.frame(X=X), data.frame(Z=Z))
#' return(res_asymp$raw_pval)
#'})
#'hist(pvals_sim)


cit_asymp <- function(Y, X, Z = NULL, space_y = FALSE, number_y = length(unique(Y))){
  
  # computations independent of Y: should be computed only once before the pbapply loop ----
  
  # no covariates Z
  if (is.null(Z)){ 
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    modelmat <- as.matrix(model.matrix(~.,data=X))
  }
  # with covariates Z
  else{
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    colnames(Z) <- sapply(1:ncol(Z), function(i){paste0('Z',i)})
    modelmat <- as.matrix(model.matrix(~.,data=cbind(X,Z)))
  }
  
  indexes_X <- which(substring(colnames(modelmat), 1, 1) == "X")

  n_Y_all <- length(Y)
  H <- n_Y_all*(solve(crossprod(modelmat)) %*% t(modelmat))[indexes_X, , drop=FALSE]
  # computing the test statistic
  # depends on Y: has to be recomputed for each gene
  Y <- as.numeric(Y) # is this really necessary ??
  
  if (space_y){
    y <- seq(from = ifelse(length(which(Y==0))==0, min(Y), min(Y[-which(Y==0)])),
             to = max(Y[-which.max(Y)]), length.out = number_y)
  }
  else{
    y <- sort(unique(Y))
  }
  p <- length(y) #number of thresholds used
  
  index_jumps <- sapply(y[-p], function(i){sum(Y <= i)})
  beta <- c(apply(X = H[, order(Y), drop=FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, ]) / n_Y_all
  test_stat <- sum(beta^2) * n_Y_all
  
  # Computing the variance ----  
  
  indi_pi <- matrix(NA, n_Y_all, (p-1))
  for (j in 1:(p-1)){ # on fait varier le seuil
    indi_Y <- 1*(Y<=y[j])
    indi_pi[,j] <- indi_Y
  }
  prop <- colMeans(indi_pi)
  

  Sigma2 <- 1/ n_Y_all * tcrossprod(H) %x% (prop - prop %x% t(prop))
  Sigma <- Sigma2*upper.tri(Sigma2, diag = TRUE) +  t(Sigma2*upper.tri(Sigma2, diag = FALSE))

  decomp <- eigen(Sigma, symmetric=TRUE, only.values=TRUE)
  
  # computing the pvalue ----
  pval <- try(survey::pchisqsum(test_stat, lower.tail = FALSE, df = rep(1, ncol(Sigma)), 
                                a = decomp$values, method = "saddlepoint"),
              silent=TRUE)
  if(inherits(pval, "try-error")){
    pval <- try(survey::pchisqsum(test_stat, lower.tail = FALSE, df = rep(1, ncol(Sigma)), 
                                  a = decomp$values, method = "satterthwaite"),
                silent=TRUE)
    if(inherits(pval, "try-error")){
      pval <- NA
    }
  }
  
  return(data.frame("raw_pval" = pval, "Stat" = test_stat))
  
}


