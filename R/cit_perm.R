#' Permutation test for conditional independence
#'
#'@param Y a numeric vector of length \code{n} to test for conditional independence 
#'with \code{X} adjusted on \code{Z}
#'
#'@param X a data frame of size \code{n x p} of numeric or factor vector(s) 
#'containing the variable(s) to be tested for conditional independence 
#'against \code{X} adjusted on \code{Z}. Multiple variables are not supported 
#'for permutation.
#' 
#'@param Z a data.frame of size \code{n x 1} of numeric or factor vector 
#'containing the covariate to condition the independence 
#'test upon. Multiple covariates are not supported for permutation.
#'
#'@param n_perm the number of permutations. Default is \code{100}.
#'
#'@param space_y a logical flag indicating whether the y thresholds are spaced out. 
#'When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and 
#'the maximum of the observations is used. Default is \code{FALSE}.
#'
#'@param number_y an integer value indicating the number of y thresholds (and therefore
#'the number of regressions) to perform the test. Default is \code{length(Y)}.
#'
#' @export
#' 
#'@return A data frame with the following elements:
#'\itemize{
#'   \item \code{score} contains the test statistic for a given gene.
#'   \item \code{raw_pval} contains the raw p-values for a given gene computed from \code{n_perm} permutations.
#' }
#' 
#' @examples
#' 
#'if(interactive()){
#'X <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
#'Y <- ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,0.5,1))
#'res_perm <- cit_perm(Y,data.frame(X=X),n_perm=10)}

cit_perm <- function(Y, X, Z = NULL, n_perm = 100, space_y = FALSE, number_y = length(Y)){
  
  stopifnot(is.vector(Y))
  stopifnot(is.data.frame(X))
  stopifnot(is.data.frame(Z) | is.null(Z))
  stopifnot(ncol(X) < 2)
  stopifnot(ncol(Z) < 2 | is.null(Z))
  
  n <- length(Y)
  stopifnot(nrow(X) == n)
  stopifnot(nrow(Z) == n | is.null(Z))
  
  

  if (is.null(Z)){ 
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    modelmat <- model.matrix(~.,data=X)
  }else{# with covariates Z
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    colnames(Z) <- sapply(1:ncol(Z), function(i){paste0('Z',i)})
    modelmat <- model.matrix(~.,data=cbind(X,Z))
  }
  
  indexes_X <- which(substring(colnames(modelmat), 1, 1) == "X")
  
  n_Y_all <- length(Y)
  H <- n_Y_all*(solve(crossprod(modelmat)) %*% t(modelmat))[indexes_X, , drop=FALSE]
  
  Y <- as.numeric(Y)
  oY <- order(Y)
  
  if (space_y){
    y <- seq(from = ifelse(length(which(Y==0))==0, min(Y), min(Y[-which(Y==0)])),
             to = max(Y), length.out = number_y)
  }else{
    y <- sort(unique(Y))
  }
  
  p <- length(y) #number of thresholds used
  
  index_jumps <- sapply(y[-p], function(i){sum(Y <= i)})
  beta <- c(apply(X = H[, oY, drop=FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, ]) / n_Y_all
  test_stat_obs <- sum(beta^2) * n_Y_all

  
  test_stat_perm <- rep(NA,n_perm)
  if (is.null(Z)){
    for (k in 1:n_perm){
      X_star <- data.frame(X=X[sample(1:nrow(X)),])
      #colnames(X_star) <- sapply(1:ncol(X), function(i){paste0('X',i)})
      modelmat_perm <- model.matrix(~.,data=X_star)
      
      H_perm <- n_Y_all*(solve(crossprod(modelmat_perm)) %*% t(modelmat_perm))[indexes_X, , drop=FALSE]
      beta_perm <- c(apply(X = H_perm[, oY, drop=FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, ]) / n_Y_all
      
      test_stat_perm[k] <- sum(beta_perm^2) * n_Y_all
    }
  }else{
    sample_X <- function(X,Z,z){
      X_sampled <- rep(NA,length(Z))
      for (zj in unique(z)){
        X_sampled[Z==zj] <- sample(X[Z==zj])
      }
      return(X_sampled)
    }
    
    #browser()
    for(k in 1:n_perm){
      X_star <- switch(class(Z[,1]),
                       "factor" = sample_X(X[,1], as.numeric(Z[,1]), unique(as.numeric(Z[,1]))), # vérifier
                       "integer" = sample_X(X[,1], as.numeric(Z[,1]), unique(as.numeric(Z[,1]))), # vérifier
                       "numeric" = perm_cont(Y = Y, X = if(is.factor(X[,1])){as.numeric(levels(X[,1]))[X[,1]]}else{as.numeric(X[,1])}, Z=Z[,1])
      )
      if (is.factor(X[,1])){
        X_star <- data.frame(X=as.factor(X_star))
      }
      
      
      #colnames(X_star) <- sapply(1:ncol(X), function(i){paste0('X',i)})
      modelmat_perm <- model.matrix(~.,data=cbind(X_star,Z))
      
      H_perm <- n_Y_all*(solve(crossprod(modelmat_perm)) %*% t(modelmat_perm))[indexes_X, , drop=FALSE]
      beta_perm <- c(apply(X = H_perm[, oY, drop=FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, ]) / n_Y_all
      
      test_stat_perm[k] <- sum(beta_perm^2) * n_Y_all
    }
  }
  
  score <- sum(test_stat_perm >= test_stat_obs)
  pval <- (score + 1) / (n_perm + 1)
  return(data.frame(score=score, raw_pval=pval))
  
}


