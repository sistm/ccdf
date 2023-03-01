#' Permutation test
#'
#'@param Y a numeric vector of size \code{n} containing the
#'preprocessed expression for a given gene from \code{n} samples (or cells).
#'
#'@param X a data frame of numeric or factor vector(s) of size \code{n}
#'containing the variable(s) to be tested (the condition(s)). Multiple variables are not allowed.
#' 
#'@param Z a data frame of numeric or factor vector(s) 
#'of size \code{n} containing the covariate(s). Multiple variables are not allowed.
#'
#'@param n_perm the number of permutations. Default is \code{100}.
#'
#'@param space_y a logical flag indicating whether the y thresholds are spaced. 
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
#'res_perm <- test_perm(Y,data.frame(X=X),n_perm=10)}

test_perm <- function(Y, X, Z = NULL, n_perm = 100, space_y = FALSE, number_y = length(Y)){
  
  
  Y <- as.numeric(Y)
  
  if (space_y){
    y <- seq(ifelse(length(which(Y==0))==0,min(Y),min(Y[-which(Y==0)])),max(Y),length.out=number_y)
  }else{
    y <- sort(unique(Y))
  }
  
  if (is.null(Z)){ 
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    modelmat <- model.matrix(~.,data=X)
  }else{# with covariates Z
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    colnames(Z) <- sapply(1:ncol(Z), function(i){paste0('Z',i)})
    modelmat <- model.matrix(~.,data=cbind(X,Z))
  }
  
  ind_X <- which(substring(colnames(modelmat),1,1) == "X")
  beta <- matrix(NA,(length(y)-1),length(ind_X))
  indi_pi <- matrix(NA,length(Y),(length(y)-1))
  
  for (i in 1:(length(y)-1)){ # on fait varier le seuil
    indi_Y <- 1*(Y<=y[i])
    indi_pi[,i] <- indi_Y
    reg <- lm(indi_Y ~ as.matrix(modelmat[,-1]))
    beta[i,] <- reg$coefficients[ind_X]
  }
  
  beta <- as.vector(beta)
  z <- sqrt(length(Y))*(beta[-length(y)])
  STAT_obs <- sum(t(z)*z)
  
  STAT_perm <- rep(NA,n_perm)
  if (is.null(Z)){
    for (k in 1:n_perm){
      X_star <- data.frame(X=X[sample(1:nrow(X)),])
      #colnames(X_star) <- sapply(1:ncol(X), function(i){paste0('X',i)})
      modelmat_perm <- model.matrix(~.,data=X_star)
      beta_perm <- matrix(NA,(length(y)-1),length(ind_X))
      indi_pi <- matrix(NA,length(Y),(length(y)-1))
      
      for (i in 1:(length(y)-1)){ # on fait varier le seuil
        indi_Y <- 1*(Y<=y[i])
        indi_pi[,i] <- indi_Y
        reg <- lm(indi_Y ~ as.matrix(modelmat_perm[,-1]))
        beta_perm[i,] <- reg$coefficients[ind_X]
      }
      
      beta_perm <- as.vector(beta_perm)
      z <- sqrt(length(Y))*(beta_perm[-length(y)])
      STAT_perm[k] <- sum(t(z)*z)
    }
  }else{
    sample_X <- function(X,Z,z){
      X_sampled <- rep(NA,length(Z))
      for (zj in unique(z)){
        X_sampled[Z==zj] <- sample(X[Z==zj])
      }
      return(X_sampled)
    }
    
    for(k in 1:n_perm){
      X_star <- switch(class(Z[,1]),
                       "factor" = sample_X(X[,1],as.numeric(Z[,1]),unique(as.numeric(Z[,1]))), # vérifier
                       "numeric" = perm_cont(Y,as.numeric(levels(X[,1]))[X[,1]],as.numeric(Z[,1])))
      
      if (is.factor(X[,1])){
        X_star <- data.frame(X=as.factor(X_star))
      }
      
      #colnames(X_star) <- sapply(1:ncol(X), function(i){paste0('X',i)})
      modelmat_perm <- model.matrix(~.,data=cbind(X_star,Z))
      beta_perm <- matrix(NA,(length(y)-1),length(ind_X))
      indi_pi <- matrix(NA,length(Y),(length(y)-1))
      
      for (i in 1:(length(y)-1)){
        indi_Y <- 1*(Y<=y[i])
        indi_pi[,i] <- indi_Y
        reg <- lm(indi_Y ~ as.matrix(modelmat_perm[,-1]))
        beta_perm[i,] <- reg$coefficients[ind_X]
      }
      
      beta_perm <- as.vector(beta_perm)
      z <- sqrt(length(Y))*(beta_perm[-length(y)])
      STAT_perm[k] <- sum(t(z)*z)
    }
  }
  
  score <- sum(STAT_perm >= STAT_obs)
  pval <- (score + 1) / (n_perm + 1)
  return(data.frame(score=score, raw_pval=pval))
  
}


