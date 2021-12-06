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
#' @param sample_group a vector of length \code{n} indicating whether the samples
#'should be grouped (e.g. paired samples or longitudinal data). Coerced
#'to be a \code{factor}. Default is \code{NULL} in which case no grouping is
#'performed.
#'
#'@param n_perm the number of permutations. Default is \code{100}.
#'
#'@param parallel a logical flag indicating whether parallel computation
#'should be enabled. Default is \code{TRUE}.
#'
#'@param n_cpus an integer indicating the number of cores to be used when
#'\code{parallel} is \code{TRUE}.
#'Default is \code{parallel::detectCores() - 1}.
#'
#'@param space_y a logical flag indicating whether the y thresholds are spaced. 
#'When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and 
#'the maximum of the observations is used. Default is \code{FALSE}.
#'
#'@param number_y an integer value indicating the number of y thresholds (and therefore
#'the number of regressions) to perform the test. Default is \code{length(Y)}.
#'
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @importFrom lme4 lmer
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

test_perm <- function(Y, X, Z = NULL, sample_group=NULL, method=c("linear regression", "mixed model"), n_perm = 100, parallel = FALSE, n_cpus = NULL, space_y = FALSE, number_y = length(Y)){

  if(parallel){
    if(is.null(n_cpus)){
      n_cpus <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(n_cpus)
    doParallel::registerDoParallel(cl)
  }
  else{
    n_cpus <- 1
  }
  if(length(method)>1) method <- method[1]
  stopifnot(method %in% c("linear regression","mixed model"))
  Y <- as.numeric(Y)
  
  if (space_y){
    y <- seq(ifelse(length(which(Y==0))==0,min(Y),min(Y[-which(Y==0)])),max(Y),length.out=number_y)
  }
  else{
    y <- sort(unique(Y))
  }
  
  
  if (is.null(Z)){ 
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    modelmat <- model.matrix(~.,data=X)
  }
  
  # with covariates Z
  else{
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    colnames(Z) <- sapply(1:ncol(Z), function(i){paste0('Z',i)})
    modelmat <- model.matrix(~.,data=cbind(X,Z))
  }
  
  ind_X <- which(substring(colnames(modelmat),1,1)=="X")
  
  beta <- matrix(NA,(length(y)-1),length(ind_X))
  indi_pi <- matrix(NA,length(Y),(length(y)-1))
  
  for (i in 1:(length(y)-1)){ # on fait varier le seuil
    indi_Y <- 1*(Y<=y[i])
    indi_pi[,i] <- indi_Y
    if(method == "linear regression"){
      reg <- lm(indi_Y ~ as.matrix(modelmat[,-1]))
      beta[i,] <- reg$coefficients[ind_X]
    }
    else if(method == "mixed model"){
      if(is.null(sample_group)) {
            warning("sample_group is null",
                "No random effects terms specified in formula")
          }
          else{
            mod_mixed <- lmer(indi_Y ~ 1 + modelmat[,-1] + (1|sample_group))
            beta[i,] <- lme4::fixef(mod_mixed)[ind_X]
      }
    }
  }
  
  beta <- as.vector(beta)
  
  z <- sqrt(length(Y))*(beta[-length(y)])
  STAT_obs <- sum(t(z)*z)
  
  if (parallel){
    if (is.null(Z)){
      
      results <- foreach(i = 1:n_perm, .combine = 'c') %dopar% {
        
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
        STAT_perm <- sum(t(z)*z)
        STAT_perm
      }
    }
    
    else{
      
      if ((ncol(X)!=1)|ncol(Z)!=1){
        stop("X and Z must be univariate. The permutation test will soon be improved to handle multivariate variable of interest and multivariate covariate.")
      }
      
      results <- foreach(i = 1:n_perm, .combine = 'c') %dopar% {
        
        sample_X <- function(X,Z,z){
          X_sample <- rep(NA,length(Z))
          for (i in 1:length(z)){
            X_sample[which(Z==z[i])] <- sample(X[which(Z==z[i])])
          }
          return(X_sample)
        }
        
        X_star <- switch(class(Z[,1]),
                         "factor" = sample_X(X[,1],as.numeric(Z[,1]),unique(as.numeric(Z[,1]))), # vérifier
                         "numeric" = perm_cont(Y,as.numeric(levels(X[,1]))[X[,1]],as.numeric(Z[,1])))
        if (is.factor(X[,1])){
          X_star <- data.frame(X=as.factor(X_star))
        }
        
        X_star <- data.frame(X=X_star)
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
        STAT_perm <- sum(t(z)*z)
        STAT_perm
      }
    }
    parallel::stopCluster(cl)
    
    score <- sum(1*(results>=STAT_obs))
    pval <- (sum(1*(results>=STAT_obs))+1)/(n_perm+1)
  }
  
  else{
    if (is.null(Z)){
      
      STAT_perm <- rep(NA,n_perm)
      
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
    }
    
    else{
      
      STAT_perm <- rep(NA,n_perm)
      
      for(k in 1:n_perm){
        
        sample_X <- function(X,Z,z){
          X_sample <- rep(NA,length(Z))
          for (i in 1:length(z)){
            X_sample[which(Z==z[i])] <- sample(X[which(Z==z[i])])
          }
          return(X_sample)
        }
        
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
    
    score <- sum(1*(STAT_perm>=STAT_obs))
    pval <- (sum(1*(STAT_perm>=STAT_obs))+1)/(n_perm+1)
  }

  return(data.frame(score=score,raw_pval=pval))
  
}


