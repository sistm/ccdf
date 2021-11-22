#' Permutation test when \code{dist_permutations} is specified
#'
#'@param Y a numeric vector of size \code{n} containing the
#'preprocessed expressions from \code{n} samples (or cells).
#'
#'@param X a data frame containing numeric or factor vector(s) of size \code{n}
#'containing the variable(s) to be tested (the condition(s) to be tested). 
#' 
#'@param Z a data frame containing numeric or factor vector(s) of size \code{n}
#'containing the covariate(s).
#'
#' @param sample_group a vector of length \code{n} indicating whether the samples
#'should be grouped (e.g. paired samples or longitudinal data). Coerced
#'to be a \code{factor}. Default is \code{NULL} in which case no grouping is
#'performed.
#'
#'@param method a character string indicating which method to use to
#'compute the CCDF, either \code{'linear regression'}, \code{'logistic regression'}
#' and  \code{'permutations'}, \code{'RF'} for Random Forests or \code{'mixed model'} for mixture model
#'Default is \code{'linear regression'} since it is the method used in the test.
#'
#'#'@param parallel a logical flag indicating whether parallel computation
#'should be enabled. Default is \code{TRUE}.
#'
#'@param n_cpus an integer indicating the number of cores to be used when
#'\code{parallel} is \code{TRUE}.
#'Default is \code{parallel::detectCores() - 1}.
#'
#'@param fast a logical flag indicating whether the fast implementation of
#'logistic regression should be used. Only if \code{'dist_permutations'} is specified.
#'Default is \code{TRUE}.
#'
#'@return A data frame with the following elements:
#'\itemize{
#'   \item \code{score} contains the test statistic for a given gene.
#'   \item \code{pval} contains the raw p-values for a given gene computed from \code{n_perm} permutations.
#' }
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export
#' 
#' @keywords internal
#' 
#' @examples
#' 
#'if(interactive()){
#'X <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
#'Y <- ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,0.5,1))
#'permut(Y,X,method="linear regression",n_perm=10,n_cpus=2)}

permut <- function(Y, X, Z = NULL, sample_group=NULL,distance = "L2", n_perm, method="logistic regression",
                     parallel = TRUE, n_cpus = NULL, fast=TRUE){

  if(parallel){
    if(is.null(n_cpus)){
      n_cpus <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(n_cpus)
    doParallel::registerDoParallel(cl)
  }else{
    cl <- NULL
  }


  if (is.null(Z)){

    res <- CDF(Y, X, sample_group=sample_group, method=method, fast=fast)
    w_init <- weights_ccdf(Y,X)

    if (distance=="L2"){
      init_dist <- sqrt(sum(w_init*(res$cdf-res$ccdf)^2))
    }
    if (distance=="L1"){
      init_dist <- sum(w_init*abs(res$cdf-res$ccdf))
    }
    if (distance=="L_sup"){
      init_dist <- max(w_init*abs(res$cdf- res$ccdf))
    }

    results <- foreach(i = 1:n_perm, .combine = 'c') %dopar% {
      X_star <- sample(X)
      res_perm <- CDF(Y=Y, X=X_star, sample_group=sample_group, method=method, fast = fast)
      w_perm <- weights_ccdf(Y,X_star)
      switch(distance,
             "L2" = sqrt(sum(w_perm*(res$cdf-res_perm$ccdf)^2)),
             "L1" = sum(w_perm*abs(res$cdf-res_perm$ccdf)),
             "L_sup" = max(w_perm*abs(res$cdf- res_perm$ccdf)))
    }

  }

  else{

    res_init <-  CDF(Y, X, Z, sample_group, method=method, fast=fast)
    w_init <- weights_ccdf(Y,X,Z)

    init_dist <- sqrt(sum(w_init*(res_init$ccdf_nox-res_init$ccdf_x)^2))

    if (distance=="L2"){
      init_dist <- sqrt(sum(w_init*(res_init$ccdf_nox-res_init$ccdf_x)^2))
    }
    if (distance=="L1"){
      init_dist <- sum(w_init*abs(res_init$ccdf_nox-res_init$ccdf_x))
    }
    if (distance=="L_sup"){
      init_dist <- max(w_init*abs(res_init$ccdf_nox-res_init$ccdf_x))
    }

    if (is.factor(Z)){
      Z_ind <- unlist(lapply(unique(Z),function(z){which(Z==z)}))
    }

    results <- foreach(i = 1:n_perm, .combine = 'c') %dopar% {

      sample_X <- function(X,Z,z){
        X_sample <- rep(NA,length(Z))
        for (i in 1:length(z)){
          X_sample[which(Z==z[i])] <- sample(X[which(Z==z[i])])
        }
        return(X_sample)
      }

      X_star <- switch(class(Z),
                    "factor" = sample_X(X,Z,unique(Z)),
                    "numeric" = perm_cont(Y,X,Z))

      res <- CDF(Y,X_star,Z, sample_group, method=method, fast=fast)
      w_perm <- weights_ccdf(Y,X_star,Z)

      switch(distance,
             L2 = sqrt(sum(w_perm*(res_init$ccdf_nox-res$ccdf_x)^2)),
             L1 = sum(w_perm*abs(res_init$ccdf_nox-res$ccdf_x)),
             L_sup = max(w_perm*abs(res_init$ccdf_nox-res$ccdf_x)))
    }
  }
  
  score <- sum(1*(results>=init_dist))+1
  pval <- (sum(1*(results>=init_dist))+1)/(n_perm+1)

  if(parallel){
    parallel::stopCluster(cl)
  }
  return(list(score=score,pval=pval))

}
