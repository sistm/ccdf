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
#'@param n_perm the number of permutations. Default is \code{100}.
#'
#'@param method a character string indicating which method to use to
#'compute the CCDF, either \code{'linear regression'}, \code{'logistic regression'}
#' and  \code{'permutations'} or \code{'RF'} for Random Forests.
#'Default is \code{'linear regression'} since it is the method used in the test.
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

permut <- function(Y, X, Z = NULL, distance = "L2", n_perm = 100, 
                   method="logistic regression", fast=TRUE){


  if (is.null(Z)){

    res <- CDF(Y, X, method=method, fast=fast)
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

    results <- sapply(1:n_perm, function(i){
      X_star <- sample(X)
      res_perm <- CDF(Y=Y, X=X_star, method=method, fast = fast)
      w_perm <- weights_ccdf(Y,X_star)
      switch(distance,
             "L2" = sqrt(sum(w_perm*(res$cdf-res_perm$ccdf)^2)),
             "L1" = sum(w_perm*abs(res$cdf-res_perm$ccdf)),
             "L_sup" = max(w_perm*abs(res$cdf- res_perm$ccdf)))
    })

  }else{

    res_init <-  CDF(Y, X, Z, method=method, fast=fast)
    w_init <- weights_ccdf(Y,X,Z)

    init_dist <- switch(distance,
           "L2" = sqrt(sum(w_init*(res_init$ccdf_nox-res_init$ccdf_x)^2)),
           "L1" = sum(w_init*abs(res_init$ccdf_nox-res_init$ccdf_x)),
           "L_sup" = max(w_init*abs(res_init$ccdf_nox-res_init$ccdf_x))
    )
    
    if (is.factor(Z)){
      Z_ind <- unlist(lapply(unique(Z),function(z){which(Z==z)}))
    }

    sample_X <- function(X,Z,z){
      X_sampled <- rep(NA,length(Z))
      for (zj in unique(z)){
        X_sampled[Z==zj] <- sample(X[Z==zj])
      }
      return(X_sampled)
    }
    
    results <- sapply(1:n_perm, function(i){
      X_star <- switch(class(Z),
                    "factor" = sample_X(X,Z,unique(Z)),
                    "numeric" = perm_cont(Y,X,Z))

      res <- CDF(Y,X_star,Z, method=method, fast=fast)
      w_perm <- weights_ccdf(Y,X_star,Z)

      switch(distance,
             L2 = sqrt(sum(w_perm*(res_init$ccdf_nox-res$ccdf_x)^2)),
             L1 = sum(w_perm*abs(res_init$ccdf_nox-res$ccdf_x)),
             L_sup = max(w_perm*abs(res_init$ccdf_nox-res$ccdf_x)))

      #AD = length(Y)*sum(((res$ccdf_nox-res$ccdf_x)^2)/(res$ccdf_nox*(1-res$ccdf_nox)))
    })
  }
  
  score <- sum(results >= init_dist)
  pval <- (score + 1)/(n_perm + 1)

  return(list(score=score, pval=pval))

}
