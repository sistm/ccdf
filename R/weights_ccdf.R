#' Computation of weights for CCDF, used in permut function when \code{dist_permuations} is specified.
#'
#'@param Y a numeric vector of size \code{n} containing the
#'preprocessed expressions from \code{n} samples (or cells).
#'
#'@param X a numeric or factor vector(s) of size \code{n}
#'containing the variable to be tested (the condition to be tested). 
#'Multiple variables are not allowed
#' 
#'@param Z a numeric vector of size \code{n}
#'containing the covariate.
#'Multiple variables are not allowed
#'
#' @param sample_group a vector of length \code{n} indicating whether the samples
#'should be grouped (e.g. paired samples or longitudinal data). Coerced
#'to be a \code{factor}. Default is \code{NULL} in which case no grouping is
#'performed.
#'
#'@export
#' 
#'@keywords internal
#' 
#'@return \code{w} a vector of weights
#'
#' @examples
#' 
#'X <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
#'Y <- ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,0.5,1))
#'res_asymp <- weights_ccdf(Y,X)


weights_ccdf <- function(Y,X,Z=NULL,sample_group=NULL){

  find_unique <- function(couple, unique_couple){
    index_list <- list()
    for (j in 1:nrow(unique_couple)){
      index <- NULL
      for (i in 1:nrow(couple)){
        if (identical(as.numeric(couple[i,]),as.numeric(unique_couple[j,]))){
          index <- c(index,i)
        }
      }
      index_list[[j]] <- index
    }
    return(index_list)
  }

  temp_order <- sort(Y,index.return=TRUE)$ix

  if (is.null(Z)){
    couple <- data.frame(Y=Y[temp_order],X=X[temp_order])
  }
  else{
    couple <- data.frame(Y=Y[temp_order],X=X[temp_order],Z=Z[temp_order])
  }
  if(!is.null(sample_group))
    couple <- cbind(couple, sample_group=sample_group[temp_order])

  unique_couple <- unique(couple)
  index_list <- find_unique(couple,unique_couple)

  w <- rep(NA, length(index_list))
  for (i in 1:length(index_list)){
    if (length(index_list[[i]])==1){
      w[index_list[[i]]] <- 1
    }
    else{
      w[index_list[[i]]] <- rep(1/length(index_list[[i]]),length(index_list[[i]]))
    }
  }
  return(w)
}
