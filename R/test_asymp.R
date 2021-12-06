#' Asymptotic test
#'
#'@param Y a numeric vector of size \code{n} containing the
#'preprocessed expression for a given gene from \code{n} samples (or cells).
#'
#'@param X a data frame of numeric or factor vector(s) of size \code{n}
#'containing the variable(s) to be tested (the condition(s))
#' 
#'@param Z a data frame of numeric or factor vector(s) 
#'of size \code{n} containing the covariate(s)
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
#'Y <- ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,0.5,1))
#'res_asymp <- test_asymp(Y,data.frame(X=X))


test_asymp <- function(Y, X, Z = NULL, space_y = FALSE, number_y = length(unique(Y))){
  
  Y <- as.numeric(Y)
  n_Y_all <- length(Y)
  
  if (space_y){
    y <- seq(ifelse(length(which(Y==0))==0,min(Y),min(Y[-which(Y==0)])),max(Y[-which.max(Y)]),length.out=number_y)
  }
  else{
    y <- sort(unique(Y))
  }
  n_y_unique <- length(y)
  
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
  
  indexes_X <- which(substring(colnames(modelmat),1,1)=="X")
  p_X <- length(indexes_X)
  
  beta <- matrix(NA, (n_y_unique-1), p_X)
  indi_pi <- matrix(NA, n_Y_all, (n_y_unique-1))
  
  Phi <- (1/n_Y_all)*(t(modelmat)%*%modelmat)
  H <- (solve(Phi)%*%t(modelmat)) # ginv
  H <- H[indexes_X, , drop=FALSE]
  
  for (i in 1:(n_y_unique-1)){ # on fait varier le seuil
    indi_Y <- 1*(Y<=y[i])
    indi_pi[,i] <- indi_Y
    modX_OLS <- modelmat[, c(1, indexes_X), drop=FALSE]
    beta[i,] <- (solve(crossprod(modX_OLS))  %*% t(modX_OLS) %*% indi_Y)[indexes_X]
  }
  
  beta <- as.vector(beta)
  prop <- colMeans(indi_pi)
  
  temp_Sigma <-  lapply(1:ncol(H), 
                        function(k){sapply(1:nrow(H), 
                                           function(s){sapply(1:nrow(H), 
                                                              function(r){H[s,k]*H[r,k]}
                                           )}
                        )}
  )
  sum_temp_Sigma <- Reduce(f = `+`, x = temp_Sigma)
  if(is.null(dim(sum_temp_Sigma))){
    sum_temp_Sigma <- matrix(sum_temp_Sigma)
  }
  ind_sig <- rep(1:(n_y_unique-1), p_X)
  
  Sigma <- matrix(data = NA, nrow = (n_y_unique-1)*p_X, 
                  ncol = (n_y_unique-1)*p_X)
  for (i in 1:((n_y_unique-1)*p_X)){
    for (j in 1:i){
      Sigma[i,j] <- Sigma[j,i] <- (1/n_Y_all)*sum_temp_Sigma[floor(i/(n_y_unique)+1),floor(j/(n_y_unique)+1)]*(prop[ind_sig[j]]-(prop[ind_sig[j]]*prop[ind_sig[i]]))
    }
  }
  
  decomp <- eigen(Sigma)
  
  z <- (sqrt(n_Y_all))*beta
  test_stat <- sum(t(z)*z)
  
  pval <- survey::pchisqsum(test_stat, lower.tail = FALSE, df = rep(1,ncol(Sigma)), 
                            a = decomp$values, method = "saddlepoint")
  
  return(data.frame("raw_pval" = pval, "Stat" = test_stat))
  
}


