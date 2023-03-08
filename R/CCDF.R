#' Estimate the empirical conditional cumulative distribution function
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
#'@param method a character string indicating which method to use to
#'compute the CCDF, either \code{'OLS'} or \code{'logistic'}.
#'Default is \code{'OLS'} for greater computational speed.
#'
#'@param fast a logical flag indicating whether the fast implementation of
#'logistic regression should be used. Only if \code{method == 'logistic'}.
#'Default is \code{TRUE}.
#'
#'@param space_y a logical flag indicating whether the y thresholds are spaced. 
#'When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and 
#'the maximum of the observations is used. Default is \code{FALSE}.
#'
#'@param number_y an integer value indicating the number of y thresholds (and therefore
#'the number of regressions) to perform the test. Default is \code{length(Y)}.
#'
#'@importFrom stats model.matrix
#' 
#'@export
#' 

#'@return A list with the following elements:\itemize{
#'   \item \code{cdf}: a vector of the cumulative distribution function of a given gene.
#'   \item \code{ccdf}: a vector of the conditional cumulative distribution function of a given gene, computed
#'   given \code{X}. Only if \code{Z} is \code{NULL}.
#'   \item \code{ccdf_nox}: a vector of the conditional cumulative distribution function of a given gene, computed
#'   given \code{Z} only (i.e. \code{X} is ignored.). Only if \code{Z} is not \code{NULL}.
#'   \item \code{ccdf_x}: a vector of the conditional cumulative distribution function of a given gene, computed
#'   given \code{X} and \code{Z}. Only if \code{Z} is not \code{NULL}.
#'   \item \code{y_sort}: a vector of the sorted expression points at which the CDF and the CCDFs are calculated.
#'   \item \code{x_sort}: a vector of the variables associated with \code{y_sort}.
#'   \item \code{z_sort}: a vector of the covariates associated with \code{y_sort}. Only if \code{Z} is not \code{NULL}.
#' }
#' 
#'@examples
#' 
#'X <- as.factor(rbinom(n=1000, size = 1, prob = 0.5))
#'Y <- ((X==1)*rnorm(n = 500,0,1)) + ((X==0)*rnorm(n = 500,0.5,1))
#'res <- ccdf(Y,data.frame(X=X),method="OLS")

ccdf <- function(Y, X, Z=NULL, method=c("OLS","logistic"), 
                 fast=TRUE, space_y=FALSE, number_y=length(Y)){
  
  if (length(method) > 1) {
    method <- method[1]
  }
  stopifnot(method %in% c("OLS","logistic"))
  
  
  if(!is.numeric(Y)){
    warning("Converting Y to a numeric vector.\n", 
            "This should have been done beforehand.")
    Y <- as.numeric(Y)
  }
  
  if (sum(is.na(Y)) > 0) {
    warning("`Y` contains ", sum(is.na(Y)), " NA values. ",
            "\nCurrently they are ignored in the computations but ",
            "you should think carefully about where do those NA/NaN ",
            "come from...")
    Y <- Y[stats::complete.cases(Y)]
  }
  
  if(space_y){
    y <- seq(from = min(Y[ - which(Y==min(Y))]), 
             to = max(Y), 
             length.out = number_y)
  }else{
    y <- sort(unique(Y))
  }
  
  output <- NULL
  
  if (is.null(Z)){
    n_Y <- length(Y)
    # temp_order <- sort(Y,index.return=TRUE)$ix
    # y <- sort(unique(Y))
    # y_sort <- sort(Y)
    # x_sort <- X[temp_order]
    #modelmat <- model.matrix(Y~X)
    
    
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    modelmat <- model.matrix(~.,data=X)
    
    
    ind_X <- which(substring(colnames(modelmat),1,1)=="X")
    
    cdf <- list()
    ccdf <- list()
    x_sort <- list()
    y_sort <- list()
    
    for (i in 1:(length(y)-1)){
      Ylow <- Y<=y[i]
      if (i==1){
        w <- Ylow
      }else{
        w <- Ylow & (Y>y[i-1])
      }
      x_sort[[i]] <- X[w,]
      y_sort[[i]] <- Y[w]
      indi_Y <- 1*Ylow
      
      # unCDF
      cdf[[i]] <- rep(sum(indi_Y)/n_Y, sum(w))
      
      if (length(unique(indi_Y))==1){
        ccdf[[i]] <- rep(1, sum(w))
      }else{
        
        if (method=="logistic"){
          # CDF
          if(fast){
            #fast
            glm_coef <- RcppNumerical::fastLR(x=modelmat, y=indi_Y,
                                              eps_f = 1e-08, eps_g = 1e-08)$coefficients
          }else{
            #safe
            glm_coef <- glm(indi_Y ~. , data=cbind.data.frame(indi_Y, X), family = binomial(link = "logit"))$coefficients
          }
          exp_predlin <- exp(-apply(t(glm_coef*t(modelmat[w, , drop=FALSE])), 1, sum))
          ccdf[[i]] <- 1/(1+exp_predlin) # exp_predlin/(1+exp_predlin)
        }else if (method=="OLS"){
          coefs_ols <- solve(crossprod(modelmat)) %*% t(modelmat) %*% indi_Y
          ccdf[[i]] <- (modelmat[w,] %*% coefs_ols)[ ,1]
        }
      }
    }
    
    
    ccdf_unlisted <- unlist(ccdf, use.names = FALSE)
    cdf_unlisted <- unlist(cdf, use.names = FALSE)
    output <- list(cdf=cdf_unlisted, ccdf=ccdf_unlisted, y=unlist(y_sort), x=unlist(x_sort))
    class(output) <- "ccdf"
    
  }else{
    n_Y <- length(Y)
    
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    colnames(Z) <- sapply(1:ncol(Z), function(i){paste0('Z',i)})
    modelmat <- model.matrix(~., data=cbind.data.frame(X, Z))
    
    ind_X <- which(substring(colnames(modelmat),1,1)=="X")
    
    x_sort <- list()
    y_sort <- list()
    z_sort <- list()
    
    ccdf_x <- list()
    ccdf_nox <- list()
    cdf <- list()
    
    for (i in 1:(length(y)-1)){
      
      #new_data <- data.frame(X[w],Z[w])
      #names(new_data) <- c("X","Z")
      
      #cdf[[i]] <- rep(sum(indi_Y)/n_Y, sum(w))
      
      Ylow <- Y<=y[i]
      
      if (i==1){
        w <- Ylow
      }
      else{
        w <- Ylow & (Y>y[i-1])
      }
      
      x_sort[[i]] <- X[w,]
      y_sort[[i]] <- Y[w]
      z_sort[[i]] <- Z[w,]
      indi_Y <- 1*Ylow
      
      # unCDF
      cdf[[i]] <- rep(sum(indi_Y)/n_Y, sum(w))
      
      if (length(unique(indi_Y))==1){
        ccdf[[i]] <- rep(1,sum(w))
      }
      
      if (method=="logistic"){
        if (fast){
          glm_coef_x <- RcppNumerical::fastLR(x=modelmat, y=indi_Y)$coefficients
          glm_coef_nox <- RcppNumerical::fastLR(x= modelmat[,-ind_X], y=indi_Y)$coefficients
          
        }
        else{
          glm_coef_x <- glm.fit(x=modelmat,y=indi_Y, family = binomial())$coefficients
          glm_coef_nox <- glm.fit(x=modelmat[,-ind_X],y=indi_Y, family = binomial())$coefficients
        }
          exp_predlin_x <- exp(-apply(t(glm_coef_x*t(modelmat[w, , drop=FALSE])), 1, sum))
          exp_predlin_nox <- exp(-apply(t(glm_coef_nox*t(modelmat[w, -ind_X, drop=FALSE])), 1, sum))
        
        ccdf_x[[i]] <- 1/(1+exp_predlin_x) # exp_predlin_x/(1+exp_predlin_x)
        ccdf_nox[[i]] <- 1/(1+exp_predlin_nox) # exp_predlin_nox/(1+exp_predlin_nox)
        
      }
      
      else if (method=="OLS"){
        coefs_ols_x <- solve(crossprod(modelmat)) %*% t(modelmat) %*% indi_Y
        ccdf_x[[i]] <-(modelmat[w, ] %*% coefs_ols_x)[ ,1]
        coefs_ols_nox <- solve(crossprod(modelmat[,-ind_X])) %*% t(modelmat[,-ind_X]) %*% indi_Y
        ccdf_nox[[i]] <- (modelmat[w,-ind_X] %*% coefs_ols_nox)[ ,1]
      }
    }
    
    ccdf_x_unlisted <- unlist(ccdf_x, use.names = FALSE)
    ccdf_nox_unlisted <- unlist(ccdf_nox, use.names = FALSE)
    cdf_unlisted <- unlist(cdf, use.names = FALSE)
    
    output <-  list(cdf=cdf_unlisted, ccdf_nox=ccdf_nox_unlisted, ccdf_x=ccdf_x_unlisted, y=unlist(y_sort), x=unlist(x_sort), z=unlist(z_sort))
    class(output) <- "ccdf"
  }
  
  
  return(output)
}

