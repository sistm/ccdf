#' Old version of CDF function
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
#'compute the CCDF, either \code{'linear regression'}, \code{'logistic regression'}
#' and  \code{'permutations'} or \code{'RF'} for Random Forests.
#'Default is \code{'linear regression'} since it is the method used in the test.
#'
#'@param fast a logical flag indicating whether the fast implementation of
#'logistic regression should be used. Only if \code{'dist_permutations'} is specified.
#'Default is \code{TRUE}.
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
#' @import RcppNumerical
#' @importFrom  randomForest randomForest
#' @import rpart
#' 
#' @export
#' 
#' @keywords internal
#' @examples
#' 
#'X <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
#'Y <- ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,0.5,1))
#'res <- CDF(Y,X,method="linear regression")
#' 


CDF <- function(Y,X,Z=NULL,method="linear regression", fast=TRUE){

  if (is.null(Z)){
    n_Y <- length(Y)
    temp_order <- sort(Y,index.return=TRUE)$ix
    y <- sort(unique(Y))
    y_sort <- sort(Y)
    x_sort <- X[temp_order]
    modelmat <- model.matrix(Y~X)

    cdf <- list()
    ccdf <- list()

    for (i in 1:length(y)){
      w <- which(Y==y[i])
      indi_Y <- 1*(Y<=y[i])

      # unCDF
      cdf[[i]] <- rep(sum(indi_Y)/n_Y, length(w))

      if (length(unique(indi_Y))==1){
        ccdf[[i]] <- rep(1,length(w))
      }

      else{

        if (method=="RF"){
          # CCDF
          rf <- randomForest(data.frame(X=X, row.names = NULL), indi_Y,ntree=100)
          ccdf[[i]] <- predict(rf)[w]
        }

        else if (method=="logistic regression"){
          # CDF
          if(fast){
            #fast
            glm_coef <- RcppNumerical::fastLR(x=modelmat, y=indi_Y,
                                              eps_f = 1e-08, eps_g = 1e-08)$coefficients
          }else{
            #safe
            glm_coef <- glm(indi_Y ~ 1 + X, family = binomial(link = "logit"))$coefficients
          }
          if (is.null(dim(modelmat[w,]))){
            exp_predlin <- exp(-sum(t(glm_coef*t(modelmat[w, ]))))
          }
          else{
            exp_predlin <- exp(-apply(t(glm_coef*t(modelmat[w, ])),1,sum))
          }
          ccdf[[i]] <- 1/(1+exp_predlin) # exp_predlin/(1+exp_predlin)
        }

        else if (method=="linear regression"){
          reg <- lm(indi_Y ~ 1 + modelmat[,2])
          ccdf[[i]] <- predict(reg)[w]
        }
      }
    }

    ccdf <- unlist(ccdf, use.names = FALSE)
    cdf <- unlist(cdf, use.names = FALSE)
    return(list(cdf=cdf, ccdf=ccdf, y=y_sort, x=x_sort))
  }

  else{

    n_Y <- length(Y)
    temp_order <- sort(Y,index.return=TRUE)$ix
    y <- sort(unique(Y))
    y_sort <- sort(Y)
    x_sort <- X[temp_order]
    z_sort <- Z[temp_order]
    modelmat <- model.matrix(Y~X+Z)

    ccdf_x <- list()
    ccdf_nox <- list()
    cdf <- list()

    for (i in 1:length(y)){
      w <- which(Y==y[i])
      indi_Y <- 1*(Y<=y[i])
      new_data <- data.frame(X[w],Z[w])
      names(new_data) <- c("X","Z")

      cdf[[i]] <- rep(sum(indi_Y)/n_Y, length(w))

      if (method=="RF"){
        # CDF
        rf_x <- randomForest(data.frame(X=X, Z=Z, row.names = NULL), indi_Y,ntree=50)
        ccdf_x[[i]] <- predict(rf_x)[w]
        # CCDF
        new_data <- data.frame(Z[w])
        names(new_data) <- c("Z")
        rf_nox <- randomForest(data.frame(Z=Z, row.names = NULL), indi_Y,ntree=50)
        ccdf_nox[[i]] <-  predict(rf_nox)[w]

      }

      else if (method=="logistic regression"){
        if (fast){
          glm_coef_x <- RcppNumerical::fastLR(x=modelmat, y=indi_Y)$coefficients
          glm_coef_nox <- RcppNumerical::fastLR(x=modelmat[,-2], y=indi_Y)$coefficients

        }
        else{
          glm_coef_x <- glm.fit(x=modelmat,y=indi_Y, family = binomial())$coefficients
          glm_coef_nox <- glm.fit(x=modelmat[,-2],y=indi_Y, family = binomial())$coefficients
        }
        if (is.null(dim(modelmat[w,]))){
          exp_predlin_x <- exp(-sum(t(glm_coef_x*t(modelmat[w, ]))))
          exp_predlin_nox <- exp(-sum(t(glm_coef_nox*t(modelmat[w,-2]))))
        }
        else{
          exp_predlin_x <- exp(-apply(t(glm_coef_x*t(modelmat[w, ])),1,sum))
          exp_predlin_nox <- exp(-apply(t(glm_coef_nox*t(modelmat[w,-2])),1,sum))
        }
        ccdf_x[[i]] <- 1/(1+exp_predlin_x) # exp_predlin_x/(1+exp_predlin_x)
        ccdf_nox[[i]] <- 1/(1+exp_predlin_nox) # exp_predlin_nox/(1+exp_predlin_nox)

      }
    }

    ccdf_x <- unlist(ccdf_x, use.names = FALSE)
    ccdf_nox <- unlist(ccdf_nox, use.names = FALSE)
    cdf <- unlist(cdf, use.names = FALSE)

    return(list(cdf=cdf, ccdf_nox=ccdf_nox, ccdf_x=ccdf_x, y=y_sort, x_sort=x_sort, z=z_sort))
  }
}

