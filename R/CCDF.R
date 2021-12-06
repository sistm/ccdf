#' Function to compute (un)conditional cumulative distribution function (CDF), used by plot_CCDF function.
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
#'@param sample_group a vector of length \code{n} indicating whether the samples
#'should be grouped (e.g. paired samples or longitudinal data). Coerced
#'to be a \code{factor}. Default is \code{NULL} in which case no grouping is
#'performed.
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
#'@param space_y a logical flag indicating whether the y thresholds are spaced.
#'When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and
#'the maximum of the observations is used. Default is \code{FALSE}.
#'
#'@param number_y an integer value indicating the number of y thresholds (and therefore
#'the number of regressions) to perform the test. Default is \code{length(Y)}.
#'
#'@importFrom stats model.matrix
#' @importFrom lme4 lmer
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
#'X <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
#'Y <- ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,0.5,1))
#'res <- CCDF(Y,data.frame(X=X),method="linear regression")
#'
#' set.seed(1)
#' Nsample_groups = 10000
#' Npersample_group = 10
#' N = Nsample_groups * Npersample_group
#' sample_group = factor(rep(1:Nsample_groups, each = Npersample_group))
#' u = rnorm(Nsample_groups, sd = .5)
#' e = rnorm(N, sd = .25)
#' X = rnorm(N)
#' Y = 2 + .5 * X + u[sample_group] + e
#'


CCDF <- function(Y,X,Z=NULL,method=c("linear regression","logistic regression","RF", "mixed model"), fast=TRUE, space_y=FALSE, number_y=length(Y),
                 sample_group=NULL){

  if (is.data.frame(Y)) Y <- Y[,1]
  if (space_y){
   y <- seq(min(Y[-which(Y==min(Y))]),max(Y),length.out=number_y)
  }
  else{
    y <- sort(unique(Y))
  }

  if (class(Z)=="NULL"){
    n_Y <- length(Y)

    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    modelmat <- model.matrix(~.,data=X)

    cdf <- list()
    ccdf <- list()
    x_sort <- NULL
    y_sort <- NULL
    sample_group_sort <- NULL

    for (i in 1:(length(y)-1)){
      if (i==1){
        w <- which(Y<=y[i])
      }
      else{
        w <- which((Y<=y[i])&(Y>y[i-1]))
      }
      x_sort <- c(x_sort,X[w,])
      y_sort <- c(y_sort,Y[w])
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
          reg <- lm(indi_Y ~ 1 + modelmat[,-1])
          ccdf[[i]] <- predict(reg)[w]
        }
        else if(method == "mixed model"){
          if(is.null(sample_group)) {
            warning("sample_group is null",
                "No random effects terms specified in formula")
          }
          else{
            mod_mixed <- lmer(indi_Y ~ 1 + modelmat[,-1] + (1|sample_group))
            ccdf[[i]] <- fitted(mod_mixed)[w]
            sample_group_sort <- c(sample_group_sort, sample_group[w])
            if (is.factor(sample_group)) sample_group_sort <- as.factor(sample_group_sort)
          }
        }
      }
    }

    if (is.factor(X[,1])){
      x_sort <- as.factor(x_sort)
    }
    ccdf <- unlist(ccdf, use.names = FALSE)
    cdf <- unlist(cdf, use.names = FALSE)
    ifelse(
      method == "mixed model",
           return(list(cdf=cdf, ccdf=ccdf, y=y_sort, x=x_sort, sample_group = sample_group_sort)),
                  return(list(cdf=cdf, ccdf=ccdf, y=y_sort, x=x_sort))
    )
  }

  else{

    n_Y <- length(Y)

    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    colnames(Z) <- sapply(1:ncol(Z), function(i){paste0('Z',i)})
    modelmat <- model.matrix(~.,data=cbind(X,Z))

    ind_X <- which(substring(colnames(modelmat),1,1)=="X")

    x_sort <- NULL
    y_sort <- NULL
    z_sort <- NULL
    sample_group_sort <- NULL

    ccdf_x <- list()
    ccdf_nox <- list()
    cdf <- list()

    for (i in 1:(length(y)-1)){

      if (i==1){
        w <- which(Y<=y[i])
      }
      else{
        w <- which((Y<=y[i])&(Y>y[i-1]))
      }

      x_sort <- c(x_sort,X[w,])
      y_sort <- c(y_sort,Y[w])
      z_sort <- c(z_sort,Z[w,])
      indi_Y <- 1*(Y<=y[i])

      # unCDF
      cdf[[i]] <- rep(sum(indi_Y)/n_Y, length(w))

      if (length(unique(indi_Y))==1){
        ccdf[[i]] <- rep(1,length(w))
      }

      if (method=="RF"){
        # CDF
        rf_x <- randomForest(data.frame(X=X, Z=Z, row.names = NULL), indi_Y,ntree=50)
        ccdf_x[[i]] <- predict(rf_x)[w]
        # CCDF
        rf_nox <- randomForest(data.frame(Z=Z, row.names = NULL), indi_Y,ntree=50)
        ccdf_nox[[i]] <-  predict(rf_nox)[w]

      }


      else if (method=="logistic regression"){
        if (fast){
          glm_coef_x <- RcppNumerical::fastLR(x=modelmat, y=indi_Y)$coefficients
          glm_coef_nox <- RcppNumerical::fastLR(x= modelmat[,-ind_X], y=indi_Y)$coefficients
        }
        else{
          glm_coef_x <- glm.fit(x=modelmat,y=indi_Y, family = binomial())$coefficients
          glm_coef_nox <- glm.fit(x=modelmat[,-ind_X],y=indi_Y, family = binomial())$coefficients
        }
        if (is.null(dim(modelmat[w,]))){
          exp_predlin_x <- exp(-sum(t(glm_coef_x*t(modelmat[w, ]))))
          exp_predlin_nox <- exp(-sum(t(glm_coef_nox*t(modelmat[w,-ind_X]))))
        }
        else{
          exp_predlin_x <- exp(-apply(t(glm_coef_x*t(modelmat[w, ])),1,sum))
          exp_predlin_nox <- exp(-apply(t(glm_coef_nox*t(modelmat[w,-ind_X])),1,sum))
        }
        ccdf_x[[i]] <- 1/(1+exp_predlin_x) # exp_predlin_x/(1+exp_predlin_x)
        ccdf_nox[[i]] <- 1/(1+exp_predlin_nox) # exp_predlin_nox/(1+exp_predlin_nox)

      }

      else if (method=="linear regression"){
        reg_x <- lm(indi_Y ~ 1 + modelmat[,-1])
        ccdf_x[[i]] <- predict(reg_x)[w]
        reg_nox <- lm(indi_Y ~ 1 + modelmat[,-c(1,ind_X)])
        ccdf_nox[[i]] <- predict(reg_nox)[w]
      }
      else if(method == "mixed model"){
        if(is.null(sample_group)) {
            warning("sample_group is null",
                "No random effects terms specified in formula")
          }
        else{
          mod_mixed_x <- lmer(indi_Y ~ 1 + modelmat[,-1] + (1|sample_group))
          ccdf_x[[i]] <- fitted(mod_mixed_x)[w]
          mod_mixed_nox <- lmer(indi_Y ~ 1 + modelmat[,-c(1,ind_X)] + (1|sample_group))
          ccdf_nox[[i]] <- fitted(mod_mixed_nox)[w]
          sample_group_sort <- c(sample_group_sort, sample_group[w])
          if(is.factor(sample_group)) sample_group_sort <- as.factor(sample_group_sort)
        }

        }
    }

    if (is.factor(X[,1])){
      x_sort <- as.factor(x_sort)
    }
    if (is.factor(Z[,1])){
      z_sort <- as.factor(z_sort)
    }
    ccdf_x <- unlist(ccdf_x, use.names = FALSE)
    ccdf_nox <- unlist(ccdf_nox, use.names = FALSE)
    cdf <- unlist(cdf, use.names = FALSE)
    ifelse(
      method == "mixed model",
           return(list(cdf=cdf, ccdf_nox=ccdf_nox, ccdf_x=ccdf_x, y=y_sort, x=x_sort, z=z_sort, sample_group = sample_group_sort)),
                  return(list(cdf=cdf, ccdf_nox=ccdf_nox, ccdf_x=ccdf_x, y=y_sort, x=x_sort, z=z_sort))
    )

  }
}

