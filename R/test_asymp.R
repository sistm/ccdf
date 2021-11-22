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
#' @param sample_group a vector of length \code{n} indicating whether the samples
#'should be grouped (e.g. paired samples or longitudinal data). Coerced
#'to be a \code{factor}. Default is \code{NULL} in which case no grouping is
#'performed.
#'
#'@param space_y a logical flag indicating whether the y thresholds are spaced. 
#'When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and 
#'the maximum of the observations is used. Default is \code{FALSE}.
#'
#'@param number_y an integer value indicating the number of y thresholds (and therefore
#'the number of regressions) to perform the test. Default is \code{length(Y)}.
#'
#'@importFrom lme4 lmer
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


test_asymp <- function(Y, X, Z = NULL,method=c("linear regression", "mixed model"), sample_group=NULL,space_y = FALSE, number_y = length(unique(Y))){
  
  Y <- as.numeric(Y)
  
  if (space_y){
    y <- seq(ifelse(length(which(Y==0))==0,min(Y),min(Y[-which(Y==0)])),max(Y[-which.max(Y)]),length.out=number_y)
  }
  else{
    y <- sort(unique(Y))
  }
  if(length(method)>1) method <- method[1]
  stopifnot(method %in% c("linear regression","mixed model"))
  # no covariates Z
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
  #nb_fact <- colSums(modelmat[,ind_X])
  
  beta <- matrix(NA,(length(y)-1),length(ind_X))
  indi_pi <- matrix(NA,length(Y),(length(y)-1))
  
  Phi <- (1/length(Y))*(t(as.matrix(modelmat))%*%as.matrix(modelmat))
  H <- (solve(Phi)%*%t(as.matrix(modelmat))) # ginv
  H <- H[ind_X,]

  for (i in 1:(length(y)-1)){ # on fait varier le seuil
    indi_Y <- 1*(Y<=y[i])
    indi_pi[,i] <- indi_Y
    #beta[i,] <- (1/length(indi_Y))*rowSums(sapply(1:length(Y),function(i){H[i]*indi_pi[i,]}))
    if(method == "linear regression"){
      reg <- lm(indi_Y ~ as.matrix(modelmat[,-1]))
      beta[i,] <- reg$coefficients[ind_X]
    }
    else if(method == "mixed model"){
      if(is.null(sample_group)) {
            warning("Some transcripts in the investigated gene sets were ",
                    "not measured:\nremoving those transcripts from the ",
                    "gene set definition...")
            break
          }
          else{
            mod_mixed <- lmer(indi_Y ~ 1 + modelmat[,-1] + (1|sample_group))
            beta[i,] <- lme4::fixef(mod_mixed)[ind_X]
      }
    }

  }
  
  beta <- as.vector(beta)
  prop <- colMeans(indi_pi)
  
  if (is.null(dim(H))){
    H_square <- sum(H^2)
    Sigma <- sapply(1:(length(y)-1), function(i){sapply(1:((length(y)-1)*length(ind_X)), function(j){
      if (i<=j){
        (prop[i]-(prop[j]*prop[i]))
      }
      else{
        (prop[j]-(prop[j]*prop[i]))
      }
    })})
    Sigma <- (1/length(Y))*(H_square*Sigma)
  }
  else{
    
    temp_Sigma <-  lapply(1:ncol(H), function(k){sapply(1:nrow(H), function(s){sapply(1:nrow(H), function(r){H[s,k]*H[r,k]})})})
    sum_temp_Sigma <- temp_Sigma[[1]]
    for (i in 2:ncol(H)){
      sum_temp_Sigma <- sum_temp_Sigma + temp_Sigma[[i]]
    }
    
    ind_sig <- rep(1:(length(y)-1),length(ind_X))
    
    Sigma <- matrix(NA,((length(y)-1)*length(ind_X)),((length(y)-1)*length(ind_X)))
    for (i in 1:((length(y)-1)*length(ind_X))){
      for (j in 1:((length(y)-1)*length(ind_X))){
        if (i<=j){
          Sigma[i,j] <- sum_temp_Sigma[floor(i/(length(y))+1),floor(j/(length(y))+1)]*(prop[ind_sig[i]]-(prop[ind_sig[j]]*prop[ind_sig[i]]))
        }
        else{
          Sigma[i,j] <- sum_temp_Sigma[floor(i/(length(y))+1),floor(j/(length(y))+1)]*(prop[ind_sig[j]]-(prop[ind_sig[j]]*prop[ind_sig[i]]))
        }
      }
    }
    Sigma <- (1/length(Y))*Sigma
  }
  decomp <- eigen(Sigma)
  A <- matrix(0,(length(ind_X)*(length(y)-1)),(length(ind_X)*(length(y)-1)))
  diag(A) <- decomp$values
  z <- (sqrt(length(Y)))*beta
  STAT <- sum(t(z)*z)
  pval <- survey::pchisqsum(STAT, lower.tail = FALSE, df = rep(1,length(diag(A))), a = diag(A), method = "saddlepoint")
  
  return(data.frame(raw_pval=pval,Stat=STAT))
  
}


