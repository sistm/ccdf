library(ccdf)

# test_perm_Modif(Y=Y[1,],X=data.frame(X=X),
#                 n_perm = 100,
#                 number_y = 10)


test_perm_Modif <- function(Y, X, Z = NULL, n_perm = 100,
                            space_y = FALSE, number_y = length(Y)){
  
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
    reg <- lm(indi_Y ~ as.matrix(modelmat[,-1]))
    beta[i,] <- reg$coefficients[ind_X]
  }
  
  beta <- as.vector(beta)
  
  z <- sqrt(length(Y))*(beta[-length(y)])
  STAT_obs <- sum(t(z)*z)
  
  STAT_perm <- rep(NA,n_perm)
  
  if (is.null(Z)){
    for (k in 1:n_perm){
      
      index_non0 <- which(Y!=0)
      index <- 1:nrow(X)
      index[index_non0] <- sample(index_non0)
      X_star <- data.frame(X=X[index, ])
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
  stperm = matrix(STAT_perm,ncol = n_perm)
  colnames(stperm) <- paste0("STAT_perm_",1:n_perm)
  
  return(cbind.data.frame(data.frame(raw_pval=pval,
                                     score=score,
                                     STAT_obs=STAT_obs),stperm))
  
}

