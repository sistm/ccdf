library(purrr)
library(MAST)
library(scDD)
library(aod)
library(arm)
library(fdrtool)
library(lars)
library(emdist)
library(scran)
library(scater)
library(SingleCellExperiment)
library(doParallel)
library(ccdf)


######## Indicators ##########


truth <- c(rep("DE",250),rep("DM",250),rep("DP",250),rep("DB",250),rep("EE",4500),rep("EB",4500))
truth_bin <- c(rep(1,1000),rep(0,9000))


FDR <- function(x){
  if (sum(x<0.05)==0){
    res <- 0
  }
  if (sum(x<0.05)!=0){
    res <- sum(x[which(truth_bin==0)]<0.05)/sum(x<0.05)
  }
  
  return(res)
}


type1_error <- function(x){
  res <- sum(x[which(truth_bin==0)]<0.05)/length(which(truth_bin==0))
  return(res)
}

TDR <- function(x){
  if (sum(x<0.05)==0){
    res <- 0
  }
  if (sum(x<0.05)!=0){
    res <- sum(x[which(truth_bin==1)]<0.05)/sum(x<0.05)
  }
  
  return(res)
}

stat_power <- function(x){
  res <- sum(x[which(truth_bin==1)]<0.05)/length(which(truth_bin==1))
  return(res)
}



res_DE <- function(res_bin){
  return(table(factor(truth_bin[truth=="DE"],levels=c("0","1")),
               factor(res_bin[truth=="DE"],levels=c("0","1")))[2,2]/
           (table(factor(truth_bin[truth=="DE"],levels=c("0","1")),
                  factor(res_bin[truth=="DE"],levels=c("0","1")))[2,2]
            +table(factor(truth_bin[truth=="DE"],levels=c("0","1")),
                   factor(res_bin[truth=="DE"],levels=c("0","1")))[2,1]))
}

res_DP <- function(res_bin){
  return(table(factor(truth_bin[truth=="DP"],levels=c("0","1")),
               factor(res_bin[truth=="DP"],levels=c("0","1")))[2,2]/
           (table(factor(truth_bin[truth=="DP"],levels=c("0","1")),
                  factor(res_bin[truth=="DP"],levels=c("0","1")))[2,2]
            +table(factor(truth_bin[truth=="DP"],levels=c("0","1")),
                   factor(res_bin[truth=="DP"],levels=c("0","1")))[2,1]))
}

res_DM <- function(res_bin){
  return(table(factor(truth_bin[truth=="DM"],levels=c("0","1")),
               factor(res_bin[truth=="DM"],levels=c("0","1")))[2,2]/
           (table(factor(truth_bin[truth=="DM"],levels=c("0","1")),
                  factor(res_bin[truth=="DM"],levels=c("0","1")))[2,2]
            +table(factor(truth_bin[truth=="DM"],levels=c("0","1")),
                   factor(res_bin[truth=="DM"],levels=c("0","1")))[2,1]))
}

res_DB <-  function(res_bin){
  return(table(factor(truth_bin[truth=="DB"],levels=c("0","1")),
               factor(res_bin[truth=="DB"],levels=c("0","1")))[2,2]/
           (table(factor(truth_bin[truth=="DB"],levels=c("0","1")),
                  factor(res_bin[truth=="DB"],levels=c("0","1")))[2,2]
            +table(factor(truth_bin[truth=="DB"],levels=c("0","1")),
                   factor(res_bin[truth=="DB"],levels=c("0","1")))[2,1]))
}

res_EE <- function(res_bin){
  return(table(factor(truth_bin[truth=="EE"],levels=c("0","1")),
               factor(res_bin[truth=="EE"],levels=c("0","1")))[1,1]/
           (table(factor(truth_bin[truth=="EE"],levels=c("0","1")),
                  factor(res_bin[truth=="EE"],levels=c("0","1")))[1,1]
            +table(factor(truth_bin[truth=="EE"],levels=c("0","1")),
                   factor(res_bin[truth=="EE"],levels=c("0","1")))[1,2]))
}

res_EB <- function(res_bin){
  return(table(factor(truth_bin[truth=="EB"],levels=c("0","1")),
               factor(res_bin[truth=="EB"],levels=c("0","1")))[1,1]/
           (table(factor(truth_bin[truth=="EB"],levels=c("0","1")),
                  factor(res_bin[truth=="EB"],levels=c("0","1")))[1,1]
            +table(factor(truth_bin[truth=="EB"],levels=c("0","1")),
                   factor(res_bin[truth=="EB"],levels=c("0","1")))[1,2]))
}


########## Matrix ##########


sample_mat_NBP <- function(n_G,n){
  
  Y <- matrix(0,n_G,n)
  
  #p <- rdunif(n_G,floor(n/8),floor(n/4))
  p <- floor(rep(n/8,n_G))
  #p1 <- rdunif(n_G,floor(n/8),floor(n/4))
  #p1[which(p1%%2==1)] <- p1[which(p1%%2==1)]+1
  #p2 <- rdunif(n_G,n/4,floor(n/2)-2)
  #p2 <- floor(0.5*p1)
  p1 <- floor(rep(0.1*(n/4),n_G))
  p3 <- floor(rep(0.3*(n/4),n_G))
  p2 <- (n/4)-p1
  p4 <- (n/4)-p3
  
  moy0 <- runif(n_G,0,0.5)
  moy1 <- rdunif(n_G,10,20)
  moy2 <- 2*moy1
  moy3 <- (moy1+moy2)/2 # 2*moy1
  moy4 <- 3*moy1
  
  
  X <- c(rep(1,n/4),rep(2,n/4),rep(3,n/4),rep(4,n/4))
  
  bimod <- rep(0,n/4)
  bimod1 <- rep(0,n/4)
  bimod2 <- rep(0,n/4)
  unimod <- rep(0,n/4)
  
  for (i in 1:n_G){
    
    if (i<=250){ # DE
      
      Y[i,] <- c(rnbinom(round(n/40),prob=0.5,size=moy0[i]),rnbinom((n/4)-round(n/40),prob=0.5,size=moy1[i]),
                 rnbinom(round(n/40),prob=0.5,size=moy0[i]),rnbinom((n/4)-round(n/40),prob=0.5,size=moy3[i]),
                 rnbinom(round(n/40),prob=0.5,size=moy0[i]),rnbinom((n/4)-round(n/40),prob=0.5,size=moy2[i]),
                 rnbinom(round(n/40),prob=0.5,size=moy0[i]),rnbinom((n/4)-round(n/40),prob=0.5,size=moy4[i]))
      
    }
    
    if (i<=500 & i>250){ # DM
      
      bimod1 <- rep(0,n/4)
      unimod1 <- rep(0,n/4)
      trimod1 <- rep(0,n/4)
      unimod2 <- rep(0,n/4)
      
      bimod1[1:p[i]] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/40),prob=0.5,size=moy1[i]))
      bimod1[(p[i]+1):(n/4)] <- rnbinom(p[i],prob=0.5,size=moy2[i])
      
      unimod1 <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(9*n/40,prob=0.5,size=moy2[i]))
      
      trimod1[1:round((n/4)/3)] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(length(1:round((n/4)/3))-floor(n/40),prob=0.5,size=moy1[i]))
      trimod1[round(((n/4)/3)+1):round(2*(n/4)/3)] <- rnbinom(length(round(((n/4)/3)+1):round(2*(n/4)/3)),prob=0.5,size=moy2[i])
      trimod1[round((2*(n/4)/3)+1):(n/4)] <- rnbinom(length(round((2*(n/4)/3)+1):(n/4)),prob=0.5,size=moy3[i])
      
      unimod2 <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(9*n/40,prob=0.5,size=moy1[i]))

      
      Y[i,] <- c(unimod1,bimod1,trimod1,unimod2)
      
    }
    
    if (i<=750 & i>500){ # DP
      
      bimod1 <- rep(0,n/4)
      bimod2 <- rep(0,n/4)
      bimod3 <- rep(0,n/4)
      bimod4 <- rep(0,n/4)
      
      bimod1[1:p1[i]] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(p1[i]-floor(n/40),prob=0.5,size=moy1[i]))
      bimod1[(p1[i]+1):(n/4)] <- rnbinom(p2[i],prob=0.5,size=moy2[i])
      
      bimod2[1:p2[i]] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(p2[i]-floor(n/40),prob=0.5,size=moy1[i]))
      bimod2[(p2[i]+1):(n/4)] <- rnbinom(p1[i],prob=0.5,size=moy2[i])
      
      bimod3[1:p3[i]] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(p3[i]-floor(n/40),prob=0.5,size=moy1[i]))
      bimod3[(p3[i]+1):(n/4)] <- rnbinom(p4[i],prob=0.5,size=moy2[i])
      
      bimod4[1:p4[i]] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(p4[i]-floor(n/40),prob=0.5,size=moy1[i]))
      bimod4[(p4[i]+1):(n/4)] <- rnbinom(p3[i],prob=0.5,size=moy2[i])
      
      Y[i,] <- c(bimod1,bimod2,bimod3,bimod4)
      
      
    }
    
    if (i<=1000 & i>750){ # DB
      
      unimod1 <- rep(0,n/4)
      bimod1 <- rep(0,n/4)
      trimod1 <- rep(0,n/4)
      quadmod1 <- rep(0,n/4)
      
      delta <- 0.9*moy3[i]
      gamma <- 0.6*moy3[i]
      epsilon <- 0.3*moy3[i]
      
      # a faire
      
      unimod1 <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(floor(9*n/40),prob=0.5,size=moy3[i])) # mean
      #unimod2 <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(floor(9*n/40),prob=0.5,size=moy3[i])) # mean

      long <- round(((n/4)-floor(n/40))/2)

      bimod1[1:((n/40)+long)] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(long,prob=0.5,size=moy3[i]-delta))
      bimod1[round((n/40)+long+1):(n/4)] <- rnbinom(length(round((n/40)+long+1):(n/4)),prob=0.5,size=moy3[i]+delta)

      long <- round(((n/4)-floor(n/40))/3)

      trimod1[1:((n/40)+long)] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(long,prob=0.5,size=moy3[i]-gamma))
      trimod1[round((n/40)+long+1):round(((n/40)+long+1)+long-1)] <- rnbinom(long,prob=0.5,size=moy3[i]+gamma)
      trimod1[round(((n/40)+long+1)+long):(n/4)] <- rnbinom(long,prob=0.5,size=moy3[i])

      long <- round(((n/4)-floor(n/40))/4)

      quadmod1[1:((n/40)+long)] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),
                                     rnbinom(long,prob=0.5,size=moy3[i]-(2*epsilon)))
      quadmod1[round((n/40)+long+1):round(((n/40)+long+1)+long-1)] <- rnbinom(long,prob=0.5,size=moy3[i]-epsilon)
      quadmod1[round(((n/40)+long+1)+long):(round(((n/40)+long+1)+long)+long-1)] <- rnbinom(long,prob=0.5,size=moy3[i]+epsilon)
      quadmod1[(round(((n/40)+long+1)+long)+long):(n/4)] <- rnbinom(length((round(((n/40)+long+1)+long)+long):(n/4)),prob=0.5,size=moy3[i]+(2*epsilon))
      
      Y[i,] <- c(unimod1,bimod1,trimod1,quadmod1)
      
    }
    
    if (i<=5500 & i>1000){
      
      Y[i,] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(floor(9*n/40),prob=0.5,size=moy1[i]),
                 rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(floor(9*n/40),prob=0.5,size=moy1[i]),
                 rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(floor(9*n/40),prob=0.5,size=moy1[i]),
                 rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(floor(9*n/40),prob=0.5,size=moy1[i]))
      
    }
    
    if (i>5500){
      
      bimod1 <- rep(0,n/4)
      bimod2 <- rep(0,n/4)
      bimod3 <- rep(0,n/4)
      bimod4 <- rep(0,n/4)
      
      bimod1[1:p[i]] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/40),prob=0.5,size=moy1[i]))
      bimod1[(p[i]+1):(n/4)] <- rnbinom(p[i],prob=0.5,size=moy2[i])
      
      bimod2[1:p[i]] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/40),prob=0.5,size=moy1[i]))
      bimod2[(p[i]+1):(n/4)] <- rnbinom(p[i],prob=0.5,size=moy2[i])
      
      bimod3[1:p[i]] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/40),prob=0.5,size=moy1[i]))
      bimod3[(p[i]+1):(n/4)] <- rnbinom(p[i],prob=0.5,size=moy2[i])
      
      bimod4[1:p[i]] <- c(rnbinom(floor(n/40),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/40),prob=0.5,size=moy1[i]))
      bimod4[(p[i]+1):(n/4)] <- rnbinom(p[i],prob=0.5,size=moy2[i])
      
      Y[i,] <- c(bimod1,bimod2,bimod3,bimod4)
      
    }
  }
  
  for (c in 1:n){
    Y[which(Y[,c]!=0),c] <- Y[which(Y[,c]!=0),c] + rnorm(length(which(Y[,c]!=0)),0,0.01)
  }
  
  return(list(Y=Y,X=X))
}



######### Results #########



size <- c(40,80,120,160,240,320,400)
fdrs <- matrix(NA,length(size),1)
type1 <- matrix(NA,length(size),1)
tdrs <- matrix(NA,length(size),1)
pwr <- matrix(NA,length(size),1)

DE <- matrix(NA,length(size),1)
DM <- matrix(NA,length(size),1)
DP <- matrix(NA,length(size),1)
DB <- matrix(NA,length(size),1)
EE <- matrix(NA,length(size),1)
EB <- matrix(NA,length(size),1)


for (n in 1:length(size)){
  print(paste0("n=", n))
  # Genes Matrix: 10000 genes of which 1000 are DE
  if (n==1){
    temp <- sample_mat_NBP(10000,40)
  }
  if (n==2){
    temp <- sample_mat_NBP(10000,80)
  }
  if (n==3){
    temp <- sample_mat_NBP(10000,120)
  }
  if (n==4){
    temp <- sample_mat_NBP(10000,160)
  }
  if (n==5){
    temp <- sample_mat_NBP(10000,240)
  }
  if (n==6){
    temp <- sample_mat_NBP(10000,320)
  }
  if (n==7){
    temp <- sample_mat_NBP(10000,400)
  }
  
  X <- temp$X
  
  genes_matrix <- temp$Y
  #Y <- calculateTPM(genes_matrix)
  Y <- genes_matrix
  
  # CCDF asymptotic test
  
  res_ccdf_asymp <- ccdf_testing(data.frame(Y=Y), data.frame(X=as.factor(X)), test="asymptotic",
                                 n_cpus = 16, space_y = TRUE, number_y = (ncol(Y)/2))$pvals$raw_pval
  
  
  # p-values data frame
  
  pvs <- data.frame(CCDF_asymp=res_ccdf_asymp)
  
  res_bin <- data.frame(CCDF_asymp=ifelse(p.adjust(pvs$CCDF_asymp,"BH")<0.05,1,0))
  
  
  DE[n,] <- apply(res_bin, 2, function(x){res_DE(x)})
  
  DM[n,] <- apply(res_bin, 2, function(x){res_DM(x)})
  
  DP[n,] <- apply(res_bin, 2, function(x){res_DP(x)})
  
  DB[n,] <- apply(res_bin, 2, function(x){res_DB(x)})
  
  EE[n,] <- apply(res_bin, 2, function(x){res_EE(x)})
  
  EB[n,] <- apply(res_bin, 2, function(x){res_EB(x)})
  
  
  # FDR
  fdrs[n,] <- apply(apply(pvs, 2, p.adjust, "BH"), 2, function(x){FDR(x)})
  # Type-1 error
  type1[n,] <- apply(pvs, 2, function(x){type1_error(x)})
  # TDR
  tdrs[n,] <- apply(apply(pvs, 2, p.adjust, "BH"), 2, function(x){TDR(x)})
  # Power
  pwr[n,] <- apply(pvs, 2, function(x){stat_power(x)})
  
}

res <- data.frame(setting=rep("4cond_NB",length(size)),n=rep(size,1),
                  fdr=matrix(fdrs,ncol=1),
                  ti_err= matrix(type1,ncol=1),
                  tpr=matrix(tdrs,ncol=1),
                  pwr= matrix(pwr,ncol=1),
                  DE= matrix(DE,ncol=1),
                  DM= matrix(DM,ncol=1),
                  DP= matrix(DP,ncol=1),
                  DB= matrix(DB,ncol=1),
                  EE= matrix(EE,ncol=1),
                  EB= matrix(EB,ncol=1))
