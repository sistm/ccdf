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
  p <- floor(rep(n/4,n_G))
  #p1 <- rdunif(n_G,floor(n/8),floor(n/4))
  #p1[which(p1%%2==1)] <- p1[which(p1%%2==1)]+1
  #p2 <- rdunif(n_G,n/4,floor(n/2)-2)
  #p2 <- floor(0.5*p1)
  p1 <- floor(rep(n/6,n_G))
  p2 <- (n/2)-p1
  
  moy0 <- runif(n_G,0,0.5)
  moy1 <- rdunif(n_G,10,20) # 1 10
  #moy2 <- rdunif(n_G,15,25) # 5 15
  #moy3 <- rdunif(n_G,10,20) # 5 10
  
  moy2 <- 3*moy1
  moy3 <- (moy1+moy2)/2
  
  X <- c(rep(1,n/2),rep(2,n/2))
  
  bimod <- rep(0,n/2)
  bimod1 <- rep(0,n/2)
  bimod2 <- rep(0,n/2)
  unimod <- rep(0,n/2)
  
  for (i in 1:n_G){
    
    if (i<=250){ # DE
      
      Y[i,] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(floor(9*n/20),prob=0.5,size=moy1[i]),
                 rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(floor(9*n/20),prob=0.5,size=moy3[i]))
      
    }
    
    if (i<=500 & i>250){ # DM
      
      bimod <- rep(0,n/2)
      unimod <- rep(0,n/2)
      
      bimod[1:p[i]] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/20),prob=0.5,size=moy1[i]))
      bimod[(p[i]+1):(n/2)] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/20),prob=0.5,size=moy2[i]))
      unimod <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(9*n/20,prob=0.5,size=moy2[i]))
      Y[i,] <- c(unimod,bimod)
      
    }
    
    if (i<=750 & i>500){ # DP
      
      bimod1 <- rep(0,n/2)
      bimod2 <- rep(0,n/2)
      
      bimod1[1:p1[i]] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p1[i]-floor(n/20),prob=0.5,size=moy1[i]))
      bimod1[(p1[i]+1):(n/2)] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p2[i]-floor(n/20),prob=0.5,size=moy2[i]))
      
      bimod2[1:p2[i]] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p2[i]-floor(n/20),prob=0.5,size=moy1[i]))
      bimod2[(p2[i]+1):(n/2)] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p1[i]-floor(n/20),prob=0.5,size=moy2[i]))
      
      Y[i,] <- c(bimod1,bimod2)
      
      
    }
    
    if (i<=2000 & i>750){ # DB
      
      bimod <- rep(0,n/2)
      unimod <- rep(0,n/2)
      
      bimod[1:p[i]] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/20),prob=0.5,size=moy1[i]))
      bimod[(p[i]+1):(n/2)] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/20),prob=0.5,size=moy2[i]))
      unimod <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(floor(9*n/20),prob=0.5,size=moy3[i]))
      Y[i,] <- c(unimod,bimod)
      
    }
    
    if (i<=5500 & i>1000){
      
      Y[i,] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(floor(9*n/20),prob=0.5,size=moy1[i]),
                 rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(floor(9*n/20),prob=0.5,size=moy1[i]))
      
    }
    
    if (i>5500){
      
      bimod1 <- rep(0,n/2)
      bimod2 <- rep(0,n/2)
      
      bimod1[1:p[i]] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/20),prob=0.5,size=moy1[i]))
      bimod1[(p[i]+1):(n/2)] <-   c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/20),prob=0.5,size=moy2[i]))
      
      bimod2[1:p[i]] <- c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/20),prob=0.5,size=moy1[i]))
      bimod2[(p[i]+1):(n/2)] <-  c(rnbinom(floor(n/20),prob=0.5,size=moy0[i]),rnbinom(p[i]-floor(n/20),prob=0.5,size=moy2[i]))
      
      Y[i,] <- c(bimod1,bimod2)
      
    }
  }
  
  for (c in 1:n){
    Y[which(Y[,c]!=0),c] <- Y[which(Y[,c]!=0),c]  + rnorm(length(which(Y[,c]!=0)),0,0.01)
  }
  
  return(list(Y=Y,X=X))
}



######### Results #########

size <- c(20,40,60,80,100,160,200)
fdrs <- matrix(NA,length(size),3)
type1 <- matrix(NA,length(size),3)
tdrs <- matrix(NA,length(size),3)
pwr <- matrix(NA,length(size),3)

DE <- matrix(NA,length(size),3)
DM <- matrix(NA,length(size),3)
DP <- matrix(NA,length(size),3)
DB <- matrix(NA,length(size),3)
EE <- matrix(NA,length(size),3)
EB <- matrix(NA,length(size),3)


for (n in 1:length(size)){
  print(n)
  # Genes Matrix: 10000 genes of which 1000 are DE
  if (n==1){
    temp <- sample_mat_NBP(10000,20)
  }
  if (n==2){
    temp <- sample_mat_NBP(10000,40)
  }
  if (n==3){
    temp <- sample_mat_NBP(10000,60)
  }
  if (n==4){
    temp <- sample_mat_NBP(10000,80)
  }
  if (n==5){
    temp <- sample_mat_NBP(10000,100)
  }
  if (n==6){
    temp <- sample_mat_NBP(10000,160)
  }
  if (n==7){
    temp <- sample_mat_NBP(10000,200)
  }

  X <- temp$X

  genes_matrix <- temp$Y
  Y <- genes_matrix


  sce <- SingleCellExperiment(list(counts=Y))
  ngenes <- nrow(genes_matrix)
  rownames(sce) <- seq_len(ngenes)
  assayNames(sce) <- "counts"
  colData(sce)$condition <- X

  # MAST

  ncells <- size[n]
  fData <- data.frame(primerid=seq_len(ngenes))
  cData <- data.frame(wellKey=seq_len(ncells))
  sca <- FromMatrix(as.matrix(sce@assays@data$counts), cData, fData,check_sanity = FALSE)
  assayNames(sca) <- "TPM"
  #
  cond <- factor(X)
  cond <- relevel(cond,"1")
  colData(sca)$condition <- cond
  zlmCond <- zlm(~condition, sca)
  #
  # only test the condition coefficient.
  summaryCond <- summary(zlmCond,doLRT="condition2")
  summaryDt <- summaryCond$datatable
  fcHurdle <- summaryDt[contrast=='condition2' & component=='H',.(`primerid`,`Pr(>Chisq)`)] #hurdle P values
  res_MAST <- fcHurdle$`Pr(>Chisq)`[sort(as.numeric(fcHurdle$primerid),index.return=TRUE)$ix]

  # SigEMD

  source("FunImpute.R")
  source("SigEMDHur.R")
  source("SigEMDnonHur.R")
  source("plot_sig.R")

  Y <- sce@assays@data$counts
  rownames(Y) <- seq_len(ngenes)
  colnames(Y) <- seq_len(ncells)
  names(X) <- colnames(Y)

  databinary<- databin(Y)
  Hur_gene<- idfyImpgene(Y,databinary,X)
  # genes_use<- idfyUsegene(Y,databinary,X,ratio = 0.5)
  #
  # datimp <- FunImpute(object = Y, genes_use = (genes_use), genes_fit = (Hur_gene),dcorgene = NULL)
  # Y<-datimp$alldat

  results<- calculate_single(data =  Y, condition =  X, Hur_gene = Hur_gene, nperm=500, binSize=0.2)
  emd <- results$emdall
  res_SigEMD <- emd[,2]

  # scDD

  sca <- sce
  rownames(sca) <- seq_len(ngenes)
  assayNames(sca) <- "normcounts"
  #colData(sce)$condition <- X
  #prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.0)
  scDatExSim <- scDD(sca, prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01),testZeroes=FALSE, categorize = FALSE)
  RES <- results(scDatExSim)
  res_scDD <- RES$nonzero.pvalue

  # CCDF asymptotic test

  res_ccdf_asymp <- ccdf_testing(data.frame(Y=Y), data.frame(X=as.factor(X)), test="asymptotic", n_cpus=16, 
                                 space_y = TRUE, number_y = (ncol(Y)/2))$pvals$raw_pval


  # p-values data frame
  
  pvs <- data.frame(MAST=res_MAST,
                    scDD=res_scDD,
                    SigEMD=res_SigEMD,
                    CCDF_asymp=res_ccdf_asymp)

  res_bin <- data.frame(MAST=ifelse(p.adjust(pvs$MAST,"BH")<0.05,1,0),
                        scDD=ifelse(p.adjust(pvs$scDD,"BH")<0.05,1,0),
                        SigEMD=ifelse(p.adjust(pvs$SigEMD,"BH")<0.05,1,0),
                        CCDF_asymp=ifelse(p.adjust(pvs$CCDF_asymp,"BH")<0.05,1,0))


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


res <- data.frame(setting=rep("2cond_NB",28),n=rep(size,4),
                  method=factor(as.character(rep(c("MAST","scDD","SigEMD","CCDF_asymp"),each=7))),
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

