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
library(DEsingle)
library(BiocParallel)


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


sample_mat_Z <- function(n_G,n){
  
  Y <- matrix(0,n_G,n)
  Z <- rnorm(n,10,2)
  X <- rep(0,n)
  X[which(Z<quantile(Z)[2])] <- rep(1,length(which(Z<quantile(Z)[2])))
  X[which((Z>quantile(Z)[3])&(Z<quantile(Z)[4]))] <- rep(1,length(which((Z>quantile(Z)[3])&(Z<quantile(Z)[4]))))
  X <- X+1
  
  for (i in 1:n_G){
    
    if (i<=1000){
      Y[i,] <- rnorm(n,5,1)*(X==1) + rnorm(n,(5+rnorm(1,1,0.5)),1)*(X==2) + rnorm(n,0,1) # H1
    }
    else{
      #Y[i,seq((n/20)+1,(n-(n/20)))] <-  rnorm(n-(n/10),10,1)*Z[seq((n/20)+1,(n-(n/20)))] + rnorm(n-(n/10),0,1) # H0
      Y[i,] <-  runif(1,0.3,0.5)*Z + rnorm(n,0,0.5) # H0
    }
  }
  
  for (j in 1:n_G){
    Y[j,which(Y[j,]<0)] <- rep(0,length(which(Y[j,]<0)))
  }
  
  return(list(Y=Y,X=X,Z=Z))
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
    temp <- sample_mat_Z(10000,20)
  }
  if (n==2){
    temp <- sample_mat_Z(10000,40)
  }
  if (n==3){
    temp <- sample_mat_Z(10000,60)
  }
  if (n==4){
    temp <- sample_mat_Z(10000,80)
  }
  if (n==5){
    temp <- sample_mat_Z(10000,100)
  }
  if (n==6){
    temp <- sample_mat_Z(10000,160)
  }
  if (n==7){
    temp <- sample_mat_Z(10000,200)
  }
  
  X <- temp$X
  Z <- temp$Z
  
  genes_matrix <- temp$Y
  Y <- genes_matrix
  rownames(Y) <- as.character(seq_len(10000))
  colnames(Y) <- as.character(seq_len(size[n]))
  
  
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
  colData(sca)$Z <- Z
  zlmCond <- zlm(~condition, sca)
  
  #only test the condition coefficient.
  summaryCond <- summary(zlmCond,doLRT="condition2")
  summaryDt <- summaryCond$datatable
  fcHurdle <- summaryDt[contrast=='condition2' & component=='H',.(`primerid`,`Pr(>Chisq)`)] #hurdle P values
  res_MAST <- fcHurdle$`Pr(>Chisq)`[sort(as.numeric(fcHurdle$primerid),index.return=TRUE)$ix]
  
  # SigEMD
  
  source("FunImpute.R")
  source("SigEMDHur.R")
  source("SigEMDnonHur.R")
  source("plot_sig.R")
  # 
  Y <- sce@assays@data$counts
  rownames(Y) <- seq_len(ngenes)
  colnames(Y) <- seq_len(ncells)
  names(X) <- colnames(Y)
  # 
  databinary<- databin(Y)
  # 
  Hur_gene<- idfyImpgene(Y,databinary,X)
  #genes_use<- idfyUsegene(Y,databinary,X,ratio = 0.5)
  # 
  #datimp <- FunImpute(object = Y, genes_use = (genes_use), genes_fit = (Hur_gene),dcorgene = NULL)
  #Y<-datimp$alldat
  # 
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
  
  res_ccdf_asymp <- cit(data.frame(Y=Y), data.frame(X=as.factor(X)), data.frame(Z=as.numeric(Z)), test="asymptotic", n_cpus = 16,
                                 space_y = TRUE, number_y = (ncol(Y)/2))$pvals
  
  # p-values data frame
  
  pvs <- data.frame(MAST=res_MAST,
    scDD=res_scDD,
    SigEMD=res_SigEMD,
    CCDF_asymp=res_ccdf_asymp)
  
  
  res_bin <- data.frame(MAST=ifelse(p.adjust(pvs$MAST,"BH")<0.05,1,0),
    scDD=ifelse(p.adjust(pvs$scDD,"BH")<0.05,1,0),
    SigEMD=ifelse(p.adjust(pvs$SigEMD,"BH")<0.05,1,0),
    CCDF_asymp=ifelse(p.adjust(pvs$CCDF_asymp,"BH")<0.05,1,0))

  
  # FDR
  fdrs[n,] <- apply(apply(pvs, 2, p.adjust, "BH"), 2, function(x){FDR(x)})
  # Type-1 error
  type1[n,] <- apply(pvs, 2, function(x){type1_error(x)})
  # TDR
  tdrs[n,] <- apply(apply(pvs, 2, p.adjust, "BH"), 2, function(x){TDR(x)})
  # Power
  pwr[n,] <- apply(pvs, 2, function(x){stat_power(x)})
  
}


res <- data.frame(setting=rep("With_Covariate",28),n=rep(size,4),
                  method=factor(as.character(rep(c("MAST","scDD","SigEMD","CCDF_asymp"),each=7))),
                  fdr=matrix(fdrs,ncol=1),
                  ti_err= matrix(type1,ncol=1),
                  tpr=matrix(tdrs,ncol=1),
                  pwr= matrix(pwr,ncol=1))
