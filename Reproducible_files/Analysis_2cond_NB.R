library(ccdf)
library(reticulate)
library(DESeq2)
library(scDD)
library(arm)
library(scde)
library(purrr)
library(MAST)
library(aod)
library(fdrtool)
library(lars)
library(emdist)
library(scran)
library(scater)
library(SingleCellExperiment)
library(doParallel)
library(CompQuadForm)
library(EBSeq)
library(DEsingle)

load("seeds_ccdf.RData")

###########
sample_mat_NBP <- function(n_G,n){
  prop_non_0 <- 0.9 # 0.1
  Y <- matrix(0,n_G,n)
  
  p <- floor(rep(n/4,n_G))
  p1 <- floor(rep(n/6,n_G))
  p2 <- (n/2)-p1
  
  moy1 <- rdunif(n_G,10,20) # 1 10
  
  moy2 <- 3*moy1
  moy3 <- (moy1+moy2)/2
  
  X <- c(rep(0,n/2),rep(1,n/2))
  
  bimod <- rep(0,n/2)
  bimod1 <- rep(0,n/2)
  bimod2 <- rep(0,n/2)
  unimod <- rep(0,n/2)
  
  for (i in 1:n_G){
    set.seed(seeds[i])
    if (i<=250){ # DE
      
      Y[i,] <- c(rbinom(n = n/2, size = 1, p=prop_non_0)*rnbinom(n/2,prob=0.5,size=moy1[i]),
                 rbinom(n = n/2, size = 1, p=prop_non_0)*rnbinom(n/2,prob=0.5,size=moy3[i]))
    }
    
    if (i<=500 & i>250){ # DM
      
      bimod <- rep(0,n/2)
      unimod <- rep(0,n/2)
      
      bimod[1:p[i]] <- c(rbinom(n = floor(length(1:p[i])), size = 1, p=prop_non_0)*rnbinom(floor(length(1:p[i])),prob=0.5,size=moy1[i]))
      bimod[(p[i]+1):(n/2)] <- c(rbinom(n = floor(length((p[i]+1):(n/2))), size = 1, p=prop_non_0)*rnbinom(floor(length((p[i]+1):(n/2))),prob=0.5,size=moy2[i]))
      unimod <- c(rbinom(n = n/2, size = 1, p=prop_non_0)*rnbinom(n/2,prob=0.5,size=moy2[i]))
      Y[i,] <- c(unimod,bimod)
      
    }
    
    if (i<=750 & i>500){ # DP
      
      bimod1 <- rep(0,n/2)
      bimod2 <- rep(0,n/2)
      
      bimod1[1:p1[i]] <- c(rbinom(n = length(1:p1[i]), size = 1, p=prop_non_0)*rnbinom(length(1:p1[i]),prob=0.5,size=moy1[i]))
      bimod1[(p1[i]+1):(n/2)] <- c(rbinom(n = length((p1[i]+1):(n/2)), size = 1, p=prop_non_0)*rnbinom(length((p1[i]+1):(n/2)),prob=0.5,size=moy2[i]))
      
      bimod2[1:p2[i]] <- c(rbinom(n = length(1:p2[i]), size = 1, p=prop_non_0)*rnbinom(length(1:p2[i]),prob=0.5,size=moy1[i]))
      bimod2[(p2[i]+1):(n/2)] <- c(rbinom(n = length((p2[i]+1):(n/2)), size = 1, p=prop_non_0)*rnbinom(length((p2[i]+1):(n/2)),prob=0.5,size=moy2[i]))
      
      Y[i,] <- c(bimod1,bimod2)
    }
    
    if (i<=2000 & i>750){ # DB
      
      bimod <- rep(0,n/2)
      unimod <- rep(0,n/2)
      
      bimod[1:p[i]] <- c(rbinom(n = length(1:p[i]), size = 1, p=prop_non_0)*rnbinom(length(1:p[i]),prob=0.5,size=moy1[i]))
      bimod[(p[i]+1):(n/2)] <- c(rbinom(n = length((p[i]+1):(n/2)), size = 1, p=prop_non_0)*rnbinom(length((p[i]+1):(n/2)),prob=0.5,size=moy2[i]))
      unimod <- c(rbinom(n = n/2, size = 1, p=prop_non_0)*rnbinom(n/2,prob=0.5,size=moy3[i]))
      Y[i,] <- c(unimod,bimod)
      
    }
    
    if (i<=5500 & i>1000){
      
      Y[i,] <- c(rbinom(n = n/2, size = 1, p=prop_non_0)*rnbinom(n/2,prob=0.5,size=moy1[i]),
                 rbinom(n = n/2, size = 1, p=prop_non_0)*rnbinom(n/2,prob=0.5,size=moy1[i]))
      
    }
    
    if (i>5500){
      
      bimod1 <- rep(0,n/2)
      bimod2 <- rep(0,n/2)
      
      bimod1[1:p[i]] <- c(rbinom(n = length(1:p[i]), size = 1, p=prop_non_0)*rnbinom(length(1:p[i]),prob=0.5,size=moy1[i]))
      bimod1[(p[i]+1):(n/2)] <-   c(rbinom(n = length((p[i]+1):(n/2)), size = 1, p=prop_non_0)*rnbinom(length((p[i]+1):(n/2)),prob=0.5,size=moy2[i]))
      
      bimod2[1:p[i]] <- c(rbinom(n = length(1:p[i]), size = 1, p=prop_non_0)*rnbinom(length(1:p[i]),prob=0.5,size=moy1[i]))
      bimod2[(p[i]+1):(n/2)] <-  c(rbinom(n = length((p[i]+1):(n/2)), size = 1, p=prop_non_0)*rnbinom(length((p[i]+1):(n/2)),prob=0.5,size=moy2[i]))
      
      Y[i,] <- c(bimod1,bimod2)
      
    }
  }
  
  return(list(Y=Y,X=X))
}


######### Results #########

size <- c(20,40,60,80,100,160,200)
methods <- c("CCDF_asymp", "CCDF_perm","MAST","scDD","SigEMD","DESingle","SCDE","D3E")

### taskid 1:56
slar_taskid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

s <- (slar_taskid - 1) %/% length(1:8) + 1
m <- (slar_taskid - 1) %% (length(1:8)) + 1

temp <- sample_mat_NBP(n_G = 10000, n = size[s])

X <- temp$X
genes_matrix <- temp$Y
Y <- genes_matrix
ptm <- proc.time()
if( m < 7){
  sce <- SingleCellExperiment(list(counts=Y))
  ngenes <- nrow(genes_matrix)
  rownames(sce) <- seq_len(ngenes)
  assayNames(sce) <- "counts"
  colData(sce)$condition <- X
  ncells <- size[s]
  
  # CCDF asymptotic test
  if(m %in% 1:2){
    Y <- sce@assays@data$counts
    rownames(Y) <- seq_len(ngenes)
    colnames(Y) <- seq_len(ncells)
    names(X) <- colnames(Y)
    if(m==1){ # CCDF ASYMP
      res_pvalue <- as.numeric(ccdf::ccdf_testing(data.frame(Y = Y), data.frame(X = as.factor(X)),
                                                  test = "asymptotic", n_cpus = 16,
                                                  space_y = TRUE, number_y = (ncol(Y) / 2))$pvals$raw_pval)
    }else{
      # CCDF PERM
      source("ccdfModif/ccdf_testing_Mod.R")
      res_pvalue <- ccdf_testing_Modif(data.frame(Y = Y), data.frame(X = as.factor(X)),
                                       n_perm=500,space_y = TRUE, 
                                       number_y = (ncol(Y) / 2))$pvals
    }
  }
  # MAST
  if(m==3){
    fData <- data.frame(primerid=seq_len(ngenes))
    cData <- data.frame(wellKey=seq_len(ncells))
    sca <- FromMatrix(as.matrix(sce@assays@data$counts), cData, fData,check_sanity = FALSE)
    assayNames(sca) <- "TPM"
    cond <- factor(X)
    cond <- relevel(cond,"1")
    colData(sca)$condition <- cond
    zlmCond <- zlm(~condition, sca)
    summaryCond <- summary(zlmCond,doLRT="condition2")
    summaryDt <- summaryCond$datatable
    fcHurdle <- summaryDt[contrast=='condition2' & component=='H',.(`primerid`,`Pr(>Chisq)`)] #hurdle P values
    res_pvalue <- fcHurdle$`Pr(>Chisq)`[sort(as.numeric(fcHurdle$primerid),index.return =TRUE)$ix]
  }
  
  # scDD
  if(m==4){
    sca <- sce
    rownames(sca) <- seq_len(ngenes)
    assayNames(sca) <- "normcounts"
    scDatExSim <- scDD(sca, prior_param=list(alpha=0.01, mu0=0, 
                                             s0=0.01, a0=0.01, 
                                             b0=0.01),
                       testZeroes=FALSE,categorize = FALSE)
    RES <- results(scDatExSim)
    res_pvalue <- RES$nonzero.pvalue
  }
  
  # SigEMD
  if(m==5){
    # File script
    source("SigEMD-master/FunImpute.R")
    source("SigEMD-master/SigEMDHur.R")
    source("SigEMD-master/SigEMDnonHur.R")
    source("SigEMD-master/plot_sig.R")
    
    # Execute ...
    Y <- sce@assays@data$counts
    rownames(Y) <- seq_len(ngenes)
    colnames(Y) <- seq_len(ncells)
    names(X) <- colnames(Y)
    
    databinary <- databin(Y)
    Hur_gene<- idfyImpgene(Y,databinary,X)
    results<- calculate_single(data =  Y, condition =  X, Hur_gene = Hur_gene, nperm=500, binSize=0.2)
    emd <- results$emdall
    res_pvalue <- emd[,2]
  }
  
  # DESingle
  if(m==6){
    resDESingle <- DEsingle(counts = sce, group = as.factor(X))
    res_pvalue <- resDESingle[,"pvalue"]
  }
}

# SCDE || D3E
if(m>6){
  Y <- as.data.frame(Y)
  groupe <- as.factor(X)
  names(groupe) <- colnames(Y)
  
  # SCDE
  if(m==7){
    n.cores <- 4
    counts<-apply(Y,2,function(x) {storage.mode(x) <- 'integer'; x})
    scdeEr <- scde.error.models(counts=counts,
                                groups=groupe,
                                n.cores=n.cores,
                                save.model.plots=F)
    
    scdePrior <- scde.expression.prior(models=scdeEr,counts=counts)
    ediff <- scde.expression.difference(scdeEr,counts,scdePrior,
                                        groups=groupe,n.cores=n.cores)
    res_pvalue <- 2*pnorm(abs(ediff$Z),lower.tail=F)
  }else{
    # Python path in curta
    pythonPath <- "/gpfs/softs/contrib/apps/python/3.9.7/bin/python3"
    
    # File script python
    pythonD3ECmd <- "D3E/D3ECmd.py"
    # Execute ...
    namesRow <- rownames(data.frame(Y))
    Y <- rbind.data.frame(GeneID=paste0("G", X), Y)
    Y <- cbind.data.frame(firstCol=c("GeneID",paste0("ENSMUSG",namesRow)),Y)
    inputData <- Y[-1,]
    colnames(inputData) <- Y[1,]
    
    
    write.table(inputData, paste0("inputFile_",slar_taskid,".txt"), dec = ".",sep = "\t",
                row.names = F, col.names = T,quote = F)
    
    args <- c(pythonD3ECmd,paste0("inputFile_",slar_taskid,".txt"), paste0("reT_",slar_taskid,".txt"), "G1", "G2",
              "-m", 1, "-t", 0, "-z", 0, "-n", 1, "-v")
    pvals <- system2(pythonPath, args=args, stdout=TRUE)
    pvals <- pvals[length(pvals)]
    pvals <- as.vector(unlist(strsplit(pvals, "[,]")))
    res_pvalue <- as.numeric(sub("[", "", sub("]", "", pvals, fixed = TRUE), fixed = TRUE))
    
  }
}
timeALL <- proc.time() - ptm

write.table(data.matrix(timeALL)[1:3],
            file=paste0("time/analysis_2cond_NB_",methods[m],"_",size[s],".txt"),
            row.names = FALSE, col.names = TRUE, sep = "\t")

write.table(res_pvalue,
            file=paste0("output/analysis_2cond_NB_",methods[m],"_",size[s],".txt"),
            row.names = FALSE, col.names = TRUE, sep = "\t")


print("DONE!!!!!!!!!!!")

