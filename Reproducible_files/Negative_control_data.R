# download files at https://doi.org/10.5281/zenodo.5701445

library(scran)
library(scater)
library(ccdf)

load("seed.RData")

incoming <- read.table("GSE54695_data_transcript_counts.txt", header=TRUE, row.names=1, sep="\t", 
                       colClasses=c("character", rep("numeric", 80), vector("list", 80), 
                                    rep("numeric", 80), vector("list", 80))) # keep only single cells.
spike.in <- grepl("^ERCC", rownames(incoming))
grouping <- sub("SC_([^_]+)_[0-9]+", "\\1", colnames(incoming))
design <- model.matrix(~grouping)
name <- "mESC (Grun)"


res <- rep(NA,10)
load("negative.RData")
negative <- negative[,1:80]
v_g <- matrixStats::rowVars(as.matrix(negative))
ind <- which(v_g==0)
Y <- negative[-ind,]


for (i in 1:10){
  print(i)
  slar_taskid <- i
  set.seed(seed[slar_taskid])
  
  index <- sample(1:80,size=80)

  X <- c(rep(1,40),rep(2,40))[index]
  Y <- calculateTPM(as.matrix(Y))
  Y <- log(Y+1,base = 2)
  
  # CCDF asymptotic test
  
  res_ccdf_asymp <- cit(data.frame(Y=Y), data.frame(X=as.factor(X)), test="asymptotic", n_cpus = 16)$pvals$raw_pval
  
  res[i] <- length(which(p.adjust(res_ccdf_asymp,"BH")<0.05))
  
}
