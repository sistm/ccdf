# download files at https://doi.org/10.5281/zenodo.5701445

library(fst)
library(citcdf)
load("genes.RData")

path <- paste0("data", ".fst")
ft <- fst::read_fst(path) # single-cell matrix, log_2(CPM) 
# "only transcripts expressed in at least 0.1\% of the cells were included in the differential analyses, yielding 10,525 genes"

# metadata
metadata <- read.delim("GSE153931_cd8_t24_processed_data_annotations.txt", header = TRUE, sep = "\t", dec = ".")

# remove patient void
void <- which(as.factor(metadata$orig.donor)=="void")
ft <- ft[,-as.numeric(void)]
metadata <- metadata[-void,]

# keep covid
illness <- rep(1,length(metadata$orig.donor))
for (k in 1:8){
  illness[which(metadata$orig.donor==unique(metadata$orig.donor)[24+k])] <- rep(0,length(which(
    metadata$orig.donor==unique(metadata$orig.donor)[24+k])))
}
covid <- metadata$orig.unit[which(illness==1)] # X
covid_cluster <- metadata$RNA_snn_res.0.2[which(illness==1)] # Z
ft <- ft[,which(illness==1)] # matrix

# arrange X
covid[which(as.factor(covid)=="void")] <- "Not_admitted"

covid_bin <- rep(1,length(covid))
covid_bin[which(as.factor(covid)=="Not_admitted")] <- rep(0,length(which(as.factor(covid)=="Not_admitted")))

# arrange Z

keep_cluster <- which(as.factor(covid_cluster)!=7)

### ccdf method (example)

# one gene
i=1
res <- cit_asymp(as.numeric(ft[i,keep_cluster]),data.frame(X=factor(covid,levels=c("ITU","Ward","Not_admitted"))[keep_cluster]), 
                      #data.frame(Z=factor(covid_cluster,levels=c("0","1","2","3","4","5","6"))[keep_cluster]),
                      space_y = TRUE, number_y = 10) 

# all genes

res <- cit(as.data.frame(ft[,keep_cluster]),data.frame(X=factor(covid,levels=c("ITU","Ward","Not_admitted"))[keep_cluster]), 
                  #data.frame(Z=factor(covid_cluster,levels=c("0","1","2","3","4","5","6"))[keep_cluster]),
                  space_y = TRUE, number_y = 10, parallel = TRUE, n_cpus = 16)

