# download files at https://doi.org/10.5281/zenodo.5701445

library(ccdf)

load("mat_mouse.RData")
top1000 <- read.table("Top1000_Jaakkola.txt")
top1000 <- as.character(top1000$V1)
genes_names <- rownames(mat_mouse)

X <- condition
Y <- mat_mouse

# CCDF asymptotic test

res_ccdf_asymp <- ccdf_testing(Y, X, test="asymptotic", n_cpus = 16)$pvals$raw_pval

nb_DE_genes <- length(which(p.adjust(res_ccdf_asymp,"BH")<0.05))
nb_inter_gold <-  length(intersect(genes_names[which(p.adjust(res_ccdf_asymp,"BH")<0.05)],top1000))
