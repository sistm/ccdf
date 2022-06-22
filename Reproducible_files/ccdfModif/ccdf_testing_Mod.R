ccdf_testing_Modif <- function(exprmat = NULL,
                               variable2test = NULL,
                               covariate = NULL,
                               n_perm = 100,
                               n_cpus = NULL,
                               space_y = FALSE,
                               number_y = ncol(exprmat)){
  source("ccdfModif/test_perm_Modif.R")
  
  # check
  stopifnot(is.data.frame(exprmat))
  stopifnot(is.data.frame(variable2test))
  stopifnot(is.data.frame(covariate) | is.null(covariate))
  stopifnot(is.numeric(n_perm))
  
  genes_names <- rownames(exprmat)
  
  if (sum(is.na(exprmat)) > 1) {
    warning("'y' contains", sum(is.na(exprmat)), "NA values. ",
            "\nCurrently they are ignored in the computations but ",
            "you should think carefully about where do those NA/NaN ",
            "come from...")
    exprmat <- exprmat[complete.cases(exprmat),]
  }
  
  
  # checking for 0 variance genes
  v_g <- matrixStats::rowVars(as.matrix(exprmat))
  if(sum(v_g==0) > 0){
    warning("Removing ", sum(v_g==0), " genes with 0 variance from ",
            "the testing procedure.\n",
            "  Those genes should probably have been removed ",
            "beforehand...")
    exprmat <- exprmat[v_g>0, ]
  }
  
  
  
  N_possible_perms <- factorial(ncol(exprmat))
  if (n_perm > N_possible_perms){
    warning("The number of permutations requested 'n_perm' is ",
            n_perm, "which is larger than the total number of ",
            "existing permutations ", N_possible_perms,
            ". Try a lower number for 'n_perm' (currently ",
            "running with 'nperm=", N_possible_perms, "').")
    n_perm <- N_possible_perms
  }
  
  
  if (space_y){
    if (is.null(number_y)){
      warning("Missing argument", number_y, ". No spacing is used.")
      space_y <- FALSE
    }
  }
  
  
  # test
  
  print(paste("Computing", n_perm, "permutations..."))
  res <- do.call("rbind",pbapply::pblapply(1:nrow(exprmat), FUN=function(i){
    test_perm_Modif(Y = exprmat[i,],
                    X = variable2test,
                    Z = covariate,
                    n_perm = n_perm,
                    space_y = space_y,
                    number_y = number_y)},cl=n_cpus))
  
  return(list(which_test = "permut",
              n_perm = n_perm, 
              pvals = res))
}
