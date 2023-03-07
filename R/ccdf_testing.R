#' Main function to perform complex hypothesis testing using (un)conditional independence test
#'
#'@param exprmat a data frame of size \code{G x n} containing the
#'preprocessed expressions from \code{n} samples (or cells) for \code{G}
#'genes. Default is \code{NULL}.
#'
#'@param variable2test a data frame of numeric or factor vector(s) 
#'of size \code{n} containing the variable(s) to be tested (the condition(s))
#' 
#'@param covariate a data frame of numeric or factor vector(s) 
#'of size \code{n} containing the covariate(s)
#'
#'@param test a character string indicating which method to use to
#'compute the test, either \code{'asymptotic'} or \code{'permutations'}.
#'Default is \code{'asymptotic'}.
#'
#'@param method a character string indicating which method to use to
#'compute the CCDF, either \code{'linear regression'}, \code{'logistic regression'}
#' and  \code{'permutations'} or \code{'RF'} for Random Forests.
#'Default is \code{'linear regression'} since it is the method used in the test.
#'
#'@param fast a logical flag indicating whether the fast implementation of
#'logistic regression should be used. Only if \code{'dist_permutations'} is specified.
#'Default is \code{TRUE}.
#'
#'@param n_perm the number of permutations. Default is \code{100}.
#'
#'@param adaptive a logical flag indicating whether adaptive permutations
#'should be performed. Default is \code{FALSE}.
#'
#'@param n_perm_adaptive a vector of the increasing numbers of 
#'adaptive permutations when \code{adaptive} is \code{TRUE}. 
#'\code{length(n_perm_adaptive)} should be equal to \code{length(thresholds)+1}. 
#'Default is \code{c(100,150,250,500)}.
#'
#'@param thresholds a vector of the decreasing thresholds to compute
#'adaptive permutations when \code{adaptive} is \code{TRUE}. 
#'\code{length(thresholds)} should be equal to \code{length(n_perm_adaptive)-1}.
#'Default is \code{c(0.1,0.05,0.01)}.
#'
#'@param distance a character string indicating which distance to use to
#'compute the test, either \code{'L2'}, \code{'L1'} or 
#'\code{'L_sup'}, when \code{method} is \code{'dist_permutations'}, 
#'Default is \code{'L2'}.
#'
#'@param parallel a logical flag indicating whether parallel computation
#'should be enabled. Default is \code{TRUE}.
#'
#'@param n_cpus an integer indicating the number of cores to be used for the computations.
#'Default is \code{parallel::detectCores() - 1}. If \code{n_cpus = 1}, then sequential 
#'computations are used without any parallelization.
#'
#'@param space_y a logical flag indicating whether the y thresholds are spaced. 
#'When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and 
#'the maximum of the observations is used. If \code{FALSE}, each unique 
#'observed expression value is used as a distinct threshold. Default is \code{TRUE}.
#'
#'@param number_y an integer value indicating the number of y thresholds (and therefore
#'the number of regressions) to perform the test. Default is 10.
#'
#'
#'@import pbapply parallel
#'@importFrom doParallel registerDoParallel
#'
#'
#'@return A list with the following elements:\itemize{
#'   \item \code{which_test}: a character string carrying forward the value of
#'   the '\code{which_test}' argument indicating which test was performed (either
#'   'asymptotic','permutations','dist_permutations').
#'   \item \code{n_perm}: an integer carrying forward the value of the
#'   '\code{n_perm}' argument or '\code{n_perm_adaptive}' indicating the number of permutations performed
#'   (\code{NA} if asymptotic test was performed).
#'   \item \code{pval}: computed p-values. A data frame with one raw for
#'   each gene, and with 2 columns: the first one '\code{raw_pval}' contains
#'   the raw p-values, the second one '\code{adj_pval}' contains the FDR adjusted p-values
#'   using Benjamini-Hochberg correction.
#' }
#' 
#'@references Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2019).
#'Distribution-free complex hypothesis testing for single-cell RNA-seq differential expression analysis, *bioRxiv* 445165.
#'[DOI: 10.1101/2021.05.21.445165](https://doi.org/10.1101/2021.05.21.445165).
#'
#'@export
#'
#'@examples
#' 
#'ncells <- 100
#'pgenes <- 500
#'X1 <- as.factor(rbinom(n=ncells, size = 1, prob = 0.5))
#'Y <- t(replicate(pgenes, ((X1==1)*rnorm(n = ncells,0,1)) + ((X1==0)*rnorm(n = ncells, 0.5, 1))))
#'X2 <- rnorm(n=ncells)
#'X3 <- rnorm(n=ncells)
#'Z1 <- rnorm(ncells)
#'Z2 <- as.factor(rbinom(n=ncells, size=1, prob = 0.5))
#'
#'res_asymp <- ccdf_testing(exprmat=data.frame(Y=Y), 
#'variable2test=data.frame(X1=X3, X2 = X2), 
#'covariate = data.frame(Z1 = Z1, Z2 = Z2),
#'test="asymptotic", n_cpus=1, parallel=FALSE)
#'hist(res_asymp$pvals$raw_pval) # asymptotic test




ccdf_testing <- function(exprmat = NULL,
                         variable2test = NULL,
                         covariate = NULL,
                         distance = c("L2","L1","L_sup"),
                         test = c("asymptotic","permutations"),
                         method = c("linear regression","logistic regression","RF"),
                         fast = TRUE,
                         n_perm = 100,
                         n_perm_adaptive = c(n_perm, n_perm, n_perm*3, n_perm*5),
                         thresholds = c(0.1,0.05,0.01),
                         parallel = TRUE,
                         n_cpus = NULL,
                         adaptive = FALSE,
                         space_y = TRUE,
                         number_y = 10){
  
  # check
  stopifnot(is.data.frame(exprmat))
  stopifnot(is.data.frame(variable2test))
  stopifnot(is.data.frame(covariate) | is.null(covariate))
  stopifnot(is.logical(parallel))
  stopifnot(is.logical(fast))
  stopifnot(is.logical(adaptive))
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
  
  
  if (length(method) > 1) {
    method <- method[1]
  }
  stopifnot(method %in% c("linear regression","logistic regression"))
  
  if (length(test) > 1) {
    test <- test[1]
  }
  stopifnot(test %in% c("asymptotic","permutations"))
  
  if (length(distance) > 1) {
    distance <- distance[1]
  }
  stopifnot(distance %in% c("L2","L1","L_sup"))
  
  if (test == "permutations"){
    
    if (adaptive){
      if ((length(n_perm_adaptive)!=(length(thresholds)+1))){
        warning("length of thresholds + 1 must be equal to length of n_perm_adaptive. \n",
                "Consider using the default parameters.")
      }
    }
    
    else{
      N_possible_perms <- factorial(ncol(exprmat))
      if (n_perm > N_possible_perms){
        warning("The number of permutations requested 'n_perm' is ",
                n_perm, "which is larger than the total number of ",
                "existing permutations ", N_possible_perms,
                ". Try a lower number for 'n_perm' (currently ",
                "running with 'nperm=", N_possible_perms, "').")
        n_perm <- N_possible_perms
      }
    }
  }
  
  if (space_y){
    if (is.null(number_y)){
      warning("Missing argument", number_y, ". No spacing is used.")
      space_y <- FALSE
    }
  }
  
  # parallel
  
  if(parallel){
    if(is.null(n_cpus)){
      n_cpus <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(n_cpus)
    doParallel::registerDoParallel(cl)
  }
  else{
    n_cpus <- 1
    cl<-1
  }
  
  # Test ----
  ## permutations ----
  if (test=="permutations"){
    
    if (adaptive==TRUE){ 
      #### adaptive ----
      message(paste("Computing", n_perm_adaptive[1], "permutations..."))
      
      res <- pbapply::pbsapply(1:nrow(exprmat), 
                               function(i){
                                 test_perm(
                                   Y = exprmat[i,],
                                   X = variable2test,
                                   Z = covariate,
                                   n_perm = n_perm_adaptive[1],
                                   space_y = space_y, 
                                   number_y = number_y)$score
                               },
                               cl=cl)
      perm <- rep(n_perm_adaptive[1],nrow(exprmat))
      
      k <- 2
      index <- which((res+1)/(perm+1) <= thresholds[k-1])
      
      while (length(index)!=0 & k<=length(n_perm_adaptive)){
        
        index <- which(((res+1)/(perm+1)) < thresholds[k-1])
        
        message(paste("Computing", sum(n_perm_adaptive[k]), "additional permutations..."))
        
        res_perm <- pbapply::pbsapply(1:length(index), 
                                      function(i){
                                        test_perm(Y = exprmat[index[i], ],
                                                  X = variable2test,
                                                  Z = covariate,
                                                  n_perm = n_perm_adaptive[k],
                                                  space_y = space_y, 
                                                  number_y = number_y)$score
                                      }, 
                                      cl = cl)
        res[index] <- res[index] + res_perm
        perm[index] <- perm[index] + rep(n_perm_adaptive[k], length(index))
        k <- k+1
        
      }
      
      pvals <- (res+1)/(perm+1)
      df <- data.frame(raw_pval = pvals,
                       adj_pval = p.adjust(pvals, method = "BH"))
      
      n_perm <- cumsum(n_perm_adaptive)
      
    }else{ 
      #### non-adaptive ----
      message(paste("Computing", n_perm, "permutations..."))
      res <- do.call("rbind",pbapply::pblapply(1:nrow(exprmat), 
                                               function(i){
                                                 test_perm(Y = exprmat[i,],
                                                           X = variable2test,
                                                           Z = covariate,
                                                           n_perm = n_perm,
                                                           space_y = space_y, 
                                                           number_y = number_y)
                                               },
                                               cl=cl))
      
      #res <- as.vector(unlist(res))
      
      df <- data.frame(raw_pval = res$raw_pval,
                       adj_pval = p.adjust(res$raw_pval, method = "BH"))
      
    }
    
    
  }
  
  ## asymptotic ----
  else if (test=="asymptotic"){
    n_perm <- NA
    
    Y <- exprmat
    X <- variable2test
    Z <- covariate
    res <- do.call("rbind",pbapply::pblapply(1:nrow(Y), 
                                             function(i){
                                               test_asymp(Y[i,], X, Z, 
                                                          space_y = space_y, 
                                                          number_y = number_y)
                                             }, 
                                             cl=cl)
    )
    df <- data.frame(raw_pval = res$raw_pval, adj_pval = p.adjust(res$raw_pval, method = "BH"), test_statistic = res$Stat)
  }
  
  if(parallel){
    parallel::stopCluster(cl)
  }
  
  #rownames(df) <- genes_names
  
  
  return(list(which_test = test,
              n_perm = n_perm, 
              pvals = df))
  
  
}
