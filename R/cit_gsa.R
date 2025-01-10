
#' Conditional independance test for gene set analysis
#'
#' @param M a \code{data.frame} or a \code{matrix} of size \code{n x r} 
#'containing the different Y variables to test for conditional independence 
#'with \code{X} adjusted on \code{Z}.
#'
#' @param X a data frame of size \code{n x p} of numeric or factor vector(s) 
#'containing the variable(s) to be tested for conditional independence 
#'against \code{X} adjusted on \code{Z}. Multiple variables (\code{p>1}) 
#'are only supported by the asymptotic test.
#'
#' @param Z a data frame of size \code{n x q} of numeric or factor vector(s) 
#'containing the covariate(s) to condition the independence 
#'test upon. Multiple covariates (\code{q>1}) are only supported by the 
#'asymptotic test. 
#'
#' @param geneset a vector, a list, a gmt file format or a BiocSet object. 
#'If the parameter is \itemize{
#'  \item a vector : corresponds to the gene name of the gene set, must be the same as those of the columns of the matrix \code{M}
#'  \item a list : each elements of the list are a gene set with the names of the genes, must be the same as those of the columns of the matrix \code{M}
#'  \item a gmt file format : the genes names of each genes set in the file, must be the same as those of the columns of the matrix \code{M}
#'  \item a BiocSet object : the genes names of each genes set in the object, must be the same as those of the columns of the matrix \code{M}
#'}
#'
#'@param test a character string indicating whether the \code{'asymptotic'} or 
#'the \code{'permutation'} test is computed.
#'Default is \code{'asymptotic'}.
#'
#'@param n_perm the number of permutations. Default is \code{100}. Only used if
#'\code{test == 'permutation'}.
#'
#'@param adaptive a logical flag indicating whether adaptive permutations
#'should be performed. Default is \code{TRUE}. Only used if
#'\code{test == 'permutation'}.
#'
#'@param n_perm_adaptive a vector of the increasing numbers of 
#'adaptive permutations when \code{adaptive} is \code{TRUE}. 
#'\code{length(n_perm_adaptive)} should be equal to \code{length(thresholds)+1}. 
#'Default is \code{c(100, 150, 250, 500)}.
#'
#'@param thresholds a vector of the decreasing thresholds to compute
#'adaptive permutations when \code{adaptive} is \code{TRUE}. 
#'\code{length(thresholds)} should be equal to \code{length(n_perm_adaptive)-1}.
#'Default is \code{c(0.1, 0.05, 0.01)}.
#'
#'
#'@param parallel a logical flag indicating whether parallel computation
#'should be enabled. Default is \code{TRUE} if \code{interactive()} is 
#'\code{TRUE}, else is \code{FALSE}.
#'
#'@param n_cpus an integer indicating the number of cores to be used for the 
#'computations. Default is \code{parallel::detectCores() - 1}. If 
#'\code{n_cpus = 1}, then sequential computations are used without any 
#'parallelization.
#'
#'@param space_y a logical flag indicating whether the y thresholds are spaced out. 
#'When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and 
#'the maximum of the observations is used. If \code{FALSE}, each unique 
#'observed expression value is used as a distinct threshold. Default is \code{TRUE}.
#'
#'@param number_y an integer value indicating the number of y thresholds (and therefore
#'the number of regressions) to perform the test. Default is 10.
#'   
#'@return A list with the following elements:\itemize{
#'   \item \code{which_test}: a character string carrying forward the value of
#'   the '\code{which_test}' argument indicating which test was performed (either
#'   'asymptotic' or 'permutation').
#'   \item \code{n_perm}: an integer carrying forward the value of the
#'   '\code{n_perm}' argument or '\code{n_perm_adaptive}' indicating the number of permutations performed
#'   (\code{NA} if asymptotic test was performed).
#'   \item \code{pval}: computed p-values. A data frame with one raw for
#'   each gene set, and with 2 columns: the first one '\code{raw_pval}' contains
#'   the raw p-values, the second one '\code{adj_pval}' contains the FDR adjusted p-values
#'   using Benjamini-Hochberg correction.
#' }
#' 
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 500
#' r <- 200
#' Z1 <- rnorm(n)
#' Z2 <- rnorm(n)#rbinom(n, size=1, prob=0.5) + rnorm(n, sd=0.05)
#' X1 <- Z2 + rnorm(n, sd=0.2)
#' X2 <- rnorm(n)
#' cor(X1, Z2)
#' Y <- replicate(r, Z2) + rnorm(n*r, 0, 0.5)
#' range(cor(Y, Z2))
#' range(cor(Y, X2))
#' res_asymp_unadj <- cit_gsa(M = data.frame(Y=Y), 
#'                            X = data.frame(X2=X1), 
#'                            geneset = paste0("Y.", 1:50),
#'                            test="asymptotic", parallel=FALSE)
#'res_asymp_unadj$pvals   
cit_gsa <- function(M,
                    X,
                    Z = NULL,
                    geneset,
                    test = c("asymptotic","permutation"),
                    n_perm = 100,
                    n_perm_adaptive = c(n_perm, n_perm, n_perm*3, n_perm*5),
                    thresholds = c(0.1,0.05,0.01),
                    parallel = interactive(),
                    n_cpus = detectCores() - 1,
                    adaptive = FALSE,
                    space_y = TRUE,
                    number_y = 10){
  
  # checks
  
  stopifnot(is.data.frame(M) | is.matrix(M))
  stopifnot(is.data.frame(X))
  stopifnot(is.data.frame(Z) | is.null(Z))
  stopifnot(is.logical(parallel))
  stopifnot(is.logical(adaptive))
  stopifnot(is.numeric(n_perm))
  stopifnot(inherits(geneset,"GSA.genesets") | inherits(geneset,"BiocSet") | is.character(geneset) | is.list(geneset))
  
  
  M_colnames <- colnames(M)
  
  
  if (sum(is.na(M)) > 1) {
    warning("'M' contains", sum(is.na(M)), "NA values. ",
            "\nCurrently they are ignored in the computations but ",
            "you should think carefully about where do those NA/NaN ",
            "come from...")
    M <- M[,complete.cases(t(M))]
  }
  
  r <- ncol(M)
  n <- nrow(M)
  stopifnot(nrow(X) == n)
  stopifnot(nrow(Z) == n | is.null(Z))
  
  if (length(test) > 1) {
    test <- test[1]
  }
  stopifnot(test %in% c("asymptotic","permutation"))
  
  if (test == "permutation"){
    
    if (adaptive){
      if ((length(n_perm_adaptive)!=(length(thresholds)+1))){
        warning("length of thresholds + 1 must be equal to length of n_perm_adaptive. \n",
                "Consider using the default parameters.")
      }
    }else{
      N_possible_perms <- factorial(n)
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
    
    if(.Platform$OS.type == "unix"){
      par_clust <- n_cpus
    }else{
      par_clust <- parallel::makeCluster(n_cpus)
    }
    
  }else{
    par_clust <- 1L
  }
  
  
  
  
  # Test ----
  ## permutations ----
  if (test=="permutation"){
    
    stopifnot(ncol(X) < 2)
    stopifnot(ncol(Z) < 2 | is.null(Z))
    
    if (adaptive==TRUE){ 
      #### adaptive ----
      message(paste("Computing", n_perm_adaptive[1], "permutations..."))
      
      res <- pbapply::pbsapply(1:r, 
                               function(j){
                                 cit_perm(
                                   Y = M[, j],
                                   X = X,
                                   Z = Z,
                                   n_perm = n_perm_adaptive[1],
                                   space_y = space_y, 
                                   number_y = number_y)$score
                               },
                               cl=par_clust)
      perm <- rep(n_perm_adaptive[1], r)
      
      k <- 2
      index <- which((res+1)/(perm+1) <= thresholds[k-1])
      
      while (length(index)!=0 & k<=length(n_perm_adaptive)){
        
        index <- which(((res+1)/(perm+1)) < thresholds[k-1])
        
        message(paste("Computing", sum(n_perm_adaptive[k]), "additional permutations..."))
        
        if(parallel & .Platform$OS.type == "unix"){
          res_perm <- pbapply::pbsapply(1:length(index), 
                                        function(i){
                                          cit_perm(Y = M[, index[i]],
                                                   X = X,
                                                   Z = Z,
                                                   n_perm = n_perm_adaptive[k],
                                                   space_y = space_y, 
                                                   number_y = number_y)$score
                                        }, 
                                        cl = par_clust, mc.preschedule=TRUE)
        }else{
          res_perm <- pbapply::pbsapply(1:length(index), 
                                        function(i){
                                          cit_perm(Y = M[, index[i]],
                                                   X = X,
                                                   Z = Z,
                                                   n_perm = n_perm_adaptive[k],
                                                   space_y = space_y, 
                                                   number_y = number_y)$score
                                        }, 
                                        cl = par_clust)
        }
        
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
      
      if(parallel & .Platform$OS.type == "unix"){
        res <- do.call("rbind", pbapply::pblapply(1:r, 
                                                  function(j){
                                                    cit_perm(Y = M[, j],
                                                             X = X,
                                                             Z = Z,
                                                             n_perm = n_perm,
                                                             space_y = space_y, 
                                                             number_y = number_y)
                                                  },
                                                  cl=par_clust, mc.preschedule=TRUE))
      }else{
        res <- do.call("rbind", pbapply::pblapply(1:r, 
                                                  function(j){
                                                    cit_perm(Y = M[, j],
                                                             X = X,
                                                             Z = Z,
                                                             n_perm = n_perm,
                                                             space_y = space_y, 
                                                             number_y = number_y)
                                                  },
                                                  cl=par_clust))
      }
      
      #res <- as.vector(unlist(res))
      
      df <- data.frame(raw_pval = res$raw_pval,
                       adj_pval = p.adjust(res$raw_pval, method = "BH"))
      
    }
    
    
  } else if (test=="asymptotic"){
    ## asymptotic ----
    n_perm <- NA
    
    
    
    # Data formatting in list format +  check column names
    if (inherits(geneset,"GSA.genesets")) { 
      geneset <- geneset$genesets
    } else if (inherits(geneset,"BiocSet")){
      geneset<- BiocSet::es_elementset(geneset)
      geneset <- lapply(X  = unique(geneset$set),
                        FUN = function(x){
                          geneset[geneset$set == x,]$element
                        }
      )
    } else if(is.vector(geneset) & !is.list(geneset)){
      geneset <- list(geneset)
    }
    
    
    
    if (is.null(Z)){ 
      colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
      modelmat <- as.matrix(model.matrix(~.,data=X))
    } else {
      colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
      colnames(Z) <- sapply(1:ncol(Z), function(i){paste0('Z',i)})
      modelmat <- as.matrix(model.matrix(~.,data=cbind(X,Z)))
    }
    
    indexes_X <- which(substring(colnames(modelmat), 1, 1) == "X")
    
    
    # Initialisation for each gene set 
    test_stat_list <- list()
    Sigma2_list <- list()
    decomp_list <- list()
    pval <- NA
    ccdf_list <- list()
    
    
    n_Y_all <- nrow(M)
    H <- n_Y_all*(solve(crossprod(modelmat)) %*% t(modelmat))[indexes_X, , drop=FALSE] 
    # length of Y, same for each genes because X and Y are the same 
    
    
    
    
    if(length(geneset)<3){
      pboptions(type="none")
    }
    
    n_cpus <- min(length(geneset), n_cpus)
    
    
    res <- pbapply::pblapply(1:length(geneset), function(k){ # 1 -- each list of gene set ----
      
      # Initialisation for each gene in the gene set
      #test_stat_gs <- NULL
      prop_gs <- list()
      indi_pi_gs <- list()
      #ccdf_gs <- list()
      
      
      measured_genes <- intersect(M_colnames, geneset[[k]])
      
      if(length(measured_genes)<1){ # check 1 : none genes of the current geneset are in M
        warning("0 genes from geneset ", k, " observed in expression data")
        pval <- NA
        test_stat_list <- NA
      }else{ 
        
        test_stat_gs <- numeric(length(measured_genes))
        
        if(length(measured_genes) < length(geneset[[k]])){ # check 2 : some genes of the current geneset are not in M
          warning(" Some genes from geneset ", k, " are not observed in expression data")
        }
        
        
        # code below if all the genes in the gene set k are in M, else pval + stat de test = NA
        
        for (i in 1:length(measured_genes)){ # 2 -- each genes in the gene set k ----
          
          Y <- M[, measured_genes[i]]
          
          
          # 1) Test statistic computation ----
          if (space_y){
            y <- seq(from = ifelse(length(which(Y==0))==0, min(Y), min(Y[-which(Y==0)])),
                     to = max(Y[-which.max(as.matrix(Y))]), length.out = number_y)
          } else{
            y <- sort(unique(Y))
          }
          p <- length(y)
          
          index_jumps <- sapply(y[-p], function(i){sum(Y <= i)}) 
          beta <- c(apply(X = H[, order(Y), drop=FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, ]) / n_Y_all # same number than thresholds
          test_stat <- sum(beta^2) * n_Y_all
          
          test_stat_gs[i] <- test_stat # test statistic for each genes in the gene set 
          
          
          # 2) Pi computation ----
          indi_pi <- matrix(NA, n_Y_all, (p-1)) 
          
          for (j in 1:(p-1)){ 
            indi_Y <- 1 * (Y <= y[j])
            indi_pi[,j] <- indi_Y
          }
          
          indi_pi_gs[[i]] <- indi_pi
          prop <- colMeans(indi_pi)
          prop_gs[[i]] <- prop # prop for each genes in the gene set 
          
          
          #ccdf_gs[[i]] <- ccdf(Y=Y, X=X, Z=Z, method="OLS", fast=TRUE, space_y=space_y, number_y=number_y)
          
          
        } 
        
      }
      
      #ccdf_list[[k]] <- ccdf_gs
      #names(ccdf_list[[k]]) <- measured_genes # not geneset, if some genes are not in the data
      # utile comme le refait plus tard ???  
      
      
      indi_pi_gs_tab <- do.call(cbind, indi_pi_gs)
      prop_gs_vec <- unlist(prop_gs)
      n_g_t <- length(prop_gs_vec)
      
      # 3) Sigma matrix creation ----
      Sigma2 <- matrix(NA, n_g_t*nrow(H), n_g_t*nrow(H)) 

      n_gs_vec <- nrow(indi_pi_gs_tab)
      temp <- indi_pi_gs_tab - matrix(prop_gs_vec, nrow=n_gs_vec, ncol=n_g_t, byrow=TRUE)
      
        # new prop/new pi computation = the one of the gene set, here it's a matrix
      new_prop <- apply(temp, MARGIN=2, function(colvec){colvec%*%temp})/n_gs_vec + 
                        sapply(prop_gs_vec, function(s){s*prop_gs_vec})
            
      Sigma2 <- 1/n * tcrossprod(H) %x%  (new_prop - prop_gs_vec %x%  t(prop_gs_vec))
      
      decomp <- eigen(Sigma2, symmetric=TRUE, only.values=TRUE)
      
      pval <- survey::pchisqsum(sum(test_stat_gs), lower.tail = FALSE, 
                                df = rep(1, ncol(Sigma2)), 
                                a =decomp$values , method = "saddlepoint")
      
      
      
      return(list("pval" = pval, "test_stat_gs" = test_stat_gs)) #,"ccdf" = ccdf_list
      
      
      
    },
    cl = n_cpus)
    pboptions(type="timer")
    
 
    
    pvals <- sapply(res, "[[", "pval")
    
    test_stat_list <- lapply(res, "[[", "test_stat_gs")
    
  
    
    #ccdf <- lapply(res, "[[", "ccdf")
    
    # for (i in 1:length(ccdf)){
    #   if(length(ccdf[[i]])>1){
    #     ccdf[[i]] <- Filter(Negate(is.null), ccdf[[i]])
    #   }
    #   ccdf[[i]] <- lapply(seq_along(ccdf[[i]][[1]]), function(j) ccdf[[i]][[1]][[j]])
    #   names(ccdf[[i]]) <- intersect(M_colnames,geneset[[i]])
    #   }
    
    # ccdf_new<- lapply(ccdf, function(x) {
    #   if (length(x) > 1) {x <- Filter(Negate(is.null), x)}
    #   lapply(seq_along(x[[1]]), function(j) x[[1]][[j]])
    #   
    # })
    
   # récup genes dans measured gene sinon ??
  }
  
  df <- data.frame(raw_pval = pvals,
                   adj_pval = p.adjust(pvals, method="BH"),
                   test_statistic = sapply(test_stat_list, sum))
  
  #rownames(df) <- M_colnames
  
  
  output <- list(which_test = test,
                 n_perm = n_perm, 
                 pvals = df) #, ccdf = ccdf
  
  class(output) <- "citcdf"
  output$type <- "gsa"
  return(output)
  
}
