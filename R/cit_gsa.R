
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
#' @param geneset 
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
#' FINIR QUAND AURA FAIT CODE POUR GMT ET BIOCSET
#' @export
#'
#' @examples
#' #FAIRE
cit_gsa <- function(M,
                      X,
                      Z = NULL,
                      test = c("asymptotic","permutation"),
                      n_perm = 100,
                      n_perm_adaptive = c(n_perm, n_perm, n_perm*3, n_perm*5),
                      thresholds = c(0.1,0.05,0.01),
                      parallel = interactive(),
                      n_cpus = detectCores() - 1,
                      adaptive = FALSE,
                      space_y = TRUE,
                      number_y = 10,
                      geneset){
  
  # checks
  if(is.matrix(M)){
    M <- as.data.frame(M)
  }
  
  stopifnot(is.data.frame(M))
  stopifnot(is.data.frame(X))
  stopifnot(is.data.frame(Z) | is.null(Z))
  stopifnot(is.logical(parallel))
  stopifnot(is.logical(adaptive))
  stopifnot(is.numeric(n_perm))
  
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
    
    if(is.vector(geneset)){ # gs vecteur (indices des gènes du gs) ----
      
      if (is.null(Z)){ 
        colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
        modelmat <- as.matrix(model.matrix(~.,data=X))
      } else {
        colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
        colnames(Z) <- sapply(1:ncol(Z), function(i){paste0('Z',i)})
        modelmat <- as.matrix(model.matrix(~.,data=cbind(X,Z)))
      }
      
      indexes_X <- which(substring(colnames(modelmat), 1, 1) == "X")
      
      
      test_stat_gs <- numeric()
      prop_gs <- list()
      indi_pi_gs <- list()

      for (i in 1:length(geneset)){ 
        
        Y <- M[,geneset[i]]
        
        n_Y_all <- length(Y)
        H <- n_Y_all*(solve(crossprod(modelmat)) %*% t(modelmat))[indexes_X, , drop=FALSE] # taille de Y , même pour chaque gène puisque X et Y ne changent pas
        
        # 1) calcule de la stat de test ----
        if (space_y){
          y <- seq(from = ifelse(length(which(Y==0))==0, min(Y), min(Y[-which(Y==0)])),
                   to = max(Y[-which.max(as.matrix(Y))]), length.out = number_y)
        } else{
          y <- sort(unique(Y))
        }
        p <- length(y)
        
        index_jumps <- sapply(y[-p], function(i){sum(Y <= i)}) 
        beta <- c(apply(X = H[, order(Y), drop=FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, ]) / n_Y_all # même nb que seuil
        test_stat <- sum(beta^2) * n_Y_all
        
        test_stat_gs[i] <- test_stat # stat de test pour chaque gène du gs
        
       # 2) calcule de pi ----
        indi_pi <- matrix(NA, n_Y_all, (p-1)) 
        for (j in 1:(p-1)){ 
          indi_Y <- 1*(Y<=y[j])
          indi_pi[,j] <- indi_Y
        }
        indi_pi_gs[[i]] <- indi_pi
        prop <- colMeans(indi_pi)
        prop_gs[[i]] <- prop # prop pour chaque gènes du gs
        
      } 
      
      indi_pi_gs_tab <- do.call(cbind, indi_pi_gs)
      prop_gs_vec <- unlist(prop_gs)

      
      # 3) création de la matrice Sigma ----
      Sigma2 <- matrix(NA,length(prop_gs_vec)*nrow(H),length(prop_gs_vec)*nrow(H)) 
      new_prop <- matrix(NA,length(prop_gs_vec),length(prop_gs_vec))
      
      if (nrow(H)>1){ # > 1 conditions
        
        for (i in 1:nrow(new_prop)){  #  calcule de la nouvelle proportion/du nouveau pi = celle du gene set : ici une matrice
          for(j in 1:ncol(new_prop)){
            
            new_prop[i,j] <- mean((indi_pi_gs_tab[,i]-prop_gs_vec[i]) * (indi_pi_gs_tab[,j]-prop_gs_vec[j])) + prop_gs_vec[i] * prop_gs_vec[j] 
          
            }
        }
        
        Sigma2 <- 1/n * tcrossprod(H) %x%  (new_prop - prop_gs_vec %x%  t(prop_gs_vec))
        
        
      } else { # 1 condition  
        
        for (i in 1:nrow(Sigma2)) {  
          for (j in 1:ncol(Sigma2)) {

            new_prop <- mean((indi_pi_gs_tab[,i]-prop_gs_vec[i]) * (indi_pi_gs_tab[,j]-prop_gs_vec[j])) + prop_gs_vec[i] * prop_gs_vec[j] 
            
            Sigma2[i, j] <- 1/n * tcrossprod(H) * (new_prop - prop_gs_vec[i] * prop_gs_vec[j])
          }
        }
      }
      
      # inutile car boucle créé sigma
      Sigma <- Sigma2*upper.tri(Sigma2, diag = TRUE) +  t(Sigma2*upper.tri(Sigma2, diag = FALSE))
      
      # isSymmetric(Sigma)
      
      
      decomp <- eigen(Sigma, symmetric=TRUE, only.values=TRUE)
      #decomp2 <- svd(Sigma)$d
      
      pval <- survey::pchisqsum(sum(test_stat_gs), lower.tail = FALSE, df = rep(1, ncol(Sigma)),a =decomp$values , method = "saddlepoint") 
      #pval <- survey::pchisqsum(sum(test_stat_gs), lower.tail = FALSE, df = rep(1, ncol(Sigma)),a =svd(Sigma)$d , method = "saddlepoint")
      
      
      # df de résultats
      df <- data.frame(raw_pval=pval,
                       adj_pval =p.adjust(pval, method = "BH"),
                       test_statistic = sum(test_stat_gs))
    
      
          
    } else if (endsWith(geneset, ".gmt") | is(geneset,"BiocSet")){ # gs format gmt ou biocset (liste) ----
      
      if (endsWith(geneset, ".gmt")) {
        geneset <- geneset$genesets
      } else if (is(geneset,"BiocSet")){
        geneset<- BiocSet::es_elementset(geneset)
        geneset_names <- unique(geneset$set)
        geneset <- lapply(X  = unique(geneset$set),
                          FUN = function(x){
                            geneset[geneset$set == x,]$element
                          }
        )
        names(geneset) <- geneset_names
      }
      
      # Colonnes de M doivent avoir le même nom que le nom des gènes dans fichiers gmt ou biocset
      # dans les test j'ai mis exatement e même nombres de gènes dans M (de colonnes) que les fichiers gmt/biocset
      
      if (is.null(Z)){ 
        colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
        modelmat <- as.matrix(model.matrix(~.,data=X))
      } else {
        colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
        colnames(Z) <- sapply(1:ncol(Z), function(i){paste0('Z',i)})
        modelmat <- as.matrix(model.matrix(~.,data=cbind(X,Z)))
      }
      
      indexes_X <- which(substring(colnames(modelmat), 1, 1) == "X")
      
      # initialisation pour chaque gs du gmt
      test_stat_list <- list()
      Sigma2_list <- list()
      decomp_list <- list()
      pval <- NA
      
      
      
      for (k in 1:length(gsa_test)){ # chaque liste de gs du gmt
        
        # initialisation pour chaque gene du gs
        test_stat_gs <- numeric()
        prop_gs <- list()
        indi_pi_gs <- list()
        
        for (i in 1:length(gsa_test[[k]])){ # chaque gène du gs
          
          Y <- M[,gsa_test[[k]][i]]
          
          n_Y_all <- length(Y)
          H <- n_Y_all*(solve(crossprod(modelmat)) %*% t(modelmat))[indexes_X, , drop=FALSE] # taille de Y , même pour chaque gène puisque X et Y ne changent pas
          
          # 1) calcule de la stat de test ----
          if (space_y){
            y <- seq(from = ifelse(length(which(Y==0))==0, min(Y), min(Y[-which(Y==0)])),
                     to = max(Y[-which.max(as.matrix(Y))]), length.out = number_y)
          } else{
            y <- sort(unique(Y))
          }
          p <- length(y)
          
          index_jumps <- sapply(y[-p], function(i){sum(Y <= i)}) 
          beta <- c(apply(X = H[, order(Y), drop=FALSE], MARGIN = 1, FUN = cumsum)[index_jumps, ]) / n_Y_all # même nb que seuil
          test_stat <- sum(beta^2) * n_Y_all
          
          test_stat_gs[i] <- test_stat # stat de test pour chaque gène du gs
          
          
          # 2) calcule de pi ----
          indi_pi <- matrix(NA, n_Y_all, (p-1)) 
          for (j in 1:(p-1)){ 
            indi_Y <- 1*(Y<=y[j])
            indi_pi[,j] <- indi_Y
          }
          indi_pi_gs[[i]] <- indi_pi
          prop <- colMeans(indi_pi)
          prop_gs[[i]] <- prop # prop pour chaque gènes du gs
          
        } 
        test_stat_list[[k]] <- test_stat_gs
        
        indi_pi_gs_tab <- do.call(cbind, indi_pi_gs)
        prop_gs_vec <- unlist(prop_gs)
        
        
        # 3) création de la matrice Sigma ----
        Sigma2 <- matrix(NA,length(prop_gs_vec)*nrow(H),length(prop_gs_vec)*nrow(H)) 
        new_prop <- matrix(NA,length(prop_gs_vec),length(prop_gs_vec))
        
        if (nrow(H)>1){ # > 1 conditions
          
          for (i in 1:nrow(new_prop)){  #  calcule de la nouvelle proportion/du nouveau pi = celle du gene set : ici une matrice
            for(j in 1:ncol(new_prop)){
              
              new_prop[i,j] <- mean((indi_pi_gs_tab[,i]-prop_gs_vec[i]) * (indi_pi_gs_tab[,j]-prop_gs_vec[j])) + prop_gs_vec[i] * prop_gs_vec[j] 
              
            }
          }
          
          Sigma2 <- 1/n * tcrossprod(H) %x%  (new_prop - prop_gs_vec %x%  t(prop_gs_vec))
          
          
        } else { # 1 condition  
          
          for (i in 1:nrow(Sigma2)) {  
            for (j in 1:ncol(Sigma2)) {
              
              new_prop <- mean((indi_pi_gs_tab[,i]-prop_gs_vec[i]) * (indi_pi_gs_tab[,j]-prop_gs_vec[j])) + prop_gs_vec[i] * prop_gs_vec[j] 
              
              Sigma2[i, j] <- 1/n * tcrossprod(H) * (new_prop - prop_gs_vec[i] * prop_gs_vec[j])
            }
          }
        }
        Sigma2_list[[k]] <- Sigma2
        
        decomp_list[[k]] <- eigen(Sigma2_list[[k]], symmetric=TRUE, only.values=TRUE)
        
        pval[k] <- survey::pchisqsum(sum(test_stat_list[[k]]), lower.tail = FALSE, df = rep(1, ncol(Sigma2_list[[k]])),a =decomp_list[[k]]$values , method = "saddlepoint") 
        
      }
      
      df <- data.frame(raw_pval=pval,
                       adj_pval =p.adjust(pval, method = "BH"),
                       test_statistic = unlist(lapply(test_stat_list,sum)))

      
    } 

  } 
    
  
  output <- df
  class(output) <- "cit_gsa"
  return(output)
  
  
  #rownames(df) <- M_colnames
  
  #output <- list(which_test = test,
  #n_perm = n_perm, 
  #pvals = df)
  
  # CHANGER SORTIE SI LISTE !!
  
  
}




















