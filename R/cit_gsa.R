
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
                    genesets){
  
  # check
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
    # res : les p-val + stat de test pour chaque gènes du geneset
    # df : p-val + stat de test pour tous le geneset
    
    if(is.vector(genesets)){
      res <- do.call("rbind", lapply(1:length(genesets),function(j){
        cit_asymp(M[, genesets[j]], X, Z, # prend la colonne du gène dans M qui correspond au nom du gene de genesets
                  space_y = space_y, 
                  number_y = number_y)
      }))
      
      # p-val ----
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
      
      for (i in 1:length(geneset)){
        Y <- M[,geneset[i]]
        
        n_Y_all <- length(Y)
        H <- n_Y_all*(solve(crossprod(modelmat)) %*% t(modelmat))[indexes_X, , drop=FALSE] # taille de Y , même pour chaque gène puisque X et Y ne changent pas
        number_y = length(unique(Y))
        
        
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
        
        
        indi_pi <- matrix(NA, n_Y_all, (p-1)) 
        for (j in 1:(p-1)){ 
          indi_Y <- 1*(Y<=y[j])
          indi_pi[,j] <- indi_Y
        }
        prop <- colMeans(indi_pi)
        prop_gs[[i]] <- prop # prop pour chaque gènes du gs
        
      } 
      
      prop_gs_vec <- unlist(prop_gs)
      
      Sigma2 <- 1/n* tcrossprod(H) %x% (prop_gs_vec - prop_gs_vec %x% t(prop_gs_vec)) 
      # utilise les proportion mise à la suite sous forme de vecteur
      Sigma <- Sigma2*upper.tri(Sigma2, diag = TRUE) +  t(Sigma2*upper.tri(Sigma2, diag = FALSE))
      
      decomp <- eigen(Sigma, symmetric=TRUE, only.values=TRUE)
      
      pval <- survey::pchisqsum(sum(test_stat_gs), lower.tail = FALSE, df = rep(1, ncol(Sigma)), #ncol(Sigma)
                                a = decomp$values, method = "saddlepoint")
      # ncol(Sigma= nombre de beta) : tester ravec length(decomp) : ne fonctione pas
      #/ nombre beta total (ne fonctionne pas  = ncol(sigma))
      # / nombre beta pour 1 gene ?
      
      # df de résultats
      df <- data.frame(raw_pval=pval,
                       adj_pval =p.adjust(pval, method = "BH"),
                       test_statistic = sum(test_stat_gs))
      
      
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




















