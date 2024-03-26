############## Fonction générale graph ccdf pour 1 gene set 


# Paramètres ----


# ccdf : object from ccdf function 
# ===> the names of the list elements have to be the names of the genes
# ===> si le passe dans la fonction direct va être trop long, mieux de le faire avant

# ligne 25
# Severity ? comment récup  nom colonne veux : pour la légende pas le choix de demander à un moment : faire plus tard
# ligne 51 : aura pb si X n'est pas une variable à levels donc attention 
# ggplot : légende titre + légende trop récis à l'appli

plot_cit_gsa <- function(ccdf, X, number_y=length(ccdf[[1]]$y)){ #method(median, mean),  space_y=FALSE mettre quand aura parallélisé la fonction, 
  # bof number_y = ccdf
  #pas besoin de Y comme calcule ccdf avant
  

  genes <- names(ccdf)
  
  # Create data frame with all the ccdf ----
  data_gene <- do.call(rbind, lapply(genes, function(name) {
    x <- ccdf[[name]]  
    df <- cbind.data.frame(x$cdf, x$ccdf, x$x, x$y)
    colnames(df) <- c("cdf", "ccdf", "x", "y") # Severity ? comment récup  nom colonne veux : pour la légende pas le choix 
    df$Gene <- rep(name, nrow(df))
    return(df)
  })) 
  
  sev <- unique(data_gene$x)
  
  # Complete data if no values for some y ----
  data_gene_sep <- data_gene%>%group_by(Gene,x)%>%group_split() # separate the table by genes and severity
  
  new_data_gene_sep <- lapply(data_gene_sep, function(df) {
    # max value of current y
    max_y_cur <- max(unique(df$y)) 
    # max value of total y
    max_y_all <- max(unique(data_gene$y))
    
    # all the y that are between max current y value and max total y value
    y_differ_data <- setdiff(data_gene$y, df$y) # y values that are not in current data 
    miss_y <- subset(y_differ_data, y_differ_data >= max_y_cur & y_differ_data <= max_y_all) # y values that are not in current data (previous) and between the 2 max
    n_miss_y <- length(miss_y)

    # add the new y in the data + affect the max value of the ccdf to these y
    # WARNING : ccdf trié pas y !
    if (n_miss_y != 0){
      df <- data.frame(cbind(rep(max(df$cdf),n_miss_y), 
                             rep(max(df$ccdf),n_miss_y), 
                             rep(df$x[1],n_miss_y), miss_y, 
                             rep(df$Gene[1],n_miss_y)))
      colnames(df) <- c("cdf","ccdf","x","y","Gene")
      level <- 
      label <- 
      df$x <- factor(df$x, levels = rep(1:length(levels(data_gene$x))) , labels = levels(data_gene$x) ) #, labels = c("ITU", "NITU", "Mild")) # AURA PB si x n'est pas une variable avec des levels 
    } else{
      df <- 0
    }
    return(df)
  })
  
  # Combine the data set
  final_data <- do.call(rbind,lapply(gene, function(name){
    # if some genes doesn't have missing y values
    no_zero <- do.call(rbind,Filter(function(x) {!is.numeric(x) || !all(x == 0)}, new_data_gene_sep)) 
    # new dataset completed
    df <- rbind(data_gene[data_gene$Gene==name,], no_zero[no_zero$Gene==name,])
    df$y <- as.numeric(df$y)
    df$cdf <- as.numeric(df$cdf)
    df$ccdf <- as.numeric(df$ccdf)
    return(df)
  }))
  
  # Thresholds ----
  Y_after <- final_data$y
  # y value for each thresholds  
  y_after <- seq(from = ifelse(length(which(Y_after==0))==0, min(Y_after), min(Y_after[-which(Y_after==0)])), 
                 to = max(Y_after[-which.max(as.matrix(Y_after))]), length.out = number_y) 
  #p_after <- length(y_after)
  # index thresholds 
  #index_jumps_after <- sapply(y_after[-p_after], function(i){sum(Y_after <= i)})  
  
  seuils <- c(0,y_after)
  
  
  # Compute the median ----
  med <- lapply(sev, function(x){
    med <- rep(NA, number_y)
    filtre_x <- final_data$x == x
    filtre_row_i0 <-  final_data$y >= seuils[1] & final_data$y < seuils[1+1]
    indices0 <- which(filtre_x & filtre_row_i0)
    med[1] <- median(final_data$ccdf[indices0])
    
    for (i in 2:length(seuils)){
      
      filtre_row_i <-  final_data$y >= seuils[i] & final_data$y < seuils[i+1]
      indices <- which(filtre_x & filtre_row_i)
      
      med[i] <- median(final_data$ccdf[indices])
      
      ref_value <- med[i-1]
      range_value <- final_data[filtre_x & filtre_row_i, ]$ccdf
      
      closest_index <- min(which(range_value > ref_value))
      closest_value <- range_value[closest_index]
      
      if (med[i] < med[i-1] & closest_index != Inf ){ # force monotony when we can, if no values > previous median leave the current median 
        med[i] <- closest_value
      }
    }
    med <- data.frame(med[1:number_y])
    med$y <- y_after 
    med$x <- x 
    names(med)[1] <- "ccdf"
    
    return(med)
  })
  browser()
  # Replace values with the max : if some values are below the max after it 
  replace_max <- lapply(med, function(m){
    ind_max <- which.max(m$ccdf)
    m[ind_max:nrow(m),]$ccdf <- max(m$ccdf)
    return(m)
  })
  
  # Step function : add the ccdf to the y values that are not a thresholds
  step_function <- lapply(replace_max, function(df){
    # y values that are not in the data of the thresholds
    new_y <- data.frame(setdiff(final_data$y, df$y)) ; colnames(new_y) <- "y" 
    # separate the values between the intervals of the thresholds
    intervalle_indices <- findInterval(new_y$y, df$y, left.open = TRUE) 
    indices_int <- intervalle_indices+1 
    
    new_y$ccdf <- df$ccdf[indices_int]
    new_y$x <- unique(df$x)
    
    return(new_y)
  })
  
  
  # Data final -----
  data_med <- rbind(do.call(rbind,step_function),do.call(rbind,med))
  
  
  
  # Data modified for the legend
  data_med$Gene <- "all"
  data_med$Legend <- "Gene set summary"  
  
  data_gene <- data_gene[-1]
  data_gene$Legend <- "Genes" 
  
  
  comb_data <- rbind(data_gene,data_med)
  
  
  # Plot -----
  ggplot(data = comb_data, aes(x = y, color=x, y = ccdf,linetype = Legend, size=Legend)) +
    geom_line(aes(group = interaction(Gene,x))) + 
    scale_size_manual(values = c(0.8,0.2)) + # change taille lignes selon variable dans aes(size) 
    #scale_color_manual(values = c("ITU" = "brown3", "NITU" = "darkgoldenrod1", "Mild"="chartreuse3")) + # pécis pour ces données
    labs(x = "Gene expression") +
    theme_minimal() 
  
}



plot_cit_gsa(ccdf_gs2_X_space, X = data.frame(X=factor(covid,levels=c("ITU","NITU","Mild"))[keep_cluster]) , number_y=60)

# faire broser pour voir pourquoi tombe 














