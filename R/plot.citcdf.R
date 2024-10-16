
#' Function for plotting the CDF of all the genes within a gene set or for one gene
#'
#' @param x an object of class \code{\link{cit_gsa}} 
#' 
#' @param ... further arguments to be passed
#'
#' @return a \code{ggplot2} of gene-wise or geneset ccdf
#' 
#'@importFrom viridisLite viridis
#'@import ggplot2
#'
#' @export
#'
#' @examples
#' #TO DO




plot.citcdf <- function(x, ...){ 
  
  if(x$type == "gsa"){
  
    ccdf_all <- x$ccdf
    
    browser()
    
    #-------------
    # Case where we have Z !! TO DO : if we don't have Z (change in cot_gsa function?)
    # Case X qualitatif (values in "sev")
    #-------------
    
    for (i in 1:length(ccdf_all)){
      ccdf <- ccdf_all[[i]]
      
      genes <- names(ccdf)
      
      # test_ccdf[[1]][[1]] # accède à Y1 Y2 names(ccdf_all[[1]][[1]])
      # test_ccdf[[2]][[1]] # accède à Y3 Y4 names(ccdf_all[[2]][[1]])
      
      # Create data frame with all the ccdf ----
      data_gene <- do.call(rbind, lapply(genes, function(name) {
        x <- ccdf[[name]]  
        df <- cbind.data.frame(x$cdf, x$ccdf_x, x$x, x$y) # not ccdf_nox et z, METTRE ?
        colnames(df) <- c("cdf", "ccdf_x", "x", "y") 
        df$Gene <- rep(name, nrow(df))
        return(df)
      })) 
      
      sev <- unique(data_gene$x)
      
      
      # Complete data if no values for some y ----
      data_gene_sep <- split(data_gene, list(data_gene$Gene, data_gene$x)) # separate the table by genes and x
      
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
        # WARNING : ccdf sorted not y !
        if (n_miss_y != 0){
          df <- data.frame(cbind(rep(max(df$cdf),n_miss_y), 
                                 rep(max(df$ccdf_x),n_miss_y), 
                                 rep(df$x[1],n_miss_y), miss_y, 
                                 rep(df$Gene[1],n_miss_y)))
          colnames(df) <- c("cdf","ccdf_x","x","y","Gene")
          df$x <- factor(df$x, levels = rep(1:length(levels(data_gene$x))) , labels = levels(data_gene$x) ) 
        } else{
          df <- 0
        }
        return(df)
      })
      
      # Combine the data set
      final_data <- do.call(rbind,lapply(genes, function(name){
        # if some genes doesn't have missing y values
        no_zero <- do.call(rbind,Filter(function(x) {!is.numeric(x) || !all(x == 0)}, new_data_gene_sep)) 
        # new dataset completed
        df <- rbind(data_gene[data_gene$Gene==name,], no_zero[no_zero$Gene==name,])
        df$y <- as.numeric(df$y)
        df$cdf <- as.numeric(df$cdf)
        df$ccdf_x <- as.numeric(df$ccdf_x)
        return(df)
      }))
      
      # Thresholds ----
      Y_after <- final_data$y
      number_y = 20
      # y value for each thresholds  
      y_after <- seq(from = ifelse(length(which(Y_after==0))==0, min(Y_after), min(Y_after[-which(Y_after==0)])), 
                     to = max(Y_after[-which.max(as.matrix(Y_after))]), length.out = number_y) 
      
      seuils <- c(0,y_after)
      
      # Compute the median ----
      med <- lapply(sev, function(x){
        med <- rep(NA, number_y)
        filtre_x <- final_data$x == x
        filtre_row_i0 <-  final_data$y >= seuils[1] & final_data$y < seuils[1+1]
        indices0 <- which(filtre_x & filtre_row_i0)
        med[1] <- median(final_data$ccdf_x[indices0])
        
        for (i in 2:length(seuils)){
          
          filtre_row_i <-  final_data$y >= seuils[i] & final_data$y < seuils[i+1]
          indices <- which(filtre_x & filtre_row_i)
          
          med[i] <- median(final_data$ccdf_x[indices])
          
          ref_value <- med[i-1]
          range_value <- final_data[filtre_x & filtre_row_i, ]$ccdf_x
          
          closest_index <- min(which(range_value > ref_value))
          closest_value <- range_value[closest_index]
          
          if (med[i] < med[i-1] & closest_index != Inf ){ # force monotony when we can, if no values > previous median leave the current median 
            med[i] <- closest_value
          }
        }
        med <- data.frame(med[1:number_y])
        med$y <- y_after 
        med$x <- x 
        names(med)[1] <- "ccdf_x"
        
        return(med)
      })
      
      # Replace values with the max : if some values are below the max after it 
      replace_max <- lapply(med, function(m){
        ind_max <- which.max(m$ccdf_x)
        m[ind_max:nrow(m),]$ccdf_x <- max(m$ccdf_x)
        return(m)
      })
      
      # Step function : add the ccdf to the y values that are not a thresholds
      step_function <- lapply(replace_max, function(df){
        # y values that are not in the data of the thresholds
        new_y <- data.frame(setdiff(final_data$y, df$y)) ; colnames(new_y) <- "y" 
        # separate the values between the intervals of the thresholds
        intervalle_indices <- findInterval(new_y$y, df$y, left.open = TRUE) 
        indices_int <- intervalle_indices+1 
        
        new_y$ccdf_x <- df$ccdf_x[indices_int]
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
      p <- ggplot(data = comb_data, aes(x = y, color=x, y = ccdf_x,linetype = Legend, size=Legend)) +
        geom_line(aes(group = interaction(Gene,x))) + 
        scale_size_manual(values = c(0.8,0.2)) + # modify row size according to the variable in aes(size) 
        labs(x = "Gene expression") +
        theme_minimal()
      
      return(p)
    }


  } 
  

  # elseif (x$type == genewise){
  
  #}
  
  
}












