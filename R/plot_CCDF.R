#' Function to plot the CCDF according to the type of X et Z
#'
#'@param Y a numeric vector of size \code{n} containing the
#'preprocessed expressions from \code{n} samples (or cells).
#'
#'@param X a numeric or factor vector of size \code{n}
#'containing the variable to be tested (the condition to be tested). 
#' 
#'@param Z a numeric or factor vector of size \code{n}
#'containing the covariate. Multiple variables are not allowed.
#'
#' @param sample_group a vector of length \code{n} indicating whether the samples
#'should be grouped (e.g. paired samples or longitudinal data). Coerced
#'to be a \code{factor}. Default is \code{NULL} in which case no grouping is
#'performed.
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
#'@param space_y a logical flag indicating whether the y thresholds are spaced. 
#'When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and 
#'the maximum of the observations is used. Default is \code{FALSE}.
#'
#'@param number_y an integer value indicating the number of y thresholds (and therefore
#'the number of regressions) to perform the test. Default is \code{length(Y)}.
#'
#'@import ggplot2 cowplot stats
#'
#'@return a \code{\link[ggplot2]{ggplot}} object
#'
#'@export
#'
#'@examples
#' 
#'X <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
#'Y <- ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,0.5,1))
#'plot_CCDF(data.frame(Y=Y),data.frame(X=X),method="linear regression")


plot_CCDF <- function(Y,X,Z=NULL,sample_group=NULL,method="linear regression",fast=TRUE,space_y=FALSE,number_y=length(Y)){
  
  
  stopifnot(is.data.frame(Y))
  stopifnot(is.data.frame(X))
  stopifnot(is.data.frame(Z) | is.null(Z))
  stopifnot(is.logical(fast))
  
  if (sum(is.na(Y)) > 1) {
    warning("'y' contains", sum(is.na(Y)), "NA values. ",
            "\nCurrently they are ignored in the computations but ",
            "you should think carefully about where do those NA/NaN ",
            "come from...")
    Y <- Y[stats::complete.cases(Y)]
  }
  
  
  if (length(method) > 1) {
    method <- method[1]
  }
  stopifnot(method %in% c("linear regression","logistic regression", "CART", "RF", "mixed model"))

  if (space_y){
    if (is.null(number_y)){
      warning("Missing argument", number_y, ". No spacing is used.")
      space_y <- FALSE
    }
  }

  res <- CCDF(as.numeric(Y[,1]),X,Z,sample_group = sample_group,method,fast,space_y,number_y)


  if (class(Z[,1])=="NULL"){
    if (method == "mixed model" & !is.null(sample_group))
      df_plot <- data.frame("y" = res$y, "x" = res$x, "sample_group" = res$sample_group,"cdf" = res$cdf, "ccdf" = res$ccdf)
    else
      df_plot <- data.frame("y" = res$y, "x" = res$x, "cdf" = res$cdf, "ccdf" = res$ccdf)

    if (is.factor(X[,1])){
      levels(df_plot$x) <- unique(X[,1])
      df_plot$x <- ordered(df_plot$x, levels = levels(df_plot$x))
      l_X <- length(unique(X[,1]))

      ggplot() + ggtitle(colnames(Y)) +
        geom_step(data = df_plot, aes_string(x = "y", y = "cdf", color = shQuote(viridis(n=(l_X+1))[1])), size = 0.5, linetype="dotted") +
        geom_step(data = df_plot, aes_string(x = "y", y = "ccdf", color = "x"), size = 0.5) +
        scale_color_manual(name = "", labels=c("CDF", paste0("CCDF X=", levels(df_plot$x)[ordered(levels(df_plot$x))])),
                           values = viridis(n=l_X+1),
                           guide = guide_legend(override.aes = list(linetype = c("dotted",rep("solid",l_X))))) + xlab("gene expression") +
        ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    }
    else{
      ggplot() + ggtitle(colnames(Y)) +
        geom_step(data = df_plot,aes_string(x = "y", y = "cdf", color = shQuote(viridis(n=2)[1])), size = 0.5, linetype="dotted") +
        geom_point(data = df_plot, aes_string(x = "y", y = "ccdf", color =  shQuote(viridis(n=2)[2])), size = 0.5) +
        scale_color_manual(name = "", labels = c("CDF", "CCDF"),
                           values = c(viridis(n=3)),
                           guide = guide_legend(override.aes = list(linetype = c("dotted","blank"),
                                                                    shape = c(NA,16)))) + xlab("gene expression") +
        ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    }
  }
  
  
  
  else{ # Z is specified
    
    res <- CCDF(as.numeric(Y[,1]),X,Z,sample_group,method,fast,space_y,number_y)
    if(method == "mixed model")
      df_plot <- data.frame("y" = res$y, "x" = res$x, "z" = res$z, "sample_group" = res$sample_group, "cdf" = res$cdf, "ccdf_nox" = res$ccdf_nox,  "ccdf_x" = res$ccdf_x, "shape1"= "20", "shape2" = "3")
    else
      df_plot <- data.frame("y" = res$y, "x" = res$x, "z" = res$z, "cdf" = res$cdf, "ccdf_nox" = res$ccdf_nox,  "ccdf_x" = res$ccdf_x, "shape1"= "20", "shape2" = "3")
    
    res_X <- CCDF(as.numeric(Y[,1]),X,Z=NULL,sample_group,method,fast,space_y,number_y)
    df_X <- data.frame("y" = res_X$y, "x" = res_X$x, "cdf" = res_X$cdf, "ccdf" = res_X$ccdf)
    
    if (is.factor(Z[,1])){ # Z factor
      if (is.factor(X[,1])){ # X factor
        
        levels(df_X$x) <- c(levels(X[,1]), "CDF")
        df_X$x <- ordered(df_X$x, levels = levels(df_X$x))
        l_X <- length(unique(X[,1]))
        
        p_cdf <- ggplot() + ggtitle(colnames(Y)) +
          geom_step(data = df_X, aes_string(x = "y", y = "cdf", color = shQuote("CDF")), size = 0.5, linetype="dotted") +
          geom_step(data = df_X, aes_string(x = "y", y = "ccdf", color = "x"), size = 0.5) +
          scale_color_manual(name = "", limits = levels(df_X$x),
                             values = c(viridis(n=(l_X+1))[-(l_X+1)],"#CC0066"),
                             guide = guide_legend(override.aes = list(linetype = c(rep("solid",l_X),"dotted")))) + xlab("gene expression") +
          ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
        
        # Z
        levels(df_plot$x) <- c(levels(X[,1]), "Marginal on Z")
        df_plot$x <- ordered(df_plot$x, levels = levels(df_plot$x))
        l_X <- length(unique(X[,1]))
        l_Z <- length(unique(Z[,1]))
        
        p_facet <- ggplot() + ggtitle(colnames(Y)) +
          geom_step(data = df_plot, aes_string(x = "y", y = "ccdf_x", color = "x")) +
          geom_step(data = df_plot, aes_string(x = "y", y = "ccdf_nox", color = shQuote("Marginal on Z")), size = 0.5, linetype = "dotted") +
          scale_color_manual(name = "CCDF", 
                             values = c(viridis(n=(l_X+1))[-(l_X+1)],"gold"), limits = levels(df_plot$x),
                             guide = guide_legend(override.aes = list(linetype = c("solid","solid","solid","dotted")))) +
          xlab("gene expression") +
          ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(.~as.factor(df_plot$z), ncol=3, scales = "fixed")
        
        ggdraw() +
          draw_plot(p_cdf, 0, 0.75, 0.75, 0.25) +
          draw_plot(p_facet, 0, 0, 1, 0.75) +
          draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.75), size = 15)
        
      }
      else{ # X continuous
        p_cdf <- ggplot() + ggtitle(colnames(Y)) +
          geom_step(data = df_X, aes_string(x = "y", y = "cdf", color=shQuote("CDF"))) +
          geom_point(data = df_X, aes_string(x = "y", y = "ccdf", color=shQuote("CCDF")), shape=16, size = 0.5) +
          scale_color_manual(name = "",
                             labels = c("CCDF", "CDF"),
                             values = c("#CC0066",viridis(n=3)[2]),
                             guide = guide_legend(override.aes = list(linetype = c("blank","solid"),
                                                                      shape = c(16,NA)))) +
          xlab("gene expression") +
          ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
        
        levels(df_plot$x) <- c(levels(X[,1]), "Marginal on Z")
        df_plot$x <- ordered(df_plot$x, levels = levels(df_plot$x))
        l_X <- length(unique(X[,1]))
        l_Z <- length(unique(Z[,1]))
        
        p_facet <- ggplot() + ggtitle(colnames(Y)) +
          geom_point(data = df_plot, aes_string(x = "y", y = "ccdf_x", color=shQuote(viridis(n=3)[1])), shape=16, size = 0.5) +
          geom_point(data = df_plot, aes_string(x = "y", y = "ccdf_nox", color=shQuote('gold')), shape=2, size = 0.5) +
          scale_color_manual(name = "CCDF",
                             labels = c("Given X and Z", "Marginal on Z"),
                             values = c(viridis(n=3)[1],'gold'),
                             guide = guide_legend(override.aes = list(linetype = c("blank","blank"),
                                                                      shape = c(16,2)))) +
          xlab("gene expression") +
          ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
          facet_wrap(.~as.factor(df_plot$z), nrow=3, scales = "fixed")
        
        
        ggdraw() +
          draw_plot(p_cdf, 0, 0.75, 0.75, 0.25) +
          draw_plot(p_facet, 0, 0, 1, 0.75) +
          draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.75), size = 15)
        
      }
      
    }
    
    else{ # Z continuous
      
      if (is.factor(X[,1])){ # X factor
        levels(df_X$x) <- c(levels(X[,1]), "CDF")
        df_X$x <- ordered(df_X$x, levels = levels(df_X$x))
        l_X <- length(unique(X[,1]))
        
        p_cdf <- ggplot() + ggtitle(colnames(Y)) +
          geom_step(data = df_X, aes_string(x = "y", y = "cdf", color = shQuote("CDF")), size = 0.5, linetype="dotted") +
          geom_step(data = df_X, aes_string(x = "y", y = "ccdf", color = "x")) +
          scale_color_manual(name = "", limits = levels(df_X$x), #labels=c("CDF", paste0("CCDF X=", levels(df_X$x)[ordered(levels(df_X$x))])),
                             values = c(viridis(n=4)[-4],"#CC0066"),
                             guide = guide_legend(override.aes = list(linetype = c(rep("solid",l_X),"dotted")))) + xlab("gene expression") +
          ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
        
        levels(df_plot$x) <- c(levels(X[,1]), "Marginal on Z")
        df_plot$x <- ordered(df_plot$x, levels = levels(df_plot$x))
        l_X <- length(unique(X[,1]))
        
        p_facet <- ggplot() + ggtitle(colnames(Y)) +
          geom_point(data = df_plot, aes_string(x = "y", y = "ccdf_x", color = "x"), shape=16, size = 0.5) +
          geom_point(data = df_plot, aes_string(x = "y", y = "ccdf_nox", color = shQuote("Marginal on Z")),  shape=16, size = 0.5) +
          scale_color_manual(name = "CCDF",
                             values = c(viridis(n=4)[-4],"gold"), limits = levels(df_plot$x),
                             guide = guide_legend(override.aes = list(linetype = c(rep("solid",l_X),"dotted")))) +
          xlab("gene expression") +
          ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
        
        ggdraw() +
          draw_plot(p_cdf, 0, 0.75, 0.75, 0.25) +
          draw_plot(p_facet, 0, 0, 1, 0.75) +
          draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.75), size = 15)
      }
      
      else{
        ggplot() + ggtitle(colnames(Y)) +
          geom_step(data = df_plot, aes_string(x = "y", y = "cdf",  color=shQuote(viridis(n=3)[2])), size = 0.7) +
          geom_point(data = df_plot, aes_string(x = "y", y = "ccdf_x", color=shQuote(viridis(n=3)[1])), shape=16, size = 0.5) +
          geom_point(data = df_plot, aes_string(x = "y", y = "ccdf_nox", color=shQuote('gold')), shape=2, size = 0.5) +
          scale_color_manual(name = "",
                             labels = c("CDF", "CCDF_X", "CCDF_noX"),
                             values = c(viridis(n=3)),
                             guide = guide_legend(override.aes = list(linetype = c("solid","blank","blank"),
                                                                      shape = c(NA,16,2)))) +
          xlab("gene expression") +
          ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
      }
      
    }
  }
}

