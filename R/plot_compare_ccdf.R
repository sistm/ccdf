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
#'@param method a character string indicating which method to use to
#'compute the CCDF, either \code{'OLS'} or \code{'logistic'}.
#'Default is \code{'OLS'} for computational speed.
#'
#'@param fast a logical flag indicating whether the fast implementation of
#'logistic regression should be used. Only if \code{method == 'logistic'}.
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
#'plot_compare_ccdf(data.frame(Y=Y),data.frame(X=X),method="OLS")


plot_compare_ccdf <- function(Y, X, Z=NULL, method=c("OLS","logistic"), 
                              fast=TRUE, space_y=FALSE, number_y=length(Y)){
  
  res <- ccdf(as.numeric(Y[,1]), X, Z, method, fast, space_y, number_y)
  
  
  if (is.null(Z)){ # Z null ----
    
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
                           guide = guide_legend(override.aes = list(linetype = c("dotted",rep("solid",l_X))))) + xlab(colnames(Y)) +
        ylab("Cumulative probability") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    }
    else{
      ggplot() + ggtitle(colnames(Y)) +
        geom_step(data = df_plot,aes_string(x = "y", y = "cdf", color = shQuote(viridis(n=2)[1])), size = 0.5, linetype="dotted") +
        geom_point(data = df_plot, aes_string(x = "y", y = "ccdf", color =  shQuote(viridis(n=2)[2])), size = 0.5) +
        scale_color_manual(name = "", labels = c("CDF", "CCDF"),
                           values = c(viridis(n=3)),
                           guide = guide_legend(override.aes = list(linetype = c("dotted","blank"),
                                                                    shape = c(NA,16)))) + xlab(colnames(Y)) +
        ylab("Cumulative probability") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    }
    
    
  }else{ # Z specified ----
    
    df_plot <- data.frame("y" = res$y, "x" = res$x, "z" = res$z, "cdf" = res$cdf, 
                          "ccdf_nox" = res$ccdf_nox,  "ccdf_x" = res$ccdf_x, 
                          "shape1"= "20", "shape2" = "3")
    
    res_X <- ccdf(as.numeric(Y[,1]), X, Z=NULL, method, fast, space_y,number_y)
    df_X <- data.frame("y" = res_X$y, "x" = res_X$x, "cdf" = res_X$cdf, "ccdf" = res_X$ccdf)
    
    if (is.factor(Z[,1])){ # Z factor
      if (is.factor(X[,1])){ # X factor
        
        levels(df_X$x) <- c(levels(X[,1]), "CDF")
        df_X$x <- ordered(df_X$x, levels = levels(df_X$x))
        l_X <- nlevels(X[,1])
        
        p_cdf <- ggplot() +
          geom_step(data = df_X, aes(x = .data[["y"]], y = .data[["cdf"]], color = "CDF"), size = 0.5, linetype="dotted") +
          geom_step(data = df_X, aes(x = .data[["y"]], y = .data[["ccdf"]], color = .data[["x"]]), size = 0.5) +
          scale_color_manual(name = "", limits = levels(df_X$x),
                             values = c(viridis(n=(l_X+1))[-(l_X+1)],"#CC0066"),
                             guide = guide_legend(override.aes = list(linetype = c(rep("solid",l_X),"dotted")))) + 
          xlab(colnames(Y)) +
          ylab("Cumulative probability") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
        
        # Z
        levels(df_plot$x) <- c(levels(X[,1]), "Marginal on Z")
        df_plot$x <- ordered(df_plot$x, levels = levels(df_plot$x))
        l_X <- length(unique(X[,1]))
        l_Z <- length(unique(Z[,1]))
        
        p_facet <- ggplot() + ggtitle(colnames(X)) +
          geom_step(data = df_plot, aes(x = .data[["y"]], y = .data[["ccdf_x"]], color = .data[["x"]])) +
          geom_step(data = df_plot, aes(x = .data[["y"]], y = .data[["ccdf_nox"]], color = shQuote("Marginal on Z")), size = 0.5, linetype = "dotted") +
          scale_color_manual(name = "CCDF", 
                             values = c(viridis(n=(l_X+1))[-(l_X+1)],"gold"), limits = levels(df_plot$x),
                             guide = guide_legend(override.aes = list(linetype = c(rep("solid",l_X),"dotted")))) +
          xlab(colnames(Y)) +
          ylab("Cumulative probability") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(.~as.factor(df_plot$z), ncol=3, scales = "fixed")
        
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
          xlab(colnames(Y)) +
          ylab("Cumulative probability") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
        
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
          xlab(colnames(Y)) +
          ylab("Cumulative probability") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
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
                             guide = guide_legend(override.aes = list(linetype = c(rep("solid",l_X),"dotted")))) + xlab(colnames(Y)) +
          ylab("Cumulative probability") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
        
        levels(df_plot$x) <- c(levels(X[,1]), "Marginal on Z")
        df_plot$x <- ordered(df_plot$x, levels = levels(df_plot$x))
        l_X <- length(unique(X[,1]))
        
        p_facet <- ggplot() + ggtitle(colnames(Y)) +
          geom_point(data = df_plot, aes_string(x = "y", y = "ccdf_x", color = "x"), shape=16, size = 0.5) +
          geom_point(data = df_plot, aes_string(x = "y", y = "ccdf_nox", color = shQuote("Marginal on Z")),  shape=16, size = 0.5) +
          scale_color_manual(name = "CCDF",
                             values = c(viridis(n=4)[-4],"gold"), limits = levels(df_plot$x),
                             guide = guide_legend(override.aes = list(linetype = c(rep("solid",l_X),"dotted")))) +
          xlab(colnames(Y)) +
          ylab("Cumulative probability") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
        
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
          xlab(colnames(Y)) +
          ylab("Cumulative probability") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
      }
      
    }
  }
}

