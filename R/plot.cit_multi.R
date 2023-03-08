#'Plot of gene-wise p-values
#'
#'Plotting the sorted exact p-values along with the Benjamini-Hochberg
#'limit and the nominal threshold
#'
#'@param x an object of class \code{\link{cit_multi}}.
#'
#'@param ... further arguments to be passed
#'
#'@param nominal_level a nominal testing level between 0 and 1. 
#'Default is 5\%: \code{0.05}.
#'
#'@return a \code{ggplot2} of sorted gene-wise p-values
#'
#'@importFrom viridisLite viridis
#'@import ggplot2
#'
#'@export
#'
#'@examples
#'ncells <- 100
#'pgenes <- 500
#'X1 <- as.factor(rbinom(n=ncells, size = 1, prob = 0.5))
#'Y <- t(replicate(pgenes, ((X1==1)*rnorm(n = ncells,0,1)) + ((X1==0)*rnorm(n = ncells, 0.5, 1))))
#'
#'res_asymp <- cit_multi(exprmat=data.frame(Y=Y),
#'                       variable2test=data.frame(X1), 
#'                       test="asymptotic", n_cpus=1, parallel=FALSE)
#'plot(res_asymp)
#'
plot.cit_multi <- function(x, ..., nominal_level = 0.05){
  
  stopifnot(isa(x, "cit_multi"))
  
  pvals <- x$pvals$raw_pval
  
  n <- length(pvals)
  t <- 1:n
  s <- t/n * nominal_level
  
  df_plot_perm <- data.frame("y" = sort(pvals), "x" = t, "s" = s)
  
  p <- ggplot()+ scale_y_log10() + ggplot2::annotation_logticks(sides="l") +
    geom_point(data = df_plot_perm, aes(x = .data[["x"]], y = .data[["y"]], color = "pvals"), size = 0.5)+
    geom_line(data = df_plot_perm, aes(y = .data[["s"]], x = .data[["x"]], color = "bhlim"), linewidth = 0.5) +
    geom_line(data = df_plot_perm, aes(y = 0.05, x = .data[["x"]], color = "nomlev"), linewidth = 0.5) +
    scale_color_manual(name = "", breaks=c("pvals", "bhlim", "nomlev") ,
                       labels = c( "Raw p-values", "Bonferroni-Hochberg threshold",
                                  paste0(round(nominal_level*100, digits = 0), "% nominal level")),
                       values = c(viridis(4)[c(2,1)], "red")) +
    guides(color=guide_legend(override.aes = list(shape = c(16,NA,NA), size = c(1,NA,NA), linewidth=c(NA,0.5,0.5)))) +
    xlab("Rank (ascending order)") +
    ylab("p-values (log10 scale)") + xlim(0, length(df_plot_perm$y)) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  return(p)
}

