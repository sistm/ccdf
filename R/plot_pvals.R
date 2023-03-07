#'Plot of gene-wise p-values
#'
#'Plotting the sorted exact p-values along with the Benjamini-Hochberg
#'limit and the nominal threshold
#'
#'@param pvals a vector of length \code{n} containing the raw p-values for
#'each gene
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
#'plot_pvals(runif(100,0,1))
plot_pvals <- function(pvals, nominal_level = 0.05){
  
  n <- length(pvals)
  t <- 1:n
  s <- t/n * nominal_level
  
  df_plot_perm <- data.frame("y" = sort(pvals), "x" = t, "s" = s)
  
  p <- ggplot()+ scale_y_log10()+
    geom_point(data = df_plot_perm, aes(x = .data[["x"]], y = .data[["y"]], color = "pvals"), size = 0.5)+
    geom_line(data = df_plot_perm, aes(y = .data[["s"]], x = .data[["x"]], color = "bhlim"), linewidth = 0.5) +
    geom_line(data = df_plot_perm, aes(y = 0.05, x = .data[["x"]], color = "nomlev"), linewidth = 0.5) +
    scale_color_manual(name = "", breaks=c("bhlim", "pvals", "nomlev") ,
                       labels = c("Bonferroni-Hochberg threshold", "Raw p-values", 
                                  paste0(round(nominal_level*100, digits = 0), "% nominal level")),
                       values = c(viridis(4)[c(1,2)], "red")) + 
    xlab("Rank") +
    ylab("P-values (log10 scale)") + xlim(0, length(df_plot_perm$y)) +
    theme_bw()
  
  return(p)
}
