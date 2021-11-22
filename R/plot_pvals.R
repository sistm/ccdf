#'Plot of gene-wise p-values
#'
#'This function prints the sorted exact p-values along with the Benjamini-Hochberg
#'limit and the 5% threshold
#'
#'@param pvals a vector of length \code{n} containing the raw p-values for
#'each gene
#'
#'@return a plot of sorted gene-wise p-values
#'@import viridisLite
#'@import ggplot2
#'
#'@return a \code{\link[ggplot2]{ggplot}} object
#'
#'@export
#'
#'@examples
#' 
#'plot_pvals(runif(100,0,1))


plot_pvals <- function(pvals){
  t <- 1:length(pvals)
  s <- (t/length(pvals))*0.05
  df_plot_perm <- data.frame("y" = sort(pvals), "x" = 1:length(pvals))
  ggplot()+ scale_y_log10()+
    geom_point(data = df_plot_perm,aes_string(x = "x", y = "y", color = shQuote(viridis(4)[1])), size = 0.5)+
    geom_line(data = df_plot_perm, aes_string(y = "s", x = "x",color = shQuote(viridis(4)[2])), size = 0.5) +
    geom_line(data = df_plot_perm, aes_string(y = 0.05, x = "x", color = shQuote("red")), size = 0.5) +
    scale_color_manual(name = "", labels = c("B-H limit", "p-values", "5% threshold"),
                       values = c(viridis(4)[c(1,2)], "red")) + xlab("rank") +
    ylab("log10 scale") + xlim(0, length(df_plot_perm$y)) +
    theme_bw()
}
