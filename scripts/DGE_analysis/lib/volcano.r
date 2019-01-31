library(ggplot2)
library(latex2exp)

volcano_plot <- function(lfc, pvals, title = NULL) {
  alpha <- 0.05
  cols <- densCols(lfc, -log10(pvals))
  p.size <- abs(lfc)
  p.size <- 3*p.size/max(p.size)
  
  plot.obj <- ggplot(data.frame(lfc = lfc, pvals = pvals), aes(x = lfc, y = -log10(pvals))) + 
              geom_point(aes(color = -log10(pvals)), size = p.size) + 
              scale_color_gradient(low = "black", high = "red") + 
              geom_vline(xintercept = c(-1,1), alpha = 0.3) + 
              geom_hline(yintercept = -log10(alpha),alpha = 0.3) + 
              xlab(TeX("$log_2(foldChange)$")) + 
              ylab(TeX("$-log_{10}(pval)$")) +
              xlim(c(-5,5)) + 
              labs(color = TeX("-log_{10}($p_{adj})"), size = TeX("$log_2(foldChange)$")) 
  
  return(plot.obj)
  
}

