library(ggplot2)
library(latex2exp)

volcano_plot <- function(deseq_res, title = NULL) {
  alpha <- 0.05
  cols <- densCols(deseq_res$log2FoldChange, -log10(deseq_res$padj))
  p.size <- abs(deseq_res$log2FoldChange)
  p.size <- 3*p.size/max(p.size)
  
  plot.obj <- ggplot(deseq_res, aes(x = deseq_res$log2FoldChange, y = -log10(deseq_res$padj))) + 
              geom_point(aes(color = -log10(deseq_res$padj)), size = p.size) + 
              scale_color_gradient(low = "black", high = "red") + 
              geom_vline(xintercept = c(-1,1), alpha = 0.3) + 
              geom_hline(yintercept = -log10(alpha),alpha = 0.3) + 
              xlab(TeX("$log_2(foldChange)$")) + 
              ylab(TeX("$-log_{10}(p_{adj})$")) +
              xlim(c(-5,5)) + 
              labs(color = TeX("-log_{10}($p_{adj})"), size = TeX("$log_2(foldChange)$")) 
  
  return(plot.obj)
  
}

