library(ggplot2)

library(ggplot2)
library(latex2exp)

volcano <- function(deseq_res, title) {
  alpha <- 0.05
  cols <- densCols(df$log2FoldChange, -log10(df$padj))
  p.size <- abs(deseq_res$log2FoldChange)
  p.size <- 3*p.size/max(p.size)
  
  print(p.size)
  
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



pth <- '/home/alma/ST-2018/CNNp/DGE/res/DGEresults/DGEres6.190119/23508.tsv'
df <- read.csv(file = pth, sep = ',', header = TRUE)
df <- df[!(is.na(df$padj)),]
print(volcano(df,"test"))