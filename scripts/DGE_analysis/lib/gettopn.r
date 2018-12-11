pth <- "/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/res/DGE_res1.csv"
n <- 100

save_top_n <- function()
df <- read.csv(pth)
topn <- order(df[['padj']])[1:n]
topframe <- df[c('X', 'padj')][topn,]
