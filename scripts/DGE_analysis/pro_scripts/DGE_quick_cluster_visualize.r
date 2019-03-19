#!/usr/bin/Rscriptquic

library(ggplot2)
library(optparse)

parser <- OptionParser()

parser <- add_option(parser, c("-f","--feature_file"))
args <- parse_args(parser)
df <- read.csv(args$feature_file, header = T, row.names = 1, sep ="\t")
df$tumor_id <- factor(df$tumor_id)

g <- ggplot(data = df, aes(x = xcoord, y = ycoord)) + 
  geom_point(aes(color = tumor_id), size = 5)


X11()
plot(g)
cat("Press [enter] to continue ...")
fil <- readLines(con="stdin", 1)
cat(fil, "\n")


