#!/usr/bin/Rscript

library(enrichR)
library(org.Hs.eg.db)
library(pathfindR)
library(optparse)
library(AnnotationDbi)


parser <- OptionParser()

parser <- add_option(parser,
                     c("-d","--dge_dir"),
                     help =""
                     )


parser <- add_option(parser,
                     c("-i","--include"),
                     help =""
)

parser <- add_option(parser,
                     c("-e","--enrichR"),
                     default = FALSE,
                     action = "store_true")

parser <- add_option(parser,
                     c("-p","--pathfindR"),
                     default = FALSE,
                     action ="store_true")

args <- parse_args(parser)

dge_dir <- args$dge_dir
include_pth <- args$include
include_names <- read.delim(include_pth, sep ="\n",header = F, stringsAsFactors = F)
files <- as.character(sapply(include_names$V1,
                             function(x) paste(c(dge_dir,x,list.files(paste(c(dge_dir,x),collapse = "/"),pattern = "*fancy*")),collapse="/")))

print(files)

for (file in files) {
  
  odir <- dirname(file)
  oname_base <- gsub("fancy\\.tsv","",basename(file))
  print(file)
  
  dgeres <- read.csv(file, header = T, row.names = 1)
  
  if (!("symbol" %in% colnames(dgeres))) {
    
    dgeres$symbol <- mapIds(org.Hs.eg.db,
                            keys=as.character(dgeres$genes),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
  }
  
  dgeres <- dgeres[!is.na(dgeres$symbol),]
  
  print(dgeres)
  
  
  if (args$enrichR) {
    dbs <- c("GO_Molecular_Function_2017b",
             "GO_Cellular_Component_2017b",
             "GO_Biological_Process_2017b",
             "KEGG_2016",
             "Cancer_Cell_Line_Encyclopedia",
             "Panther_2016",
             "NCI-60_Cancer_Cell_Lines",
             "Disease_Signatures_from_GEO_up_2014")
    
    enriched <- enrichr(genes = as.character(dgeres$symbol),dbs)
    adj.enriched <- list()
  
    for (db in names(enriched)) {
      adj.enriched[[db]] <- enriched[[db]][which(enriched[[db]]$Adjusted.P.value < 0.01),]
      oname_spec <- paste(c(oname_base,"enrichR.",db,".csv"),collapse = "")
      write.csv(adj.enriched[[db]], file = paste(c(odir,oname_spec),collapse="/"), row.names = T, quote = F)
    }
  }
  
  if (args$pathfindR) {
    try(pf_res <- run_pathfindR(dgeres[c("symbol","logFC","FDR")],
                            output = paste(c(odir,"pathfindR"),collapse="/"),
                            p_val_threshold = 0.01))
  }

}