#!/usr/bin/Rscript


sh <- suppressPackageStartupMessages
sh(library(enrichR))
sh(library(org.Hs.eg.db))
sh(library(pathfindR))
sh(library(ReactomePA))
sh(library(optparse))
sh(library(AnnotationDbi))
sh(library(gridExtra))


parser <- OptionParser()

parser <- add_option(parser,
                     c("-d","--dge_dir"),
                     default = NULL,
                     help =""
                     )
parser <- add_option(parser,
                     c("-f","--file"),
                     default = NULL,
                     help = ""
                     )

parser <- add_option(parser,
                     c("-i","--include"),
                     default = NULL,
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

parser <- add_option(parser,
                     c("-r","--reactomePA"),
                     default = FALSE,
                     action ="store_true")


args <- parse_args(parser)


if (!(is.null(args$include))) {
  include_pth <- args$include
  include_names <- read.delim(include_pth, sep ="\n",header = F, stringsAsFactors = F)
  files <- as.character(include_names$V1)
} else if (!is.null((args$file))){
  files <- args$file
} else {
    quit(status = 1)  
}
  

for (file in files) {
  bname <- gsub("^\\.\\/","",basename(file))
  odir <- paste(c(dirname(normalizePath(file)),gsub("\\.fancy\\.tsv","\\_EA\\_",bname)),collapse = '/')
  print(odir)
  #oname_base <- gsub("fancy\\.tsv","",basename(file))
  if (!(dir.exists(odir))) {
    dir.create(odir)
  }

  dgeres <- read.csv(file, header = 1, row.names = 1)
  
  if (!("symbol" %in% colnames(dgeres))) {
    
    dgeres$symbol <- mapIds(org.Hs.eg.db,
                            keys=as.character(dgeres$genes),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
  }
  
  dgeres <- dgeres[!is.na(dgeres$symbol),]
  
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
    pdf(file = paste(c(odir, "enrichR.res.pdf"),collapse = "/"), width = 40, height = 40)
    for (db in names(enriched)) {
      adj.enriched[[db]] <- enriched[[db]][which(enriched[[db]]$Adjusted.P.value < 0.01),]
      oname_spec <- paste(c("enrichR.",db,".csv"),collapse = "")
      write.csv(adj.enriched[[db]], file = paste(c(odir,oname_spec),collapse="/"), row.names = T, quote = F)
    }
    dev.off()
  }
  
  if (args$pathfindR) {
    try(pf_res <- run_pathfindR(dgeres[c("symbol","logFC","FDR")], visualize_pathways  = F,
                            output = paste(c(odir,'pathfindR'),collapse="/"),
                            p_val_threshold = 0.05))
  }

}

if (args$reactomePA) {
  entrez_genes <- mapIds(org.Hs.eg.db, keys=rownames(dgeres),column="ENTREZID", keytype="ENSEMBL", multiVals="first")
  entrez_genes <- entrez_genes[!(is.na(entrez_genes))]
  repares <- try( x <-enrichPathway(gene=entrez_genes,pvalueCutoff=0.05, readable=T))
  if(!(class(repares) == "try-error")) {
       pdf(file = paste(c(odir,"ReactomePA_plots.pdf"),collapse = '/'), height = 20, width = 20)
       try(print(barplot(x)))
       try(print(emapplot(x)))
       try(print(dotplot(x)))
       try(print(cnetplot(x, categorySize = 'pvalue', foldChange = dgeres$logFC)))
       dev.off()
       pdf(file = paste(c(odir,"ReactomePA_table.pdf"),collapse = '/'), height = 30, width = 60)
       try(print(grid.table(as.data.frame(x))))
       dev.off()
       }
  }