pre_design_formulas <- function(dfid) {
    #' Provide easy pre-written design formulas
    #' If negative value is passed, information about
    #' available design formulas is printed. Otherwise 
    #' string with chosen design formula will be passed
    
    formulas <- list()
    contrasts <- list()
  
    formulas <- append(formulas,list('~id:replicate+replicate + num.id '))
    contrasts <- append(contrasts,list(c('num.id',"P1","P2")))
    formulas <- append(formulas,list('~replicate:tumor + tumor'))
    contrasts <- append(contrasts,list(c("tumor","tumor","non")))
    formulas <- append(formulas,list('~id:pseudo_replicate + tumor + tumor:replicate + tumor:id + id'))
    contrasts <- append(contrasts,list(c("num.id","P1","P2")))
    
    
    if (dfid >= 1 & dfid <= length(formulas)) {
        out <- list(design = formulas[[dfid]], contrast = contrasts[[dfid]])

    }  else {
        help_text <- paste(c(" 1.", formulas[[1]],
                              "\n Compare between different samples accounting",
                              "for replicate differences. Assumes that only one type",
                              "of spots are present (tumor or non-tumor\n",
                              "2.", formulas[[2]],
                              "\n Compare tumor vs. non-tumor within one",
                              "patient accounting for replicate differences\n",
                              "3.", formulas[[3]],
                              "\n Compare between patients and control for replicate",
                             "as well as tumor\n"
                        ), collapse = " ")


        cat(help_text)
        out <- NULL
        return(out)
    }
}
