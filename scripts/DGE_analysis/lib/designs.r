pre_design_formulas <- function(dfid) {
    #' Provide easy pre-written design formulas
    #' If negative value is passed, information about
    #' available design formulas is printed. Otherwise 
    #' string with chosen design formula will be passed
    
    f1 <- '~id:pseudo.replicate+pseudo.replicate + id '
    f2 <- '~pseudo.replicate:tumor + tumor'
    f3 <- '~id:psuedo_replicate + tumor + tumor:pseudo.replicate + tumor:id + id'
    formulas <- list(f1,f2,f3)
    
    if (dfid >= 1 & dfid <= length(formulas)) {
        out <- formulas[[dfid]]

    }  else {
        help_text <- paste(c(" 1.", f1,
                              "\n Compare between different samples accounting",
                              "for replicate differences. Assumes that only one type",
                              "of spots are present (tumor or non-tumor\n",
                              "2.", f2,
                              "\n Compare tumor vs. non-tumor within one",
                              "patient accounting for replicate differences",
                              "3.", f3,
                              "\n Compare between patients and control for replicate",
                             "as well as tumor\n"
                        ), collapse = " ")


        cat(help_text)
        out <- NULL
        return(out)
    }
}
