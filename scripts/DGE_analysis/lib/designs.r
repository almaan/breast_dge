#Desgin formulas for DESeq2


pre_design_formulas <- function(dfid) {


    f1 <- as.formula('~id:pseudo.replicate+pseudo.replicate+id')
    f2 <- as.formula('~tumor+pseduo.replicate:tumor+pseudo.replicate')

    #TODO:Need to turn this into fstring in order for nice representation
    if (dfid >= 1) {
        formulas <- list(f1,f2)
        out <- formulas[dfid]

    }  else {
        help_text <- paste(c("1. ~ id:pseudo.replicate + pseudo.replicate + id\n",
                        "\t used to compare between different samples accounting",
                        "for replicate differences",
                    "2. ~ tumor + pseduo.replicate:tumor pseudo.replicate\n",
                        "\t used to compare tumor vs. non-tumor within one",
                        "patient accounting for replicate differences"
                        ), collapse = " ")


        print(help_text)
        out <- NULL
        return(out)
    }
}
