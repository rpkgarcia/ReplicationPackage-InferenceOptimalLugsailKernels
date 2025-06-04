library(xtable)

fitted_cv <- "fitted_CV.csv"
the_fits <- read.csv(fitted_cv)
the_fits <- the_fits[the_fits$d==5, ]

the_kernel <- "Bartlett"
latex_table <- the_fits[the_fits$kernel == the_kernel,]

# Change order, remove some columns 
latex_table <- latex_table[, c(9,1:5)]
# 
# for(i in 2:6){
#   latex_table[, i] <- as.character(formatC(as.numeric(latex_table[,i]), 3))
# }
# 
# latex_table$alpha <- as.character(latex_table$alpha)


latex_table <- xtable(latex_table, digits = 3)


colnames(latex_table) <- c("$\\\alpha$", "$a_0$", "$a_1$", 
                           "$a_2$", "$a_3$", "$R^2$")
#addtorow <- list()
#addtorow$pos <- list(-1)
#addtorow$command <- c(" \\hline &  & \\multicolumn{5}{c|}{\textit{Mother}}& \\multicolumn{2}{c|}{\textit{Zero}}  \\\\\n")

#rownames(latex_table) <- c("Bartlett", "", "  ", "   ")

rownames(latex_table) <- c("Mother", "", " ", "  ", 
                           "Zero", "   ", "    ", "     ", 
                           "Adpative", "      ", "       ", "        ", 
                           "Over", "         ", "          ", "           ")

italic <- function(x){
  paste0('{\\emph{ ', x, '}}')
}

align(latex_table) <- "cc||ccccc"

print(latex_table, sanitize.rownames.function = italic,
      #add.to.row = addtorow,
      hline.after = c(-1,  4, 8, 12), 
      include.colnames = T)
