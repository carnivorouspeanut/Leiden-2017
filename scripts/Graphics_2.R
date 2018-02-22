args <- commandArgs(TRUE)
polarity_colour <- args[1]
barpl <- args[2]
boxpl <- args[3]
colvec <- args[4]
out1 <- args[5]
out2 <- args[6]

#reading the table which was created manually!!!!
pol_col <- read.csv(polarity_colour, header = TRUE, sep = '\t', stringsAsFactors = FALSE)


load(barpl)
num_of_prots <- forplotvector[pol_col$group]
pol_col <- cbind(pol_col, num_of_prots)

#making the vector of colors for taxons
colorvector <- pol_col$color
names(colorvector) <- pol_col$group
save(colorvector, file = colvec)



#group barplot from "Graphics.R" and selecting color for every group
cols <- c()
for(i in(1:nrow(pol_col))){
  pol <- as.character(pol_col[i, ]$polarity)
  if(identical(pol, "ds") == TRUE){
    cols <- append(cols, "lightcoral")
  }
  else if(identical(pol, "minus") == TRUE){
    cols <- append(cols, "turquoise")
  }
  else if(identical(pol, "plus") == TRUE){
    cols <- append(cols, "cadetblue")
  }
  else if(is.na(pol)){
    cols <- append(cols, "grey")
  }
}


#pdf file with barplot about number of proteins in dataset
pdf(out1)
par(mar = par()$mar + c(4, 2, 2, 2))
pol_col <- cbind(pol_col, cols)
pol_col <- pol_col[order(pol_col$polarity, -pol_col$num_of_prots), ]
barplot(pol_col$num_of_prots, col = as.character(pol_col$cols), 
        names.arg=as.character(pol_col$group), las = 2, cex.axis = 0.7, cex.names = 0.7)
legend("topleft", legend = c("dsRNA", "ssRNA-", "ssRNA+", "unclassified RNA"), 
       col = c(unique(as.character(pol_col$cols))), lwd = 5, cex = 0.7)
dev.off()


load(boxpl)
#setting colors for boxplot
colors_for_boxplot <- c()
for(name in(names(bplist))){
  m <- which(as.character(pol_col$group) == name)
  colors_for_boxplot <- append(colors_for_boxplot, as.character(pol_col$cols)[m])
}


#pdf file with boxplot about protein length in dataset
pdf(out2)
par(mar = par()$mar + c(4, 2, 2, 2))
boxplot(bplist, las = 2, col = colors_for_boxplot, cex.axis = 0.7)
legend("topright", legend = c("dsRNA", "ssRNA-", "ssRNA+", "unclassified RNA"), 
       col = c(unique(as.character(pol_col$cols))), lwd = 5, cex = 0.7)
dev.off()
