
#vectors with P for the first hit and color for scatterplot
P_iter <- c()
P_whole <- c()
colors_for_plot <- c()
hits <- c()
#load a table with color for each taxon and turn it into a named vector
taxon_colors <- read.table(file = "../Tables/groups_of_RNA_viruses_in_dataset_POLARITY_COLOR.txt", header = TRUE, stringsAsFactors = FALSE)
tax_cols <- taxon_colors$color
names(tax_cols) <- taxon_colors$group

#load a table with taxon inf about each entry and turn it into a named vector
taxon_inf <- read.table(file = "../Tables/Merged_table_of_RNAviruses.txt", header = TRUE, stringsAsFactors = FALSE)
tax_inf <- taxon_inf$Family
#fix all NAs in tax_inf
for(i in(1:nrow(taxon_inf))){
  if(is.na(tax_inf[i])){
    if(identical(taxon_inf[i, 8], "ssRNA")){
      tax_inf[i] = "ssUnclassified"
    }
    if(identical(taxon_inf[i, 8], "dsRNA")){
      tax_inf[i] = "dsUnclassified"
    }
    else{
      tax_inf[i] = "Unclassified"
    }
  }
}
names(tax_inf) <- taxon_inf$GenomeSegmentID


directories <- list.dirs(path = "../Results/out", recursive = FALSE)
TAB <- data.frame(NULL, stringsAsFactors=FALSE)

for(each in directories){
  if(grepl("results_2.0.16_iterations", each)){
    each <- strsplit(each, "/")[[1]]
    each <- each[4]
    #parse the name of file to get IDs and foldernames
    filename = strsplit(each, "_")[[1]]
    GenomeSegmentID <- paste(filename[4], filename[5], sep = "_")
        #add color of taxon to the vector
    #taxon <- tax_inf[GenomeSegmentID]
    #color <- tax_cols[taxon]
    #colors_for_plot <- c(colors_for_plot, color)
    HH_foldername <- paste("results_2.0.16_whole_pp", filename[4], filename[5], filename[6], filename[7], sep = "_")
  
    #extract table with results of iHHsearch
    iter_name <- paste(filename[4], filename[5], filename[6], filename[7], "annotation_table.tsv", sep = "_")
    ipath <- paste("../Results/out", each, iter_name, sep = "/")
    iHHsearch_table = read.table(file = ipath, header = TRUE, stringsAsFactors = FALSE, dec = ".")
    iHHsearch_table <- subset(iHHsearch, iteration > 1)
    tab <- iHHsearch_table[, c("Hit", "cl_from", "cl_to", "Prob")]    # after we obtain the data with the new version of the package, "q_id" should also be extracted
    names(tab)[names(tab)=="Prob"] <- "Prob_ihh"

    tab$Prob_hh <- NA
    tab$Eval_hh  <- NA
    tab$Len_hh  <- NA
    print(tab)


    
    if(nrow(tab) >= 1){
      for(each in (1:nrow(tab))){
        iteration_number <- tab[each, ]$iteration
        folder <- paste("iteration", iteration_number, sep = "_")
        HHpath <- paste("../Results/out", HH_foldername, folder, sep = "/")
        HHsearch <- data.frame(NULL, stringsAsFactors = FALSE)
        cumulative_hits <-  list.files(path = HHpath, pattern = "cumulative_hits.tsv$")
        for(table in cumulative_hits){
          table <- paste(HHpath, table, sep = "/")
          cum_table <- read.table(table, header = TRUE, stringsAsFactors = FALSE)
          HHsearch <- rbind(HHsearch, cum_table)
        }
    }
    for (i in 1:nrow(tab)) {
      idx <- which( HHsearch_table$Hit==tab$Hit[i] )
      # check that the hit of interest is mapped on the same region of polyprotein, by both iterative and regular HHsearch
      map_ihh <- tab$cl_from[i]:tab$cl_to[i]
      idx <- idx[ sapply(idx, function (j) any(HHsearch_table$h_from[j]:HHsearch_table$h_to[j] %in% map_ihh)) ]
      if (length(idx) > 0) { tab[i, c("Prob_hh", "Eval_hh", "Len_hh")] <- HHsearch_table[idx[1], c("Prob", "E.value", "h_len")] }
    }
                              

    TAB <- rbind(TAB, tab)
}
}
 
print(TAB)


#idx <- which(!is.na(TAB$Prob_hh))

 

#plot(x=TAB$Prob_hh[idx], y=TAB$Prob_ihh[idx], pch=20, col=ifelse(TAB$Len_hh[idx] <= 50 | TAB$Eval_hh[idx] >= 10, "red", "green"))

 

#if (length(idx) < nrow(TAB)) {

 

      #points(x=rep(0, nrow(TAB)-length(idx)), y=TAB$Prob_ihh[-idx], pch=20, col="black")

 

                             #}

#print(colors_for_plot)
#count limits to get a good-looking plot
#maximum <- max(P_iter, P_whole)
#minimum <- min(P_iter, P_whole) 
#pdf("P_trial.pdf")
#par(xpd = T, mar = par()$mar + c(4, 4, 4, 4))
#plot(P_whole, P_iter, xlim = c(minimum, maximum), ylim = c(minimum, maximum), pch = 16, xlab = "HHsearch", ylab = "iHHsearch", main = "Hit probability")
#lines(x=par()$usr[1:2], y=par()$usr[3:4], col = "red")
#legend(x=100, y=100, legend=taxon_colors$group, col=taxon_colors$color, lwd=5, cex = 0.69)
#dev.off()

#print(hits)
#save(hits, file = "best_hits.Rdata")









