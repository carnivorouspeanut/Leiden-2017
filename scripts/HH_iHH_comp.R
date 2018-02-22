directories <- list.dirs(path = "../Results/out", recursive = FALSE)
TAB <- data.frame(NULL, stringsAsFactors=FALSE)
prot_id_segment_id <- c()

for(each in directories){
  #do the following things only for folders with iterative HHsearch
  if(grepl("results_2.0.16_iterations", each)){
    each <- strsplit(each, "/")[[1]]
    each <- each[4]

    #parse the name of file to get IDs and foldernames
    filename = strsplit(each, "_")[[1]]
    GenomeSegmentID <- paste(filename[4], filename[5], sep = "_")
    HH_foldername <- paste("results_2.0.16_whole_pp", filename[4], filename[5], filename[6], filename[7], sep = "_")
    prot_segment <- paste(filename[4], filename[5], filename[6], filename[7], sep = "_")
    iter_name <- paste(filename[4], filename[5], filename[6], filename[7], "annotation_table.tsv", sep = "_")
    ipath <- paste("../Results/out", each, iter_name, sep = "/")
    
    #extract table with results of iHHsearch  
    iHHsearch_table = read.table(file = ipath, header = TRUE, stringsAsFactors = FALSE, dec = ".")
    tab <- iHHsearch_table[, c("index", "iteration", "Hit", "cl_from", "cl_to", "Prob")]    # after we obtain the data with the new version of the package, "q_id" should also be extracted
    names(tab)[names(tab)=="Prob"] <- "Prob_ihh"
    tab <- subset(tab, iteration > 1)
    

    #do the following only if we have more than 1 iteration (with every iteration but first)
    if(nrow(tab) >= 1){
      tab$Prob_hh <- NA
      tab$Eval_hh  <- NA
      tab$Len_hh  <- NA


      HHpath <- paste("../Results/out", HH_foldername, "iteration_1", sep = "/")
      #HHpath <- paste(HHpath, "/", sep = "")
      #HHsearch <- data.frame(NULL, stringsAsFactors = FALSE)
      cumulative_hits <-  list.files(path = HHpath, pattern = "cumulative_hits.tsv$")
      cumulative_hits <- paste(HHpath, cumulative_hits, sep = "/")
      HHsearch_table <- read.table(cumulative_hits, header = TRUE, stringsAsFactors = FALSE)

      for (i in 1:nrow(tab)) {
        idx <- which( HHsearch_table$Hit==tab$Hit[i] )
        # check that the hit of interest is mapped on the same region of polyprotein, by both iterative and regular HHsearch
        map_ihh <- tab$cl_from[i]:tab$cl_to[i]
        idx <- idx[ sapply(idx, function (j) any(HHsearch_table$h_from[j]:HHsearch_table$h_to[j] %in% map_ihh)) ]
        prot_id_segment_id <- c(prot_id_segment_id, prot_segment)
        if (length(idx) > 0) { tab[i, c("Prob_hh", "Eval_hh", "Len_hh")] <- HHsearch_table[idx[1], c("Prob", "E.value", "h_len")] }
      }
                              

    TAB <- rbind(TAB, tab)        


      
    }

  }
}


#create a plot
pdf("P_comp.pdf")
idx <- which(!is.na(TAB$Prob_hh))
maximum <- max(TAB$Prob_hh[idx], TAB$Prob_ihh[idx])
minimum <- min(TAB$Prob_hh[idx], TAB$Prob_ihh[idx])
plot(x=TAB$Prob_hh[idx], y=TAB$Prob_ihh[idx], pch=20, col=ifelse(TAB$Len_hh[idx] <= 50 | TAB$Eval_hh[idx] >= 10, "red", "green"), xlim = c(minimum, maximum), ylim = c(minimum, maximum), 
	xlab = "HHsearch", ylab = "iHHsearch", main = "Probability of hit")
if (length(idx) < nrow(TAB)) {
      points(x=rep(0, nrow(TAB)-length(idx)), y=TAB$Prob_ihh[-idx], pch=20, col="black")
                             }
legend('bottomleft', legend=c(expression('Len > 50 aa and E-val < 10, both hits'), expression('Len <= 50 aa or E-val >= 10, HHsearch hit'), 'no HHsearch hit'), pch=20, col=c('green','red','black'))

dev.off()

#white table with iterations of interest
TAB <- cbind(TAB, prot_id_segment_id)
write.table(TAB, file = "hits_correspondence_iHHsearch_HHsearch.txt", sep = '\t', quote = FALSE, row.names = FALSE)
