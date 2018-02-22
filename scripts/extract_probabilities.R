tab <- read.table("hits_correspondence_iHHsearch_HHsearch.txt", header = TRUE, stringsAsFactors = FALSE)

new_tab <- subset(tab, select = c("index", "prot_id_segment_id", "Hit", "Prob_ihh", "Prob_hh"))
#this is the vector of hits probabilities on subsequences
Prob_seq <- c()

for(i in (1:nrow(new_tab))){
  #find a file containing information about hits found by HHsearch on current FASTA subsequence
  foldername1 <- paste(new_tab[i, ]$prot_id_segment_id, new_tab[i, ]$Hit, new_tab[i, ]$index, sep = "_")
  foldername1 <- paste("../FASTA_best_hit/out", foldername1, "iteration_1", sep = "/")
  foldername2 <- list.files(path = foldername1, pattern = "cumulative_hits.tsv$")
  tablename <- paste(foldername1, foldername2, sep = "/")
  #check if we even have this file
  if(!file.exists(tablename)){
    print(tablename)
  }
  tab_PFAM_HH <- read.table(tablename, header = TRUE, stringsAsFactors = FALSE)
  #extract probability of a hit we are interested in
  Hit <- data.frame(tab_PFAM_HH[tab_PFAM_HH$Hit == new_tab[i, ]$Hit, ]$Prob)
  Hit <- Hit[1, ] # we take only the first hit, if there is more than one
  Prob_seq <- c(Prob_seq, Hit)
}

Prob_seq <- as.numeric(unlist(Prob_seq))
new_tab <- cbind(new_tab, Prob_seq)

write.csv(new_tab, file = "Probabilities_HH_iHH_partHH.txt")
#plot
pdf("iHH_HHpart.pdf")
lim_min = min(new_tab$Prob_seq, new_tab$Prob_ihh)
lim_max = max(new_tab$Prob_seq, new_tab$Prob_ihh)
plot(new_tab$Prob_ihh, new_tab$Prob_seq, pch = 16, xlab = "iHHsearch", ylab = "HHsearch_on_hit_sequence", xlim = c(lim_min, lim_max), ylim = c(lim_min, lim_max))
abline(0, 1, col = "red")
dev.off()
