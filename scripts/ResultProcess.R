args <- commandArgs(TRUE)
t <- args[1]
v <- args[2]
f_whole <- args[3]
f_iter <- args[4]
pplen <- args[5]
threshold <- args[6]
o1 <- args[7]
o1.1 <- args[8]
o2 <- args[9]
o3 <- args[10]
o4 <- args[11]
dir <- args[12]
database <- args[13]
per1 <- args[14]
per2 <- args[15]
dom1 <- args[16]
dom2 <- args[17]

pplen <- as.numeric(pplen)
threshold <- as.numeric(threshold)

#a function which searches for "two last" iterations
last_iter <- function(f){
  if(!file.exists(f)){
    return(FALSE)
  }
  tab <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
  for(i in(1:nrow(tab))){
    len = (tab[i, 2] - tab[i, 1] + 1)
    if(len > threshold){
      return(TRUE)
    }
  }
  return(FALSE)
}


#a function which returns a list: number of overlapping domans, summ. lehgth of all overlaps in one annotation table
overlap <- function(tab){
  minus = 0
  summoverlaps = 0
  un <- unique(tab$Hit)
  if(length(un) < nrow(tab)){
    for(k in(1:(nrow(tab)-1))){
      if(tab[k, ]$cl_to >= tab[k+1, ]$cl_from){
        minus = minus + 1
        summoverlaps = summoverlaps + (tab[k, ]$cl_to - tab[k+1, ]$cl_from)
      }
    }
    
  }
  result <- c(minus, summoverlaps)
  return(result)
}



#getting only long proteins from merged table
colourtable <- read.table(t, stringsAsFactors = FALSE, header = TRUE)
idx <- which(colourtable$ProteinLength >= pplen)
colourtable_5000 <- colourtable[idx, ]

#creating lists for plots
domains_whole <- c()
domains_iter <- c()
colors <- c()
unique1 <- c()
unique2 <- c()
percentage1 <- c()
percentage2 <- c()
domains_iter_notBF <- c()
domains_iter_BF <- c()

names_BF <- c()
names_notBF <- c()


#vector of colors for each major taxon created by Set_colors_for_taxons.R
load(v)

#loading a PFAMdb
PFAM <- read.csv(database, stringsAsFactors = FALSE, header = FALSE)
PFAM <- na.omit(PFAM)
colnames(PFAM) <- c("name")


fl1 <- unique(colourtable_5000$The_major_taxon_range)
fl2 <- c()
for(i in(1:nrow(colourtable_5000))){ 
  #adding a color for taxon
  add = colorvector[colourtable_5000[i, 16]]
  colors <- append(colors, add)
  #creating a little colorvector for legend
  fl2 <- append(fl2, colorvector[fl1[i]])
  
  
  
  #path
  segID <- colourtable_5000[i, 2]
  prID <- colourtable_5000[i, 1]
  ID <- paste(c(segID, prID), collapse = "_")
  path_folder_whole <- path_whole <- paste(c(f_whole, ID), collapse = "/")
  path_folder_iter <- path_whole <- paste(c(f_iter, ID), collapse = "/")
  path_whole <- paste(c(path_folder_whole, "/", ID, "_annotation_table.tsv"), collapse = "")
  path_iter <- paste(c(path_folder_iter, "/", ID, "_annotation_table.tsv"), collapse = "")
  
  #check if there is a table with domains; if so...
  if(!file.exists(path_whole)){
    domains_whole <- append(domains_whole, 0)
    percentage1 <- append(percentage1, 0)
  }


  else if(file.exists(path_whole)){
    tab1 <- read.table(path_whole, header = TRUE, stringsAsFactors = FALSE)
	#add column with domain names
    Domain_name <- c()
    for(j in (1:nrow(tab1))){
      PFAMid <- tab1[j, ]$Hit
      DBsubset <- PFAM[grep(PFAMid, PFAM$name), ]
      if(identical(length(DBsubset), 0) == FALSE){
        DBsubset <- strsplit(DBsubset, ";")
        DBsubset <- DBsubset[[1]]
        Domain_name <- append(Domain_name, tail(DBsubset, n=1))
      }

    }
    #tab1 <- cbind(tab1, Domain_name)
    #write.table(tab1, file = path_whole, sep = '\t')

    domains_whole <- append(domains_whole, nrow(tab1))
    unique1 <- append(unique1, tab1$Hit)
    percentage1 <- append(percentage1, (sum(tab1$cl_len))/colourtable_5000[i, 5])
  }


  if(!file.exists(path_iter)){
    domains_iter <- append(domains_iter, 0)
    percentage2 <- append(percentage2, 0)
    domains_iter_notBF <- append(domains_iter_notBF, 0)
    domains_iter_BF <- append(domains_iter_BF, 0)

    names_notBF <- append(names_notBF, prID)
    names_BF <- append(names_BF, prID)
    
  }
  
  
  else if(file.exists(path_iter)){
    tab2 <- read.table(path_iter, header = TRUE, stringsAsFactors = FALSE)
	#add column with domain names
    Domain_name <- c()
    for(j in (1:nrow(tab2))){
      PFAMid <- tab2[j, ]$Hit
      DBsubset <- PFAM[grep(PFAMid, PFAM$name), ]
      if(identical(length(DBsubset), 0) == FALSE){
        DBsubset <- strsplit(DBsubset, ";")
        DBsubset <- DBsubset[[1]]
        Domain_name <- append(Domain_name, tail(DBsubset, n=1))
      }

    }
    #tab2 <- cbind(tab2, Domain_name)
    #write.table(tab2, file = path_iter, sep = '\t')


    unique2 <- append(unique2, tab2$Hit)
    ol <- overlap(tab2)
    minus <- ol[1]
    summoverlaps <- ol[2]
    domains_iter <- append(domains_iter, (nrow(tab2)-minus))
    percentage2 <- append(percentage2, (sum(tab2$cl_len) - summoverlaps)/colourtable_5000[i, 5])
    
    #check if we have BF iterations; if so, split the table ("normal" search/BF search)
    withouthits <- paste(c(path_folder_iter, "regions_without_hits.tsv"), collapse = "/")
    if(isTRUE(last_iter(withouthits))){
      files_it <- list.files(path_folder_iter, pattern = "query_regions_coo")
      iterations <- c()
	#here we extract number of iterations from filenames
      for(m in(1:length(files_it))){
        gr <- gsub('.*_([0-9]+).*', '\\1', files_it[m])
        iterations <- append(iterations, gr)
      }
      second_iteration <- sort(iterations, partial=length(iterations)-1)[length(iterations)-1]
      tab_notBF <- subset(tab2, iteration < second_iteration)
      tab_BF <- subset(tab2, iteration >= second_iteration)
      
      ol_notBF <- overlap(tab_notBF)
      minus_notBF <- ol_notBF[1]
      summoverlaps_notBF <- ol_notBF[2]
      domains_iter_notBF <- append(domains_iter_notBF, (nrow(tab_notBF)-minus_notBF))
      names_notBF <- append(names_notBF, prID)
      
      ol_BF <- overlap(tab_BF)
      minus_BF <- ol_BF[1]
      summoverlaps_BF <- ol_BF[2]
      domains_iter_BF <- append(domains_iter_BF, (nrow(tab_BF)-minus_BF))
      names_BF <- append(names_BF, prID)
      
      
    }
    
  } 


}
names(domains_iter_BF) <- names_BF

#setting "heatmap" colors for the following scatterplot
forheatmap <- rep(0, length(domains_iter_BF))
for (i in(1:length(forheatmap))){
  vec <- c(domains_iter_BF[i], domains_iter_notBF[i])
  for(j in(1:length(forheatmap))){
    cmp <- c(domains_iter_BF[j], domains_iter_notBF[j])
    if(all(vec == cmp)){
      forheatmap[i] = forheatmap[i] + 1
    }
  }
}

heatmap <- colorRampPalette(c("blue", "red"))
heatmap <- heatmap(length(unique(forheatmap)))
sorted_forheatmap <- sort(forheatmap, decreasing = FALSE)
names(heatmap) <- unique(sorted_forheatmap)
heatcolors <- heatmap[as.character(forheatmap)]



#scatterplot domains biologically meaningful/BF
pdf(o1)
par(xpd = T, mar = par()$mar + c(4, 4, 4, 4))
maxx <- max(domains_iter_BF)
maxy <- max(domains_iter_notBF)
plot(domains_iter_BF, domains_iter_notBF, col = heatcolors, pch=16,cex=1.3, xlim = c(-1, maxx + 1), ylim = c(-1, maxy + 1), xaxs = 'i', yaxs = 'i', xlab="Number of domains found by brute force", ylab="Number of domains found by biologically \nmeaningful iterations iterations")
names(domains_iter_BF) <- names_BF
names(domains_iter_notBF) <- names_notBF
#text(x = domains_iter_BF, y = domains_iter_notBF, labels = names_BF, cex = 0.5, col = "black", pos = 3)
legend(x=maxx + 1, y=maxy + 1, legend=unique(sorted_forheatmap), col=heatmap, lwd=5, cex = 0.69)
dev.off()

#scatterplot percentage of domains biologically meaningful/BF
biol_mean <- domains_iter_notBF/(domains_iter_BF + domains_iter_notBF)
brute_force <- domains_iter_BF/(domains_iter_BF + domains_iter_notBF)


pdf(o1.1)
par(xpd = T, mar = par()$mar + c(4, 4, 4, 4))
biol_mean <- biol_mean[!is.na(biol_mean)]
brute_force <- brute_force[!is.na(brute_force)]
limit <- pmax(biol_mean, brute_force)
maximum <- max(limit)
plot(brute_force, biol_mean, col = colors, pch=16,cex=1.3, xlim = c(-0.2, maximum + 0.2), ylim = c(-0.2, maximum + 0.2), xaxs = 'i', yaxs = 'i', xlab="Percentage of domains found by brute force", ylab="Percentage of domains found by \nbiologically meaningful iterations")
lines(x=par()$usr[1:2], y=par()$usr[3:4], col = "red")
legend(x=maximum + 0.2, y=maximum + 0.2, legend=fl1, col=fl2, lwd=5, cex = 0.69)
dev.off()

#scatterplot regular iterative
for_num_of_domains <- (as.numeric(domains_iter) - as.numeric(domains_whole))/as.numeric(domains_iter)
pdf(o2)
par(xpd = T, mar = par()$mar + c(4, 4, 4, 4))
maximum = max(domains_whole)
plot(domains_whole, for_num_of_domains, col = colors, pch=16, xlim = c(-1, maximum + 1), ylim = c(-0.1, 1.1), xaxs = 'i', yaxs = 'i', xlab="HHSearch", ylab="iHHsearch-HHsearch/iHHsearch", main = "Number of domains")
#x <- seq(-1, (maximum+1))
#y <- rep(1, length(x))
#lines(x, y, col = "red")
legend(x=maximum + 1.01, y=1.1, legend=fl1, col=fl2, lwd=5, cex = 0.69)
dev.off()
#write both "number of domains" vectors in Rdata to compare them with "genbank"-annotated rates (script Percentage_comparison.R)
names(domains_whole) <- colourtable_5000$ProteinID
names(domains_iter) <- colourtable_5000$ProteinID
save(domains_whole, file=dom1)
save(domains_iter, file=dom2)



#% annotated
for_per_of_domains <- (as.numeric(percentage2) - as.numeric(percentage1))/as.numeric(percentage2)
pdf(o3)
par(xpd = T, mar = par()$mar + c(4, 4, 4, 4))
plot(percentage1, for_per_of_domains, col = colors, pch=16,  xlim = c(-0.03, 1), ylim = c(-0.2, 1.2), xaxs = 'i', yaxs = 'i', xlab="HHSearch", ylab="iHHsearch-HHsearch/iHHsearch", main = "Polyprotein coverage")
#x <- seq(-0.03, 1.03, by=0.001)
#y <- rep(1, length(x))
#lines(x, y, col = "red")
legend(x=1, y=1.2, legend=fl1, col=fl2, lwd=5, cex = 0.69)
dev.off()
#write both "percentage" vectors in Rdata to compare them with "genbank"-annotated rates (script Percentage_comparison.R)
names(percentage1) <- colourtable_5000$ProteinID
names(percentage2) <- colourtable_5000$ProteinID
save(percentage1, file=per1)
save(percentage2, file=per2)


#barplot with number of unique PFAM identifiers
pdf(o4)
par(xpd = T, mar = par()$mar + c(4, 4, 4, 4))
bp1 <- unique(unique1)
bp2 <- unique(unique2)
forbp <- c(length(bp1), length(bp2))
names(forbp) <- c("regular HHSearch", "iterative HHSearch")
barplot(forbp, col=c("turquoise3", "tomato1"), ylab="Number of unique PFAM identifiers")
dev.off()



#tables with domains not found by iterative search (in spetial folder in Results)
for(i in(1:nrow(colourtable_5000))){
  if(domains_whole[i] > (domains_iter[i])){
    if(!file.exists(dir)){
       dir.create(dir)
       filename <- paste(c(dir, colourtable_5000[i, 2], "_", colourtable_5000[i, 1], ".txt"), collapse = "")
       t1name <- paste(c(f, "results_2.0.16_whole_pp_", colourtable_5000[i, 2], "_", colourtable_5000[i, 1], "/", colourtable_5000[i, 2], "_", colourtable_5000[i, 1], "_annotation_table.tsv"), collapse = "")
       t2name <- paste(c(f, "results_2.0.16_iterations_", colourtable_5000[i, 2], "_", colourtable_5000[i, 1], "/", colourtable_5000[i, 2], "_", colourtable_5000[i, 1], "_annotation_table.tsv"), collapse = "")
       t1 <- read.table(t1name, header=TRUE)
       t2 <- read.table(t2name, header=TRUE)
       diff <- setdiff(t1$Hit, t2$Hit)
       diff <- as.vector(diff)
       t3 <- t1[which(t1$Hit %in% diff), ] 
       write.table(t3, file = filename, row.names = FALSE)
    }
    else{
       filename <- paste(c(dir, colourtable_5000[i, 2], "_", colourtable_5000[i, 1], ".txt"), collapse = "")
       t1name <- paste(c(f, "results_2.0.16_whole_pp_", colourtable_5000[i, 2], "_", colourtable_5000[i, 1], "/", colourtable_5000[i, 2], "_", colourtable_5000[i, 1], "_annotation_table.tsv"), collapse = "")
       t2name <- paste(c(f, "results_2.0.16_iterations_", colourtable_5000[i, 2], "_", colourtable_5000[i, 1], "/", colourtable_5000[i, 2], "_", colourtable_5000[i, 1], "_annotation_table.tsv"), collapse = "")
       t1 <- read.table(t1name, header=TRUE)
       t2 <- read.table(t2name, header=TRUE)
       diff <- setdiff(t1$Hit, t2$Hit)
       diff <- as.vector(diff)
       t3 <- t1[which(t1$Hit %in% diff), ] 
       write.table(t3, file = filename, row.names = FALSE, sep = '\t')
    }
  }
}
