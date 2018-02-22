args <- commandArgs(TRUE)
result <- args[1]
out <- args[2]
g1 <- args[3]
g2 <- args[4]

t1 <- read.csv (result, header=TRUE, sep='\t', stringsAsFactors=FALSE)
colnames(t1) <- c('GenomeSegmentID', 'GenomeSegmentLength', 'ProteinID', 'ProteinName', 'ProteinLength', 'Locations', 'Taxonomy', 'group', 'genome', 'name', 'len')
newtab <- t1[order(t1$ProteinLength, decreasing=TRUE), ]

#taxonomy parsing and adding "unclassified" inf if needed
df = data.frame(matrix(vector(), 0, 6, dimnames=list(c(), c("ProteinID", "Order", "Family", "Subfamily", "Genus", "The_major_taxon_range"))), stringsAsFactors=F)

taxonomy = newtab[ ,7]
for(i in (1:nrow(newtab))){
  a <- strsplit(taxonomy[i], ", ")[[1]]
  for(j in (1:length(a))){
    if(grepl("virales", a[j]) == TRUE){
      df[i, 2] = a[j]
    }
    if(grepl("viridae", a[j]) == TRUE){
      df[i, 3] = a[j]
    }
    if(grepl("virinae", a[j]) == TRUE){
      df[i, 4] = a[j]
    }
    if(grepl("[:alpha:]+virus", a[j]) == TRUE & j!=1){
      df[i, 5] = a[j]
    }
  }
  df[i, 1] = newtab[i, 3]
}


for(l in (1:nrow(newtab))){
  for(k in (2:5)){
    if(is.na(df[l, k]) == FALSE){
      df[l, 6] = df[l, k]
      l = l + 1
    }
    else {
      if(is.na(df[l, 5]) == TRUE & identical(newtab[l, 8], "ssRNA") == TRUE){
      df[l, 6] = "ssUnclassified"
      }
      if(is.na(df[l, 5]) == TRUE & identical(newtab[l, 8], "dsRNA") == TRUE){
        df[l, 6] = "dsUnclassified"
      }
      if(is.na(df[l, 5]) == TRUE & identical(newtab[l, 8], "unclassifRNA") == TRUE){
        df[l, 6] = "Unclassified"
      }
    }
      
  }
}

#adding parsed taxonomy to the initial table
final <- merge(newtab, df, by='ProteinID', sorted=FALSE)
write.table(final, file=out, sep = "\t")



#barplot vector
forplot <- table(final$The_major_taxon_range)
forplot <- as.data.frame(forplot)
colnames(forplot) <- c("taxons", "number_of_prots")
forplotvector <- forplot[, 2]
names(forplotvector) <- forplot[, 1]
save(forplotvector, file = g1)

bp <- subset(final, select=c("The_major_taxon_range", "ProteinLength"))

#making vectors for boxplot
maximums <- c()
forboxplot=list()
listofnames <- as.character(forplot[, 1])
for(i in (1:nrow(forplot))){
  vec <- c()
  for(j in(1:nrow(bp))){
    if(identical(bp[j,1], listofnames[i]) == TRUE){
      vec <- append(vec, bp[j,2])
    }

  }
  vec <- as.numeric(vec)
  maximums <- append(maximums, max(vec))
  forboxplot[[length(forboxplot)+1]]<-vec
}


#sorting vectors in list "forboxplot" by maximum - to get a beautiful boxplot, and naming it
names(forboxplot) <- listofnames
names(maximums) <- listofnames
maximums <- sort(maximums, decreasing = TRUE)


bplist=list()
for(j in(1:length(maximums))){
  for(k in(1:length(forboxplot))){
    if(identical(names(maximums[j]), names(forboxplot[k])) == TRUE){
      bplist[[length(bplist)+1]]<-forboxplot[[names(maximums[j])]]
    }
  }
}

names(bplist) <- names(maximums)

#boxplot list
save(bplist, file = g2)
