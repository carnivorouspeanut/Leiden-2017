args <- commandArgs(TRUE)
t1 <- args[1]
matp <- args[2]
cdss <- args[3]
pplen <- args[4]
pe1 <- args[5]
pe2 <- args[6]
co <- args[7]
o1 <- args[8]
o2 <- args[9]

pplen <- as.numeric(pplen)


#function which searches for nested proteins
overlaps <- function(t1, t2){
  minus <- c()
  for(j in(1:nrow(t1))){
    for(k in(1:nrow(t2))){
      s1 <- t1[j, ]$start
      s2 <- t2[k, ]$start
      e1 <- t1[j, ]$end
      e2 <- t2[k, ]$end
      if((s1 >= s2) && (e1 < e2) || ((s1 > s2) && (e1 <= e2))){
        minus <- append(minus, t1[j, ]$polyprotein)
      }
    }
  }
  return(unique(minus))
}


#function which searches the seqment length covered by all peptides in table
coverage <- function(tab, prot){
  num_of_peptides = 0
  result <- list()
  vstart <- prot$start
  vend <- prot$end
  vec <- logical(length = (vend - vstart + 1))
  if(nrow(tab) == 0){
    result[[1]] = 0
    result[[2]] = 0
    return(result)
  }
  else{
    for(i in (1:nrow(tab))){
      cur_start <- tab[i, ]$start
      cur_end <- tab[i, ]$end
      if(cur_start >= vstart & cur_end <= vend){
        vec_start <- cur_start - vstart
        vec_end <- vec_start + (cur_end - cur_start + 1)
        num_of_peptides = num_of_peptides + 1
        for(j in(vec_start:(vec_end - vec_start + 1))){
          vec[j] = TRUE
        }
      }
    }
    result[[1]] = sum(vec)
    result[[2]] = num_of_peptides
    return(result)
  }
}


#load colorvector
load(co)
#load vectors for plot with percentage of prots annotated by HHsearch(percentage1, percentage2)
load(pe1)
load(pe2)

#load vectors for plot with number of domains annotated by HHsearch(percentage1, percentage2)
#load("dom1.Rdata")
#load("dom2.Rdata")



#read table with taxonomy to set colors for points on scatterplots
merged_table <- read.csv(t1, header = TRUE, stringsAsFactors = FALSE, sep = '\t')

merged_table <- subset(merged_table, merged_table$ProteinLength >= pplen)
fl1 <- unique(merged_table$The_major_taxon_range)
fl2 <- colorvector[fl1]


#open tables with mat peptides and CDSs
mat <- read.csv(matp, header = FALSE, stringsAsFactors = FALSE, sep = '\t')
cds1 <- read.csv(cdss, header = FALSE, stringsAsFactors = FALSE, sep = '\t')
colnames(cds1) <- c("genome_segment", "polyprotein", "start", "end")
cds1 <- cds1[, c("polyprotein", "genome_segment", "start", "end")]


#subset > 1000 aa
cds <- subset(cds1, (cds1$end - cds1$start + 1 >= pplen*3))
colnames(mat) <- c("genome_segment", "polyprotein", "mat_product", "start", "end")
mat <- mat[, c("genome_segment", "polyprotein", "mat_product", "start", "end")]
mylist <- split(cds, cds$genome_segment)




#remove nested CDSs
minus <- c()
for(i in(1:length(mylist))){
  tab <- mylist[[i]]
  m <- overlaps(tab, tab)
    minus <- append(minus, m)
}
cds <- cds[!grepl(paste(minus, collapse = "|"), cds$polyprotein), ]
mylist <- split(cds, cds$genome_segment)



#count an annotated rate for every protein
perc <- c()
prots <- c()
for(i in(1:nrow(cds))){
  tab <- subset(mat, genome_segment == cds[i, ]$genome_segment)
  percentage <- coverage(tab, cds[i, ])
  perc <- append(perc, (percentage[[1]]/(cds[i, ]$end - cds[i, ]$start + 1)))
  prots <- append(prots, cds[i, ]$polyprotein)
}

names(perc) <- prots
print(names(perc))
#extracl all names present in perc vector from "HH-calculated" percentage vectors
#print(grepl(paste(names(percentage1), collapse = "|"), names(perc)))

p1 <- perc[grepl(paste(names(percentage1), collapse = "|"), names(perc))]
p2 <- perc[grepl(paste(names(percentage2), collapse = "|"), names(perc))]

#vice versa + setting colors
p_regular <- percentage1[grepl(paste(names(p1), collapse = "|"), names(percentage1))]
colors_regular <- c()
for(i in(1:length(p_regular))){
  merged_str <- merged_table[grepl(paste(names(p_regular[i])), merged_table$ProteinID), ]
  taxon <- merged_str$The_major_taxon_range
  colors_regular <- append(colors_regular, colorvector[taxon])
}


p_iterative <- percentage2[grepl(paste(names(p2), collapse = "|"), names(percentage2))]
colors_iterative <- c()
for(i in(1:length(p_iterative))){
  merged_str <- merged_table[grepl(paste(names(p_iterative[i])), merged_table$ProteinID), ]
  taxon <- merged_str$The_major_taxon_range
  colors_iterative <- append(colors_iterative, colorvector[taxon])
}

#plots
pdf(o1)
par(xpd=T, mar = par()$mar + c(5, 5, 5, 5))
plot(p1, p_regular, col = colors_regular, pch = 16, 
     xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1), xaxs = "i", yaxs = "i", 
     xlab="Percentage of polyprotein annotated in GB", ylab="Percentage of polyprotein annotated \n by regular HHsearch")
lines(x = par()$usr[1:2], y = par()$usr[3:4], col = "red")
legend(x = par()$usr[2], y = par()$usr[4], legend = fl1, col = fl2, lwd=5, cex = 0.69)
dev.off()


pdf(o2)
par(xpd=T, mar = par()$mar + c(5, 5, 5, 5))
plot(p2, p_iterative, col = colors_regular, pch = 16, 
     xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1), xaxs = "i", yaxs = "i", 
     xlab="Percentage of polyprotein annotated in GB", ylab="Percentage of polyprotein annotated \n by iterative HHsearch")
lines(x = par()$usr[1:2], y = par()$usr[3:4], col = "red")
legend(x = par()$usr[2], y = par()$usr[4], legend = fl1, col = fl2, lwd=5, cex = 0.69)
dev.off()
