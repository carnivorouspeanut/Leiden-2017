args <- commandArgs(TRUE)
input <- args[1]
output <- args[2]

colourtable <- read.table(input, stringsAsFactors = FALSE, header = TRUE)
namesforvector <- unique(colourtable$The_major_taxon_range)

colorvector <- c()

for(i in(1:length(namesforvector))){
  index = (i + 1)*12
  colorvector[i] <- colors()[index]
}

names(colorvector) <- namesforvector


save(colorvector, file=output)

