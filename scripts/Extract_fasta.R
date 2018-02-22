TAB <- read.table("hits_correspondence_iHHsearch_HHsearch.txt", header = TRUE, stringsAsFactors = FALSE)

for(i in(1:nrow(TAB))){
  #path to HHsearch FASTA
  fasta_name <- paste(TAB[i, ]$prot_id_segment_id, "fasta", sep = ".")
  #number of index and hit
  index <- TAB[i, ]$index
  hit <- TAB[i, ]$Hit
  #coordinates of current domain
  start = TAB[i, ]$cl_from
  end = TAB[i, ]$cl_to
  fasta_path <- paste("/home/eocheredko/Desktop/parse_ncbi_table/fasta_files", fasta_name, sep = "/")
  fasta = read.table(file = fasta_path, header = FALSE, stringsAsFactors = FALSE)
  #merge all strings of FASTA file into one sequence
  fs <- ""
  for(j in(2:nrow(fasta))){
	fs <- paste(fs, fasta[j, ], sep = "")
      }
  #a name for FASTA with the current hit
  newfasta <- paste(TAB[i, ]$prot_id_segment_id, hit, index, sep = "_") 
  newfasta <- paste(newfasta, "fasta", sep = ".")
  newfastaname <- paste("../FASTA_best_hit_extended", newfasta, sep = "/")
  heading = head(fasta, n = 1)[[1]]
  if((start > 30) & (end < (nchar(fs) - 30))){
    sequence = substr(fs, start, end)
  }
  else if((start < 30) & (end < (nchar(fs) - 30))){
    sequence = substr(fs, 1, end)
  }
  else if((start < 30) & (end > (nchar(fs) - 30))){
    sequence = substr(fs, o, nchar(fs))
  }
  else if((start > 30) & (end > (nchar(fs) - 30))){
    sequence = substr(fs, start, nchar(fs))
  }
  write(heading, file = newfastaname)
  write(sequence, file = newfastaname, append = TRUE)

}
