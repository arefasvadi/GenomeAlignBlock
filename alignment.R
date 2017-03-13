rm(list=ls())
setwd("~/projects/GenomeAlignBlock/R scripts")
included_libs <- c("data.table","Biostrings","BiocParallel")
source("general_functions.R")

include_libs(included_libs)
input.genomes <- "genomes-nodup.fasta"

raw.genomes <- readDNAStringSet(input.genomes)
#aligned_bag <- c()
#select maximum size sequence as reference
#geneome.max <- raw.genomes[which.min(raw.genomes@ranges@width)]
extended_genome <- list()
for(i in 1:length(raw.genomes)){
  extended_genome[[i]] <- list(raw.genomes[i],raw.genomes[i]@ranges@NAMES)
}
for(l in 1:length(raw.genomes)){
  #mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
  results <- bplapply(X = extended_genome,FUN =  par_align,raw.genomes[l],BPPARAM = MulticoreParam())
  #aligntest <- pairwiseAlignment(raw.genomes[1],raw.genomes[2],gapOpening=1, type="global", gapExtension=1,substitutionMatrix = mat)
  edits <- bplapply(X =  results,FUN =  calc_nedits,BPPARAM = MulticoreParam())
  edits <- rbindlist(edits)
  colnames(edits) <- c("query_id","data_id","distance")
  edits <- edits[base::order(edits$distance),]
  query_id <- as.character(edits[1,1])
  write.csv(edits,file = paste(c("results/query_id_",query_id,".csv"),collapse = ""))
}