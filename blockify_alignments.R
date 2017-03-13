rm(list=ls())
setwd("~/projects/GenomeAlignBlock/R scripts")
included_libs <- c("data.table","Biostrings","BiocParallel","hash")
source("general_functions.R")
BLOCK_SIZE <- 14
h.a.s.h <- hash()

include_libs(included_libs)
input.genomes <- "genomes-nodup.fasta"

raw.genomes <- readDNAStringSet(input.genomes)

extended_genome <- list()
for(i in 1:length(raw.genomes)){
  extended_genome[[i]] <- list(raw.genomes[i],raw.genomes[i]@ranges@NAMES)
}
#select maximum size sequence as reference
genome.max <- raw.genomes[which.max(raw.genomes@ranges@width)]
#genome.max <- raw.genomes[1]
results <- bplapply(X = extended_genome,FUN=par_align,genome.max ,BPPARAM = MulticoreParam())
blocked_ref <- break_ref_into_blocks(pattern(results[[1]]))

#blocked_seqs <- list()
#for(i in 1:length(results))
#{
#  blocked_seqs[[i]] <- blockify_sequence(aligned_ref = aligned(pattern(results[[i]])),aligned_sub = aligned(subject(results[[i]])))
#}
blocked_seqs <- bplapply(X = results,FUN = blockify_sequence)
#rm(results)

for(l in 1:length(raw.genomes)){
#for(l in 122:length(raw.genomes)){
  #l <- 1 
  edits <- bplapply(X = blocked_seqs,FUN =  calc_blocked_nedits, blocked_seqs[[l]],BPPARAM = MulticoreParam())
  edits <- rbindlist(edits)
  colnames(edits) <- c("query_id","data_id","distance")
  edits <- edits[base::order(edits$distance),]
  query_id <- as.character(edits[1,1])
  write.csv(edits,file = paste(c("results-blocked-ref-max/blocksize_14/query_id_",query_id,".csv"),collapse = ""))
}