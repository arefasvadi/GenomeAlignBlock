rm(list=ls())
setwd("~/projects/GenomeAlignBlock/R scripts")
inputpath <- "result-35-3477-max.csv"
included_libs <- c("ggplot2","ggiraph","dplyr","tidyr","stringr","rgl","plotly","data.table","Biostrings")

include_libs(included_libs)

read_aligned_csv <- function(inputpath)
{
   return (read.csv(inputpath))
  
}

aligned_table <- read_aligned_csv(inputpath)
seq_size <- rep(c(0),times=length(aligned_table[,1]))
cbind(aligned_table,seq_size)
ref_size <- rep(c(0),times=length(aligned_table[,1]))
cbind(aligned_table,ref_size)

for(i in 1:length(aligned_table[,1])){
  aligned_table[i,"seq_size"] = nchar(as.vector(aligned_table[i,"seq_aligned"]))
  aligned_table[i,"ref_size"] = nchar(as.vector(aligned_table[i,"ref_aligned"]))
}

for(i in 1:length(aligned_table[,1])){
  if(aligned_table[i,"seq_size"] != aligned_table[i,"ref_size"]){
    print(paste(c("Discrepencey at row ",i)))
  }
}

grp_seq <- aligned_table %>% group_by(seq_aligned) %>% count() %>% arrange(desc(n))
grp_ref <- aligned_table %>% group_by(ref_aligned) %>% count() %>% arrange(desc(n))

write.csv(aligned_table,file = "max2-aligned-table.csv")
write.csv(grp_seq,file = "max2-grp-seq.csv")
write.csv(grp_ref,file = "max2-grp-ref.csv")