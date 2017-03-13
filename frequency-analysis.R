#library(plotly)

stuff <- list()
for(i in 1:length(blocked_ref)){
  mycol = c()
  for(j in 1:length(blocked_seqs)){
    mycol <- c(mycol,as.character(blocked_seqs[[j]][[1]][[i]]))
  }
  #stuff[[i]] <- mycol
  #print(paste(c("Field: ",i),collapse = ""))
  stuff[[i]] <- base::as.data.frame(base::table(mycol))
  colnames(stuff[[i]]) <- c("sequence_type","count_s")
  p <- ggplot(data = stuff[[i]],aes(x=sequence_type,y =count_s,fill=sequence_type)) +
    geom_bar(stat = "identity",width = 0.5)+xlab("Sequence types")+ylab("counts") +
    ggtitle(base::paste(c("Frequency of different types in field ",i),collapse = "")) +
    theme_bw()
    ggsave(filename = base::paste(
              c("results-blocked-ref-max/blocksize_14/analysis/barplots/","field_",i,"_barplot.png"),collapse = ""
           ),
         plot = p)
}
#barplot(dt$count_s,xlab="Sequence Type",ylab = "Counts",names.arg = dt$sequence_type)
msizes <- c()
for(i in 1:length(stuff)){
 msizes <- c(msizes,nrow(stuff[[i]])) 
}

as.data.frame(table(msizes))

#compute mean absolute deviation from mode
compute_mode <- function(x){
  return(base::which.max(x$count_s))
}
compute_madm <- function(x,n_size,m_mode)
{
  sum_edit <- 0
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
  #Issues with: TGTCACTTCAAATC and TGTCACTTCAAAT - library miscalculation
  for(i in 1:n_size){
    num_edit <- nedit(pairwiseAlignment(DNAString(as.character(x[m_mode,1])),DNAString(as.character(x[i,1])),gapOpening=1, type="global", gapExtension=1,substitutionMatrix = mat))
    num_edit <- num_edit * x[i,2]
    sum_edit <- sum_edit + num_edit
  }
  return((sum_edit*(1/n_size)))
}
varriability <- c()
for(j in 1:length(stuff)){
  mmode <-  compute_mode(stuff[[j]])
  vr <- compute_madm(x = stuff[[j]],n_size = nrow(stuff[[j]]),m_mode =mmode)
  varriability <- c(varriability,vr)
}

as.data.frame(table(varriability))












