include_libs <- function(lib_vectors)
{
  for(i in 1:length(lib_vectors))
  {
    library(lib_vectors[i],character.only = TRUE)
  }
}
par_align <- function(gene,gene_ref)
{
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
  #aligned_bag <<- c(aligned_bag,c(pairwiseAlignment(gene_ref,gene,gapOpening=1, type="global", gapExtension=1,substitutionMatrix = mat)))
  aligned_now <- pairwiseAlignment(gene_ref,gene[[1]],gapOpening=1, type="global", gapExtension=1,substitutionMatrix = mat)
  aligned_now@subject@unaligned@ranges@NAMES <- gene[[2]]
  return (aligned_now)
}
calc_nedits <- function(aligned_seq)
{
  return(list(aligned_seq@pattern@unaligned@ranges@NAMES,aligned_seq@subject@unaligned@ranges@NAMES,nedit(aligned_seq)))
}

calc_blocked_nedits <- function(aligned_seq,aligned_ref)
{
  #aligned_ref is query
  sum_edit <- 0
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
  #print ("Out of first for loop")
  for(i in 1:length(blocked_ref))
  {
    hash.query <- paste(c(as.character(aligned_ref[[1]][i][[1]]),"|",as.character(aligned_seq[[1]][i][[1]])),collapse = "")
    check.in.hash <- has.key(hash.query,h.a.s.h)
    if(check.in.hash==FALSE){
      new.edit.distance <- nedit(pairwiseAlignment(DNAString(aligned_ref[[1]][i][[1]]), 
                                                   DNAString(aligned_seq[[1]][i][[1]]),
                                                   gapOpening=1, type="global", gapExtension=1,substitutionMatrix = mat))
      h.a.s.h[[hash.query]] <- new.edit.distance
      sum_edit <- sum_edit + new.edit.distance
      
    }
    else{
      sum_edit <- sum_edit + h.a.s.h[[hash.query]]
    }
    
  }
  
  return(list(aligned_ref[[2]],aligned_seq[[2]],sum_edit))
}
break_ref_into_blocks <- function(var.seq)
{
  blocks <- list()
  #var.seq <- var.seq[! var.seq %in% c("-")]
  var.seq <- unaligned(var.seq)
  var.seq <- as.character(var.seq)
  
  comp.blocks <- floor(nchar(var.seq)/BLOCK_SIZE)
  incomp.blocks <- nchar(var.seq) - BLOCK_SIZE * comp.blocks
  
  #BLOCKED_LENGTH <- (nchar(var.seq[[1]]) %% BLOCK_SIZE) + 1
  list_block_index <- 1
  
  for(i in 1:comp.blocks){
    blocks[[i]] <- substr(var.seq[[1]],((i-1)*BLOCK_SIZE+1),(i*BLOCK_SIZE))
  }
  if(incomp.blocks != 0){
    blocks[[comp.blocks+1]] <- substr(var.seq[[1]],(comp.blocks*BLOCK_SIZE+1),nchar(var.seq))
  }
  return (blocks)
}

# send aligned(pattern(aligned_seq)) and aligned(subject(aligned_seq))
blockify_sequence <- function(aligned_seq)
{
  aligned_ref <- aligned(pattern(aligned_seq))
  aligned_sub <- aligned(subject(aligned_seq))
  #ref_offset <- 0
  blocked_seq <- list()
  
  previous.pos <- 1
  current.pos <- 1
  
  for(i in 1:length(blocked_ref))
  {
    b <- blocked_ref[[i]]
    #ref_len <- 1
    for(j in 1:nchar(b))
    {
      #new_index <- ref_offset+ref_len
      while (substr(b,j,j) != as.character(aligned_ref[1][[1]][current.pos]))
      {
        current.pos <- current.pos + 1
        #ref_len <- ref_len + 1
        #new_index <- ref_offset+ref_len
      }
      current.pos <- current.pos + 1
      #ref_len <- ref_len + 1
    }
    #new_ref_offset <- ref_offset+ref_len - 1
    #ref_offset <- ref_offset + 1
    blocked_seq[[i]] <- DNAString(gsub(pattern = '-',replacement = '',x = aligned_sub[1][[1]][previous.pos:(current.pos-1)]))
    #ref_offset <-new_ref_offset + 1
    previous.pos <- current.pos
  }
  return(list(blocked_seq,aligned_seq@subject@unaligned@ranges@NAMES))
  
}








