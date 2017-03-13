#rm(list = ls())
rankify <- function(path.file)
{
    dat <- read.csv(path.file)
    #dat <- dat[-1,-1]
    dat <- dat[,-1]
    dat <- dat[which(dat[,1] != dat[,2]),]
    dat$rank <- rank(dat$distance,ties.method = "min")
    return(dat)
}

path.file <- "results-gt/"
tempfiles <- list.files(path = path.file,pattern = "*.csv")
expr.file <- "results-blocked-ref-max/blocksize_14/"
expr.tempfiles <- list.files(path = expr.file,pattern = "*.csv")
ground.truth.files <- list()
new.approach.files <- list()
for (i in 1:length(tempfiles))
{
  ground.truth.files[[i]] <- rankify(paste(c(path.file,tempfiles[i]),collapse = ""))
  new.approach.files[[i]] <- rankify(paste(c(expr.file,expr.tempfiles[i]),collapse = ""))
}
if(length(intersect(as.vector(tempfiles),
                    as.vector(expr.tempfiles)))
   !=max(length(tempfiles),length(expr.tempfiles)))
{
  print("Files are not the same!")
  stop()
}
#Top K Processing
output.frame <- data.frame(k=integer(0),query_id=integer(0),t_r_size=integer(0),a_r_size=integer(0),num.correct=integer(0),
                           recall.ratio.k=numeric(0),recall.ratio.t=numeric(0),recall.ratio.e=numeric(0),
                          stringsAsFactors = FALSE)
colnames(output.frame) <- c("K","Query ID","True Size","Estimated Size", "#Correct", "Recall K","Recall T_Size","Recall E_Size")

for(k in 1:20)
{
  for(i in 1:length(ground.truth.files))
  {
    ground.truth <- ground.truth.files[[i]]
    new.approach <- new.approach.files[[i]]
    
    if(!all.equal(ground.truth[,1],new.approach[,1]))
    {
      print("Query IDs are not the same!")
      stop()
    }
     
    #pick first k from ground truth
    ground.truth.k <- ground.truth[1:k,]
    l <- 1;
    while(ground.truth[(k+l),4] == ground.truth[(k+l-1),4])
    {
      ground.truth.k[(k+l),] <- ground.truth[(k+l),]
      l <- l+1
    }
    
    #pick first k from ground truth
    new.approach.k <- new.approach[1:k,]
    l <- 1;
    while(new.approach[(k+l),4] == new.approach[(k+l-1),4])
    {
      new.approach.k[(k+l),] <- new.approach[(k+l),]
      l <- l+1
    }
    
    #now check for common results
    num.intersect <- length(intersect(as.vector(new.approach.k[,2]),as.vector(ground.truth.k[,2])))
    recall.ratio.k <- num.intersect / k
    recall.ratio.esize <- num.intersect / nrow(new.approach.k)
    recall.ratio.tsize <- num.intersect / nrow(ground.truth.k)
    
    new.row <- nrow(output.frame)+1
    output.frame[new.row,"K"] <- k
    output.frame[new.row,"Query ID"] <- ground.truth[1,1]
    output.frame[new.row,"True Size"] <- nrow(ground.truth.k)
    output.frame[new.row,"Estimated Size"] <- nrow(new.approach.k)
    output.frame[new.row,"Recall K"] <- recall.ratio.k
    output.frame[new.row,"Recall T_Size"] <- recall.ratio.tsize
    output.frame[new.row,"Recall E_Size"] <- recall.ratio.esize
    output.frame[new.row,"#Correct"] <- num.intersect
  }
}
#write.csv(output.frame[,1:base::ncol(output.frame)],file = "results-blocked-ref-857/blocksize_10/analysis/analysis.csv")
graph.framee <- data.frame(b=integer(0), k=integer(0),recall.ratio.e=numeric(0),
                           stringsAsFactors = FALSE)
colnames(graph.framee) <- c("b","k","r")

for(i in 1:length(base::unique(output.frame[,1])))
{
  #print(paste(c(i,",",mean(output.frame[which(output.frame$K==i),"Recall E_Size"])),collapse = ""))
  rn <- nrow(graph.framee)+1
  graph.framee[rn,"b"] = "18"
  graph.framee[rn,"k"] = i
  graph.framee[rn,"r"] = mean(output.frame[which(output.frame$K==i),"Recall K"])
  #plt <- 
  #ggplotly(plt)
}
 p <- ggplot(data = graph.framee,aes(x=k,y=r,group=b,colour=b)) +
  geom_line() +
  geom_point() +
  theme_bw()
#(which(output.frame$K==5 & output.frame$`#Correct`==0))
#length((which(output.frame$K==11 & output.frame$`#Correct` >= 10 )))

#outlist <- list()
#outlist[["blocksize10"]] <- output.frame

























