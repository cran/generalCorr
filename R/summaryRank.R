#' Compute ranks of rows of matrix and summarize them into a choice suggestion.
#' 
#' This function allows getting out the choice (of column stock) from four
#' rows of numbers quantifying the four orders of exact stochastic dominance comparisons.
#' If the choice row has 1 then that column is to be chosen with the largest
#' (portfolio) weight. If the original matrix row names exist, same are repeated
#' for the extra rows representing the ranks.  The row name for sum of ranks is
#' sumRanks and the ranking of the sumRanks providing the choice is named choice
#' along the bottom row of the output matrix called "out." 
#' 
#' @param mtx {matrix to be ranked by row and summarized}
#' @return a matrix appending 4 more rows of ranks and a couple of rows
#' summarizing ranks and also choice priority number for each stock (asset)
### @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
### @seealso \code{\link{}}
#' 
#' @export

summaryRank=function(mtx){
  n=NROW(mtx)
  p=NCOL(mtx)
  mtxrank=mtx #place to store
  for (i in 1:n){mtxrank[i,]=(p+1)-rank(mtx[i,],
               na.last = TRUE,  ties.method = "average")}
  sumRanks=apply(mtxrank,2,sum)
  choice=rank(sumRanks,na.last=TRUE,ties.method = "average")
  out=rbind(mtx,mtxrank,sumRanks,choice)
  return(out)
}#end summaryRank function
