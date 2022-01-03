#' Compute ranks of rows of matrix and summarize them into a choice suggestion.
#' 
#' This function allows getting out the choice (of a column representing a stock) from four
#' rows of numbers quantifying the four orders of exact stochastic dominance comparisons.
#' If the last or 10-th row for ``choice" has 1 then the stock representing
#' that column is to be chosen. That is it should get the largest
#' (portfolio) weight. If the original matrix row names are SD1 to SD4, 
#' the same names are repeated for the extra rows representing their ranks.  
#' The row name for ``sum of ranks" is
#' sumRanks. Finally, the ranks associated with sumRanks provide the row named choice
#' along the bottom (10-th) row of the output matrix called "out." 
#' 
#' @param mtx {matrix to be ranked by row and summarized}
#' @return a matrix called `out' having 10 rows and p columns (p=No.of stocks). 
#' Row Numbers 1 to 4 have SD1 to SD4 evaluation of areas over ECDFs. 
#' There are 6 more rows. Row No.5= SD1 ranks,
#' Row No.6= SD2 ranks, Row No.7= SD3 ranks, Row No.8= SD4 ranks
#' Row No.9= sum of the ranks in earlier four rows for ranks of SD1 to SD4
#' Row No.10= choice rank based on all four (SD1 to SD4) added together
#' Thus, the tenth row yields choice priority number for each stock (asset)
#' after combining the all four criteria.
### @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso \code{\link{exactSdMtx}}
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
