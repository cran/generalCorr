#' Function compares nine deciles of stock return distributions. 
#' 
#' The first step computes a minimum reference return and  nine deciles.
#' The input x must be a matrix having p columns (with a name for each column)
#' and n rows as in the data.  If data are missing for some columns, insert NA's.
#' Thus x has p column of the data matrix ready for comparison and ranking. 
#' For example, x has a matrix of stock returns.
#' The output matrix produced by this function also has p columns for each 
#' column (i.e., for each stock being compared). The output matrix has
#' nineteen rows. The top nine rows have the magnitudes of deciles.
#' Rows 10 to 18 have respective ranks of the decile magnitudes. 
#' The next (19-th) row  of the output reports a weighted sum
#' of ranks.  
#' Ranking always gives the smallest number 1 to the most desirable outcome.  
#' We suggest that a higher portfolio weight
#' be given to the column having smallest rank value (along the 19th line).
#' The 20-th row further ranks the weighted sums of ranks in row 19. Investor
#' should choose the stock (column) representing the smallest rank 
#' value along the last (20th) row of the `out' matrix.
#' 
#' @param mtx { (n X p) matrix of data. For example, returns on p stocks n months}
#' @param howManySd {used to define `fixmin'= imaginary lowest return defined by going
#'   howManySd=default=0.1 maximum of standard deviations of all stocks below 
#'  the minimum return for all stocks in the data}
### @importFrom moments skewness kurtosis
#' @return out is a matrix with p columns (same as in the input matrix) and
#' twenty rows. Top nine rows have 9 deciles, next nine rows have their ranks.
#' The 19-th row of `out' has a weighted sum of 9 ranks. All columns refer to
#' one stock. The weighted sum for each stock is then ranked. A
#' portfolio manager is assumed to prefer higher return represented by
#' high decile values represented by the column with the largest weighted sum. 
#' can give largest weight to the column with the smallest bottom line.
#' The bottom line (20-th) labeled ``choice" of the `out' matrix is
#' defined so that choice =1 suggests the stock deserving 
#' the highest weight in the portfolio. The portfolio manager will
#' generally give the lowest weight (=0?) to the stock representing column 
#' having number p as the choice number. The manager may want to sell this stock. 
#' Another output of the `decileVote' function is `fixmin' representing the
#' smallest possible return of all the stocks in the input `mtx' of returns.
#' It is useful as a reference stock. We compute stochastic dominance numbers
#' for each stock with this imaginary stock yielding fixmin return for all time
#' periods.
#' 
### @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
### @seealso \code{\link{}}
#' @examples
#' x1=c(1,4,7,2,6)
#' x2=c(3,4,8,4,7)
#' decileVote(cbind(x1,x2))
#' 
#' @export

decileVote=function(mtx, howManySd=0.1){
  n=NROW(mtx)
  p=NCOL(mtx)
  #compare deciles to zero return investment
  pr=1:9/10
  allmin=apply(mtx,2,min,na.rm=TRUE)
  allsd=apply(mtx,2,sd,na.rm=TRUE)
  maxsd=max(allsd,na.rm=TRUE)
  fixmin=min(allmin,na.rm=TRUE)-howManySd*maxsd
# print(c("reference.min for all p columns",fixmin))
#  q1=quantile(fixmin,probs=pr)
#  print(q1)
  #q1=rep(fixmin,p)
  q2=apply(mtx,2,quantile,probs=pr,na.rm=TRUE)
#  print(q1)
#  print(q2)
  decdiff=q2 #place to store
  decRank=q2 #place to store
  for (j in 1:p) {decdiff[,j]=q2[,j]-fixmin}
  for (i in 1:9){decRank[i,]=(p+1)-rank(decdiff[i,],
                na.last = TRUE,  ties.method = "average")}
  colsum=apply(decRank,2,sum)
  ord1=rank(colsum,na.last=TRUE)
  choice=ord1
  out=rbind(q2,decRank,colsum,choice)
  list(out=out,fixmin=fixmin)
}
