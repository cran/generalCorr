#' Function compares nine deciles of a matrix to reference minimum (eg p stock returns) 
#' 
#' The first step computes a minimum reference return and  nine deciles.
#' The input x must be a matrix having p columns (col.names recommended).
#' and n rows as in the data.  If data are missing for some columns, insert NA's.
#' Thus x has p column of data matrix ready for comparison and ranking. 
#' For example, x has a matrix of stock returns.
#' The output matrix produced by this function also has p columns for each 
#' column (i.e. for each stock being compared). The output matrix has
#' nineteen rows. Top nine rows have the magnitudes of deciles.
#' rows 10 to 18 have respective ranks of the decile magnitudes. 
#' The output final row  reports a weighted sum
#' of ranks.  Ranking always gives the smallest number 1 to the most desirable outcome.  
#' The 19th line of the
#' output matrix has weighted sum of ranks and we suggest higher portfolio weight
#' be given to the column having smallest value (in the 19th line)
#' or chooing in the order of numbers in the last (20th) line of out matrix.
#' 
#' @param mtx { (n X p) matrix of data. For example, returns on p stocks n months}
#' @param howManySd {used to define fixmin= imaginary lowest return defined by going
#'   howManySd=default=0.1 maximum standard deviations of all stockss below 
#'  the minimum return for all stocks in the data}
### @importFrom moments skewness kurtosis
#' @return out is a matrix with p columns (same as in the input matrix x) and
#' twenty rows. Top nine rows have deciels quantities, next nine are their ranks.
#' The 19-th row of out has weighted sum of ranks. All compared to the minimum rank
#' portfolio manager can give largest weight to the column with smalles bottom line.
#' The bottom line (20-th) labelled choice of output matrix suggests 
#' portfolio manager give the highest weight to the investment option in
#' the column having number 1 and the lowest weight (=0?) to the column 
#' having number p as the choice number. 
#' Another output of the function is fixmin calculated in this function.
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
