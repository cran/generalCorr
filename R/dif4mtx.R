#' order four differencing of a matrix of time series 
#' 
#' This is for momentum traders who focus on growth, acceleration, its growth
#' and further acceleration. The diff function of R seems to do 
#' recycling of available numbers, not wanted for our purposes. Hence, this
#' function is needed in portfolio studies based on time series.
#' 
#' @param mtx {(n X p) matrix of p time series (market returns) with n items each}
#' @return out matrix having 12 rows, (data, D1 to D4
#' and ranks of D1 to D4
#' The column names of out are those of input matrix mtx.
### @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
### @seealso \code{\link{}}
#' @examples
#' x=c(2,8,3,5,1,8,19,22,23)
#' y=c(3,11,2,6,7,9,20,25,21)
#' dif4mtx(cbind(x,y))
#' 
#' @export

dif4mtx=function(mtx){
n=NROW(mtx)
p=NCOL(mtx)
out=matrix(NA,nrow=5, ncol=p)
 for (j in 1:p){
  d5=dif4(mtx[,j])  
  a1=apply(d5,2,sum)
  out[,j]=a1
 }#end j loop
  rownames(out)=colnames(d5)
  colnames(out)=colnames(mtx)
  s1=summaryRank(out)  
return(s1)}


