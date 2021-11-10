#' order 4 differencing of a time series vector
#' 
#' This is for momentum traders who focus on growth, acceleration, its gorwth
#' and further acceleration. The diff function of R seems to do 
#' recycling of available numbers, not wanted for our purposes.
#' 
#' @param x {(n X 1) vector of time series (market returns) with n items each}
#' @return ou2 matrix having five columns, first for x, the next four
#' columns have diff(x), diff-squared(x), diff-cubed(x) and diff-fourth(x)
### @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
### @seealso \code{\link{}}
#' @examples
#' x=c(2,8,3,5,1,8,19,22,23)
#' dif4(x)
#' 
#' @export


dif4=function(x){
  xa=c(rep(NA,4),x)
d1=diff(xa)
d2=diff(d1)
d3=diff(d2)
d4=diff(d3)
out=cbind(xa[5:length(xa)],d1[4:length(d1)],
          d2[3:length(d2)],
          d3[2:length(d3)],
          d4[1:length(d4)])
ou2=out[5:NROW(out),]
colnames(ou2)=c("x","D1","D2","D3","D4")
return(ou2)}


