#' Two sequences: starting+ending values from n and blocksize (internal use)
#' 
#' This is an auxiliary function for gmcmtxBlk. It gives sequences of starting
#' and ending values 
#' 
#' @param n {length of the range}
#' @param blksiz {blocksize}
## @importFrom stats seq
#' @return two vectors sqLO and sqUP
### @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso \code{\link{gmcmtxBlk}}
#' @examples
#' 
#' getSeq(n=99, blksiz=10)
#' 
 
#' @export


getSeq=function(n, blksiz){
  k=floor(n/blksiz)
  sqLO=seq(from=1, length.out=k, by=blksiz)
  sqUP=c(seq(from=blksiz, by=blksiz, length.out=k-1),n)
  list(sqLO=sqLO, sqUP=sqUP)
}
