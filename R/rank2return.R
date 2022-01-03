#' Compute the portfolio return knowing the rank of a stock in the input `mtx'.
#' 
#' This function computes the return earned knowing the rank of a stock in
#' the input mtx of stock returns. For example, mtx has p=28 Dow Jones stocks
#' over n=169 monthly returns. Portfolio weights are assumed to be linearly
#' declining. If maxChosen=4, the weights are 1/10, 2/10, 3/10 and 4/10, which add
#' up to unity. These portfolio weights are assigned in reverse order
#' in the sense that first chosen stock (choice rank =1) gets portfolio weight=4/10.
#' The function computes return from the stocks by the `myrank' argument.
#' 
#' @param mtx {a matrix with n rows (number of returns) p columns (number of stocks)}
#' @param myrank {vector of p integers listing the rank of each stock, 1=best}
#' @param maxChosen {number of stocks in the portfolio (with nonzero weights)
#' default=0. When maxChosen=0, we let pctChoose determine the maxChosen}
#' @param pctChoose {percent of p stocks chosen inside the portfolio, default=20}
#' @param verbo {logical if TRUE, print, default=TRUE}
## @importFrom stats seq
#' @return average return from the linearly declining
#'  portfolio implied by the myrank vector.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso \code{\link{outOFsamp}}
#' @export

rank2return = function(mtx,myrank,maxChosen=0, pctChoose=20, verbo=FALSE){
  nn=NROW(mtx)
  pp=NCOL(mtx)
  p=pp
if(verbo) print(c("rank2return: mtx dimensions n,p",nn,pp))
if(maxChosen==0) maxChosen=round(pctChoose*p/100,0)
if(verbo){
    print("inptus to rank2return function")
    print(c("maxChosen",maxChosen))
    print(c("myrank",myrank))  }#end verbo if
if(length(myrank)<p) stop("invalid myrank input")
#mtxw=mtx #initial place holder
  wtx=rep(0,p) #initially all zero weights
if(maxChosen>p) print(c("maxChosen exceeds p",maxChosen,p))
if(maxChosen>p) stop("bad maxChosen")
  wt=(maxChosen:1)/sum(1:maxChosen) #nonzero weights
  wtx=c(wt,rep(0,p-maxChosen)) #length(wtx)=p
  roundMyrank=round(myrank[1:p],0) #force integer ranks
  wtx2=wtx[roundMyrank]#replace wtx by (linear) weight
#  print(wtx2)
  sumWtx2=sum(wtx2)
  if (sumWtx2 != 1) wtx2=wtx2/sumWtx2  #sum(wtx2)=1 forced
  if(verbo) print(c("nonzero weights",wt))
  if(verbo) print(c("weights wtx replaced",wtx2))
  if (pp>=2) avg=apply(as.matrix(mtx)%*%as.matrix(diag(wtx2)),2,mean,na.rm=TRUE)
  if (pp==1) avg=mean(mtx*as.matrix(wtx2),na.rm=TRUE)
  return(sum(avg))
}#end rank2return function
