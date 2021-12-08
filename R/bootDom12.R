#' bootstrap confidence intervals for exact stochastic dominance SD1 to SD4.
#' 
#' This calls the meboot package to create J=999 replications of portfolio return
#' matrices and compute 95\% confidence intervals on x1, x2 and their
#' difference (x2-x1).  If the interval on (x2-x1) contains zero the choice
#' between the two can reverse due to sampling variation
#' 
#' @param x1 {a vector of n portfolio returns}
#' @param x2 {a vector of n portfolio returns}
#' @param confLevel {confidene level confLevel=95 is default}
#' @param reps {number of bootstrap resamples, default is reps=999}
## @importFrom stats seq
#' @return a matrix with six columns. First two Low1 and Upp1
#' are confidence interval limits for x1. Next two columns
#' have analogous limits for x2. The last but first columns entitled
#' Lowx2mx1 means lower confidence limit for (x2-x1), where m=minus.
#' The last column entitled
#' Uppx2mx1  means upper confidence limit for (x2-x1). 
#' 
#' For strong stochastic dominance of x2 over x1
#' dominance beyond sampling variability, zero should not be inside
#' the confidence interval in the last two columns. 
#' 
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso see \code{\link{exactSdMtx}}
#' @importFrom meboot meboot
#' @importFrom stats quantile
#' 
#' @export

bootDom12=function(x1,x2,confLevel=95,reps=999){
#  require(meboot)
  alp=(100-confLevel)/100
  x1boot=meboot(x1,reps=reps)$ens
  x2boot=meboot(x2,reps=reps)$ens
  x1conf=matrix(NA,nrow=reps,ncol=4)
  x2conf=matrix(NA,nrow=reps,ncol=4)
  for(j in 1:reps){
    x1b=x1boot[,j]
    x2b=x2boot[,j]
    s12=exactSdMtx(cbind(x1b,x2b))
    s12x=s12$out
  #  if (j<3) print(s12x)
    x1conf[j,1:4]=s12x[,1]
    x2conf[j,1:4]=s12x[,2]
  }#end j loop over reps
  alpby2=alp/2
  oneminus=1-alpby2
  Fn1=function(x) quantile(x,probs=alpby2)
  Fn2=function(x) quantile(x,probs=oneminus)
  Low1=apply(x1conf,2,Fn1)
  Upp1=apply(x1conf,2,Fn2)
  CI1=cbind(Low1,Upp1)
  #print("confidence interval for x2")
  Low2=apply(x2conf,2,Fn1)
  Upp2=apply(x2conf,2,Fn2)
  CI2=cbind(Low2,Upp2)
  #print("confidence interval for difference x2-x1")
  Lowx2mx1=apply((x2conf-x1conf),2,Fn1)
  Uppx2mx1=apply((x2conf-x1conf),2,Fn2)
  CIx2mx1=cbind(Lowx2mx1,Uppx2mx1)
  #print(CIx2mx1)
  CIall=cbind(CI1,CI2,CIx2mx1)
  rownames(CIall)=c("SD1","SD2","SD3","SD4")
  return(CIall)  
}#end function bootDom12 