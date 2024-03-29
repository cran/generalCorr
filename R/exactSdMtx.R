#' Exact stochastic dominance computation from areas above ECDF pillars.
#' 
#' ECDF=empirical cumulative distribution functions. These are sufficient
#' statistics representing probability density functions
#' defined by observable finite data (e.g., stock returns). The exact computation
#' of stochastic dominance orders SD1 to SD4 needs areas between two ECDFs,
#' since such areas represent integrals. Higher-order SDs with continuous 
#' variables involve repeated integrals. Our quantification needs
#' areas of ECDFs defined from areas of lower-order ECDFs. We argue that these
#' computations are convenient if there is an ECDF of an imaginary
#' reference minimum (x.ref) return, whose ECDF is a rectangle 
#' common for all stock comparisons. A common (x.ref) avoids having to compute
#' all possible pairs of p stocks. Choosing a common reference as SP500 index
#' stock cannot avoid a slower trapezoidal approximation for integrals,
#' since its returns vary over time. We want exact areas of rectangles and fast.
#' 
#' The \code{exactSdMtx} function inputs `mtx' (n X p) matrix data 
#' (e.g., n monthly returns on p stocks).
#' Its output has four matrices SD1 to SD4, each with dimension (n X p). They measure
#' exact dominance areas between empirical CDF for each column to the ECDF of
#' (x.ref), an artificial stock with minimal return in all time periods. A fifth
#' output matrix called `out' produced by \code{exactSdMtx}
#' has 4 rows and p columns containing column sums of SD1 to SD4. 
#' We intend that this
#' `out' matrix produced by \code{exactSdMtx} is then input to another
#' function \code{summaryRank()} in the package designed for practitioners.
#' For example, it indicates the best and the worst columns 
#' representing (the best stock to buy and best stock to sell)
#' from the input data `mtx' for investment based on a sophisticated computation
#' of their ranks. 
#' 
#' 
#' 
#' 
#' @param mtx {(n X p) matrix of data. For example, returns on p stocks
#' over n months}
#' @param howManySd {used to define (x.ref)= lowest return number.
#'  If the grand minimum of all returns in `mtx' is denoted GrMin, then
#'  howManySd equals the number of max(sd) (maximum standard deviation for data
#'  columns) below the GrMin used to define (x.ref). Thus,
#'  (x.ref)=GrMin-howManySd*max(sd). default howManySd=0.1 }
## @importFrom stats seq
#' @return five matrices. SD1 to SD4 contain four orders of stochastic 
#' dominance areas using the ECDF pillars and 
#' a common (x.ref). The fifth "out" matrix is another output with 4 rows for
#' SD1 to SD4, and p columns (p=No. of columns in data matrix mtx) having a 
#' summary of ranks using all four, SD1 to SD4.
### @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
### @seealso \code{\link{}}
#' @examples
#' x1=c(2,5,6,9,13,18,21)
#' x2=c(3,6,9,12,14,19,27) 
#' st1=exactSdMtx(cbind(x1,x2))
#' 
#' 
#' @export

exactSdMtx=function(mtx, howManySd=0.1){
  p=NCOL(mtx)
  n=NROW(mtx)
if (n<5) stop("stop n<5 input matrix to exactStDomMtx")
if (p<2) stop("stop p<2 input matrix to exactStDomMtx")
    maxsd=max(apply(mtx,2,sd,na.rm=TRUE),na.rm=TRUE)#largest sd
    minval=min(apply(mtx,2,min,na.rm=TRUE),na.rm=TRUE)#min value in data
    refmin=minval-howManySd*maxsd
#    print("howManySd,maxsd,minval,refmin, next line")
#    print(c(howManySd,maxsd,minval,refmin))
    SD1=mtx #place to store
    SD2=mtx
    SD3=mtx
    SD4=mtx #place to store  SD1 computations next
  for (j in 1:p){
    x=mtx[,j] #pick jth col of data
    sx=sort(c(refmin,x),na.last=TRUE)
    height=rev(1:n)/n  #should have n heights of pillars above ECDFs
    wid=diff(sx)
    ara=height*wid
    SD1[,j]=ara
  }#end j loop for SD1,   SD2 computations next
    for (j in 1:p){
      x= SD1[,j]#pick jth col of SD1
      sx=sort(c(0,x),na.last=TRUE)
      height=rev(1:n)/n  #should have n heights of pillars above ECDFs
      wid=diff(sx)
      ara=height*wid
      SD2[,j]=ara
    }#end j loopfor SD2
    for (j in 1:p){
      x= SD2[,j]#pick jth col of SD2
      sx=sort(c(0,x),na.last=TRUE)
      height=rev(1:n)/n  #should have n heights of pillars above ECDFs
      wid=diff(sx)
      ara=height*wid
      SD3[,j]=ara
    }#end j loop, 
    for (j in 1:p){
      x= SD3[,j]#pick jth col of SD3
      sx=sort(c(0,x),na.last=TRUE)
      height=rev(1:n)/n  #should have n heights of pillars above ECDFs
      wid=diff(sx)
      ara=height*wid
      SD4[,j]=ara
    }#end j loop for SD4 using SD3
    out=matrix(NA,nrow=4, ncol=p)
    out[1,]=apply(SD1,2,sum,na.rm=TRUE)
    out[2,]=apply(SD2,2,sum,na.rm=TRUE)
    out[3,]=apply(SD3,2,sum,na.rm=TRUE)
    out[4,]=apply(SD4,2,sum,na.rm=TRUE)
    rownames(out)=c("SD1","SD2","SD3","SD4")
    list(SD1=SD1,SD2=SD2,SD3=SD3,SD4=SD4,out=out)
} #end function
