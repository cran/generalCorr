#' Exact stochastic dominance computation from areas above ECDFs compared to 
#' common reference minimum (refmin) for order SD1 to SD4 
#' 
#' This function inputs 'mtx' (n X p) matrix data (e.g. monthly returns on p stocks)
#' Its output has four matrices SD1 to SD4 each with dimension (n X p). 
#' 
#' 
#' @param mtx {(n X p) matrix of data. For example, returns on p stocks n months}
#' @param howManySd {used to define (x.ref)= lowest return number by going
#'  default=0.1 howManySd =number of maximum standard deviations in the data below 
#'  the all-row all-column minimum return in the data}
## @importFrom stats seq
#' @return four matrices SD1 to SD4 for four order of stochastic dominance 
#' over a common x.ref.  The "out" matrix is another output having 4 rows for
#' SD1 to SD4 and p columns (p=No. of columns in data matrix mtx).
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
#    mtxmax=max(apply(mtx,2,max))
#    fixmax=mtxmax+howManySd*mtxsd
    refmin=minval-howManySd*maxsd
    # always print  
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
    }#end j loop
    out=matrix(NA,nrow=4, ncol=p)
    out[1,]=apply(SD1,2,sum,na.rm=TRUE)
    out[2,]=apply(SD2,2,sum,na.rm=TRUE)
    out[3,]=apply(SD3,2,sum,na.rm=TRUE)
    out[4,]=apply(SD4,2,sum,na.rm=TRUE)
    rownames(out)=c("SD1","SD2","SD3","SD4")
    list(SD1=SD1,SD2=SD2,SD3=SD3,SD4=SD4,out=out)
} #end function
