#' Replace asymmetric matrix by max of abs values of [i,j] or [j,i] elements.
#' 
#' It is useful in symmetrizing the gmcmtx0 matrix containing a non-symmetric
#' generalized correlation matrix.
#'
#'
#' @param mtx {non-symmetric matrix}
#' @return 
#' \item{mtx2}{replace [i,j] and [j,i] by the max of absolute values
#' with common sign}
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @examples 
#' 
#' \dontrun{
#' example 
#'mtx=matrix(1:16,nrow=4)
#'symmze(mtx)
#'}#' 
#' @export
 
symmze=function(mtx){
  nr=NROW(mtx)  
  nc=NCOL(mtx)
  mtx2=mtx
  if(nc!=nr) stop("stop: symmze(mtx) error: must input a square matrix")
  for (i in 1:nr) {
    for (j in 1:nc){
      if(i!=j) {
        aij=abs(mtx[i,j])
        aji=abs(mtx[j,i])
        if (aij>=aji) mtx2[i,j]=sign(mtx[i,j])*aij
        if (aji>aij) mtx2[i,j]=sign(mtx[j,i])*aji
      }#eindif
    }#end j loop    
  }# end i loop
  return(mtx2)
}# end function
