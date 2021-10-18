#' Replace asymmetric matrix by max of abs values of [ij] or [ji] elements
#' useful in symmetrizing gmcmtx0 general correlation matrix
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
