#' Function to do matched deletion of missing rows from x, y and control variable(s). 
#' 
#' The aim in three-way deletions is to retain only the largest 
#' number of available data triplets with all non-missing data.
#' 
#' @param x Vector of x data
#' @param y Vector of y data
#' @param ctrl {Data matrix on the control variable(s) kept beyond causal path determinations}
#' @return 
#' \item{newx}{A new vector x after removing triplet-wise missing data} 
#' \item{newy}{A new vector or matrix y after removing triplet-wise missing data} 
#' \item{newctrl}{A new vector or matrix ctrl after removing triplet-wise missing data} 
## @note %% ~~further notes~~
#' @importFrom stats complete.cases
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See \code{\link{napair}}.
#' @examples
#' 
#' \dontrun{
#' x=sample(1:10);y=sample(1:10);x[2]=NA; y[3]=NA
#' w=sample(2:11)
#' naTriplet(x,y,w)}
#' 
#' @export


naTriplet = function(x, y, ctrl) {
  # ctrl is a matrix of control variables
  p=NCOL(y)
  pc=NCOL(ctrl)
  len=length(ctrl)
  if(len==1) {  newctrl=0
    ok=complete.cases(x,y)
    newx = x[ok]
    if(p==1)newy = y[ok] 
    if(p>1)newy=y[ok,]
    }  #delete NAs from x and y   
    if(len>1) {
    ok=complete.cases(x,y,ctrl)
    newx = x[ok]
    if(p==1)newy = y[ok] 
    if(p>1) newy=y[ok,]
    if (pc==1) newctrl=ctrl[ok]
    if(pc>1) newctrl=ctrl[ok,] }  
  #delete NAs from x, y, ctrl   
  list(newx = newx, newy = newy, newctrl = newctrl)  #delete NAs
}

