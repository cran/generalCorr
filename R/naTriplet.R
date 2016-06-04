#' Function to do matdched deletion of missing rows from x, y and control variable(s). 
#' 
#' The aim in three-way deletions is to retain only the largest 
#' number of available data triplets with all non-missing data.
#' 
#' @param x Vector of x data
#' @param y Vector of y data
#' @param ctrl {Data matrix on the control variable(s) beyond causal path issues}
#' @return 
#' \item{newx}{A new vector x after removing pairwise missing data} 
#' \item{newy}{A new vector y after removing pairwise missing data} 
## @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See \code{\link{napair}}.
#' @examples
#' 
#' \dontrun{
#' x=sample(1:10);y=sample(1:10);x[2]=NA; y[3]=NA
#' napair(x,y)}
#' 
#' @export


naTriplet = function(x, y, ctrl) {
  # ctrl is a mtrix of control variables
  nk = NCOL(ctrl)  #number of control variables in ctrl matrix
  ava.x = which(!is.na(x))  #ava means available
  ava.y = which(!is.na(y))  #ava means non-missing
  ava3 = intersect(ava.x, ava.y)
  if (nk == 1) 
    ava.k = which(!is.na(ctrl))
  if (nk > 1) {
    for (k in 1:nk) {
      ava.k = which(!is.na(ctrl[, k]))
    }  #non-missing from kth control variable
    ava3 = intersect(ava3, ava.k)
  }
  if (nk == 1) 
    newctrl = ctrl[ava3]
  if (nk > 1) 
    newctrl = ctrl[ava3, ]
  newx = x[ava3]
  newy = y[ava3]
  list(newx = newx, newy = newy, newctrl = newctrl)  #delete NAs
}

