#' Function to do matched deletion of missing rows from x, y 
#' and z variable(s). 
#' 
#' The aim in three-way deletions is to retain only the largest 
#' number of available data triplets with all non-missing data.
#' This works where naTriplet fails (e.g.parcorVecH()). This is
#' called by parcorHijk
#' 
#' @param x Vector of x data
#' @param y Vector of y data
#' @param z vector or a matrix of additional variable(s) 
#' @return 
#' \item{newx}{A new vector x after removing triplet-wise missing data} 
#' \item{newy}{A new vector or matrix y after removing triplet-wise missing data} 
#' \item{newz}{A new vector or matrix ctrl after removing triplet-wise missing data} 
## @note %% ~~further notes~~
#' @importFrom stats complete.cases
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See \code{\link{napair}} \code{\link{naTriplet}}.
#' @examples
#' 
#' \dontrun{
#' x=sample(1:10);y=sample(1:10);x[2]=NA; y[3]=NA
#' w=sample(2:11)
#' naTriple(x,y,w)}
#' 
#' @export


naTriple=
function (x, y, z) 
{
  p = NCOL(y)
  pz = NCOL(z)
  ok = complete.cases(x, y, z)
    newx = x[ok]
    if (p == 1) 
      newy = y[ok]
    if (p > 1) 
      newy = y[ok, ]
    if (pz == 1) 
      newz = z[ok]
    if (pz > 1) 
      newz = z[ok, ]
  list(newx = newx, newy = newy, newz = newz)
}
