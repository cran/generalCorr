#' Compute cofactor of a matrix based on row r and column c.
#' 
#' 
#' @param x {matrix whose cofactor is desired to be computed}
#' @param r {row number}
#' @param c {column number}
#' @return cofactor of x,  w.r.t. row r and column c.
#' @note needs the function `minor'' in memory. attaches sign (-1)^(r+c) to the minor.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso \code{minor(x,r,c)}
#' @concept cofactor of a matrix
#' @examples
#' 
#' ## The function is currently defined as
#' function (x, r, c) 
#' {
#'     out = minor(x, r, c) * ((-1)^(r + c))
#'     return(out)
#'   }
#' @export

cofactor <- function(x, r, c) {
    out = minor(x, r, c) * ((-1)^(r + c))
    return(out)
}
