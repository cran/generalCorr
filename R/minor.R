#' Function to do compute the minor of a matrix defined by row r and column c.
#' 
#' @param x {The input matrix}
#' @param r {The row number}
#' @param c {The column number}
#' @return The appropriate `minor' matrix defined from the input matrix.
#' 
#' @note This function is needed by the cofactor function.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @concept minor of a matrix
#' @examples
#' 
#' \dontrun{
#'  x=matrix(1:20,ncol=4)
#' minor(x,1,2)}
#' 
#' @export

minor <- function(x, r, c) {
    # x is n by p matrix we want its minor after eliminating rth row and cth column
    n = nrow(x)
    p = ncol(x)
    myn = 1:n
    myp = 1:p
    if (n < r) 
        stop("n<r, minor undefined")
    if (p < c) 
        stop("p<c, minor undefined")
    if (c <= 0) 
        stop("c<=0, minor undefined")
    newr = myn[-r]
    newc = myp[-c]
    out = x[newr, newc]
    return(out)
} 
