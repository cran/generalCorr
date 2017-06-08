#' Sort all columns of matrix x with respect to the j-th column.
#'
#' This function can use the sort.list function in R. The reason
#' for using it is that one wants the sort to carry along all columns.
#'
#' @param  x {An input matrix with several columns}
#' @param j {The column number with reference to which one wants to sort}
#' @return {A sorted matrix}
#' @examples
#'
#' set.seed(30)
#' x=matrix(sample(1:50),ncol=5)
#' y=sort_matrix(x,3);y
#' @export

sort_matrix <- function(x, j) 
{
  x[order(x[,j]),]
}
