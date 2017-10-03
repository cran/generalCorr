#' Partial correlation coefficient between Xi and Xj after removing the linear
#' effect of all others.
#'
#' This function uses a symmetric correlation matrix R as input to compute
#' usual partial correlations between \eqn{X_i} and \eqn{X_j}
#' where j can be any one of the remaining
#' variables. Computation removes the effect of all other variables in the matrix.
#' The user is encouraged to remove all known irrelevant rows and columns 
#' from the R matrix before submitting it to this function.
#'
#' @param x {Input a p by p matrix R of symmetric correlation coefficients.}
#' @param i {A column number identifying the first variable.}
#' @param j {A column number identifying the second variable.}
#' @return 
#' \item{ouij}{Partial correlation Xi with Xj after removing all other X's}
#' \item{ouji}{Partial correlation Xj with Xi after removing all other X's}
#' \item{myk}{A list of column numbers whose effect has been removed}
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{parcor_ijk}} for generalized partial
#'  correlation coefficients useful for causal path determinations.
#' @note This function calls \code{\link{minor}}, and \code{\link{cofactor}} 
#' @examples 
#' 
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c('V1', 'v2', 'V3')
#' c1=cor(x)
#' parcor_linear(c1, 2,3)
#' }
#' 
#' @export

parcor_linear = function(x, i, j) {
    n = NROW(x)
    p = NCOL(x)
    if (n < i) 
        stop("n<i, parcor undefined")
    if (p < j) 
        stop("p<j, parcor undefined")
    if (i <= 0 | j <= 0) 
        stop("i OR j <=0, parcor undefined")
    myn = 1:n
    myp = 1:p
    myk = myp[c(-i, -j)]
    numij = det(cofactor(x, i, j))
    numji = det(cofactor(x, j, i))
    deni = abs(det(cofactor(x, i, i)))
    denj = abs(det(cofactor(x, j, j)))
    ouij = (numij)/sqrt(deni * denj)
    ouji = (numji)/sqrt(deni * denj)
    list(ouij = ouij, ouji = ouji, myk = myk)
}
