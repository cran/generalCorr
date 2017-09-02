#' Generalized partial correlation coefficient between Xi and Xj after removing the
#' effect of all others.
#'
#' This function uses a generalized correlation matrix R* as input to compute
#' generalized partial correlations between \eqn{X_i} and \eqn{X_j}
#' where j can be any one of the remaining
#' variables. Computation removes the effect of all other variables in the matrix.
#' The user is encouraged to remove all known irrelevant rows and columns 
#' from the R* matrix before submitting it to this function.
#'
#' @param x {Input a p by p matrix R* of generalized correlation coefficients.}
#' @param i {A column number identifying the first variable.}
#' @param j {A column number identifying the second variable.}
#' @return 
#' \item{ouij}{Partial correlation Xi with Xj (=cause) after removing all other X's}
#' \item{ouji}{Partial correlation Xj with Xi (=cause) after removing all other X's}
#' \item{myk}{A list of column numbers whose effect has been removed}
#' @note This function calls \code{\link{minor}}, and \code{\link{cofactor}} and is called 
#'   by \code{parcor_ridge}.
#' @examples 
#' 
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c('V1', 'v2', 'V3')
#' gm1=gmcmtx0(x)
#' parcor_ijk(gm1, 2,3)
#' }#' 
#' @export

parcor_ijk = function(x, i, j) {
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
    numij = -det(cofactor(x, i, j))
    numji = -det(cofactor(x, j, i))
    deni = abs(det(cofactor(x, i, i)))
    denj = abs(det(cofactor(x, j, j)))
    ouij = (numij)/sqrt(deni * denj)
    ouji = (numji)/sqrt(deni * denj)
    list(ouij = ouij, ouji = ouji, myk = myk)
}
