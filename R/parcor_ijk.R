#' Generalized partial correlation coefficients between Xi and Xj after removing the
#' effect of xk via nonparametric regression residuals.
#'
#' This function uses data on two column vectors, xi, xj and 
#' xk which can be a vector or a matrix usually of the remaining 
#' variables in the model including optional control
#' variables. It works with kernel regression (xi on xk) and (xj on xk)  residuals,
#' removes missing data from input variables before proceeding.
#'
#' @param xi {Input vector of data for variable xi}
#' @param xj {Input vector of data for variable xj}
#' @param xk {Input data for variables in xk, usually control variables}
#' @return 
#' \item{ouij}{Generalized partial correlation Xi with Xj (=cause) after removing xk}
#' \item{ouji}{Generalized partial correlation Xj with Xi (=cause) after removing xk}
#' allowing for control variables. 
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{parcor_linear}}.
#' @note This function calls \code{\link{kern}}, 
#' @examples 
#' 
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c('V1', 'v2', 'V3')
#' parcor_ijk(gm1, 2,3)
#' }#' 
#' @export

parcor_ijk =function (xi, xj, xk) 
{
  na2 = naTriplet(x=xi, y=xj, ctrl=xk)
  xi = na2$newx
  xj = na2$newy
  xk = na2$newctrl
  uik = kern(dep.y = xi, reg.x = xk, residuals = TRUE)$resid
  ujk = kern(dep.y = xj, reg.x = xk, residuals = TRUE)$resid
  sgn = sign(cor(uik, ujk))
  ouij = sgn * kern(dep.y = uik, reg.x = ujk)$R2
  ouji = sgn * kern(dep.y = ujk, reg.x = uik)$R2
  list(ouij = ouij, ouji = ouji)
}
