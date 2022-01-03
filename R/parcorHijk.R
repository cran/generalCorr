#' Generalized partial correlation coefficients between Xi and Xj, after removing the
#' effect of Xk, via OLS regression residuals.
#'
#' This function uses data on two column vectors, xi, xj, and a third set
#' xk, which can be a vector or a matrix. xk usually has the remaining 
#' variables in the model, including control variables, if any. This function
#' first removes missing data from all input variables. Then,
#' it computes residuals of OLS (no kernel) regression (xi on xk) and (xj on xk). 
#' This hybrid version uses both OLS and then generalized correlation among
#' OLS residuals. This solves the potential problem of having too little
#' information content in kernel regression residuals, since kernel fits are
#' sometimes too close, especially when there are many variables in xk.
#'
#' @param xi {Input vector of data for variable xi}
#' @param xj {Input vector of data for variable xj}
#' @param xk {Input data for all variables in xk, usually control variables}
#' @importFrom stats cov lm resid
#' @return 
#' \item{ouij}{Generalized partial correlation Xi with Xj (=cause) after removing xk}
#' \item{ouji}{Generalized partial correlation Xj with Xi (=cause) after removing xk}
#' allowing for control variables. 
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{parcor_ijk}}.
#' @note This function calls \code{\link{kern}}, 
#' @examples 
#' 
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' options(np.messages=FALSE)
#' parcorHijk(x[,1], x[,2], x[,3])
#' }#' 
#' @export


parcorHijk=function (xi, xj, xk) 
{
  na2 = naTriplet(x = xi, y = xj, ctrl = xk)
  xi = na2$newx
  xj = na2$newy
  xk = na2$newctrl
  uik = resid(lm(xi~xk))
  ujk = resid(lm(xj~xk))
  sgn = sign(cov(uik, ujk))
  ouij = sgn * kern(dep.y = uik, reg.x = ujk)$R2
  ouji = sgn * kern(dep.y = ujk, reg.x = uik)$R2
  list(ouij = ouij, ouji = ouji)
}
