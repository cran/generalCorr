#' Block version of generalized partial correlation coefficients between Xi 
#' and Xj, after removing the
#' effect of xk, via nonparametric regression residuals.
#'
#' This function uses data on two column vectors, xi, xj and a third
#' xk which can be a vector or a matrix, usually of the remaining 
#' variables in the model, including control variables, if any.
#' It first removes missing data from all input variables. Then,
#' it computes residuals of kernel regression (xi on xk) and (xj on xk). 
#' This is a block version of parcor_ijk.
#'
#' @param xi {Input vector of data for variable xi}
#' @param xj {Input vector of data for variable xj}
#' @param xk {Input data for variables in xk, usually control variables}
#' @param blksiz {block size, default=10, if chosen blksiz >n, where n=rows in matrix
#'      then blksiz=n. That is, no blocking is done}
#' @importFrom stats cov
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
#' parcorBijk(x[,1], x[,2], x[,3], blksi=10)
#' }#' 
#' @export

parcorBijk =function (xi, xj, xk, blksiz=10) 
{
  na2 = naTriplet(x=xi, y=xj, ctrl=xk)
  xi = na2$newx
  xj = na2$newy
  xk = na2$newctrl
  n=NROW(xi)
  q=NCOL(xk)
  if (blksiz>n) blksiz=n
  ge=getSeq(n,blksiz=blksiz)
  uuik=rep(NA,n)
  uujk=rep(NA,n)
  LO=ge$sqLO
  UP=ge$sqUP
  #   print(cbind(LO,UP))
  k=length(LO)
  
  for (ik in 1:k){
    L1=LO[ik]  
    U1=UP[ik]
 
    xxi= xi[L1:U1] 
    xxj= xj[L1:U1]
    if (q==1) xxk=xk[L1:U1]
    if (q>1) xxk=xk[L1:U1,]
#print(length(L1:U1))

  uik = kern(dep.y = xxi, reg.x = xxk, residuals = TRUE)$resid
  ujk = kern(dep.y = xxj, reg.x = xxk, residuals = TRUE)$resid
#print(length(uik))
#print(length(ujk))
  uuik[L1:U1] = uik
  uujk[L1:U1] = ujk
  } # end ik loop
  sgn = sign(cov(uik, ujk))
  ouij = sgn * kern(dep.y = uuik, reg.x = uujk)$R2
  ouji = sgn * kern(dep.y = uujk, reg.x = uuik)$R2
  list(ouij = ouij, ouji = ouji)
}
