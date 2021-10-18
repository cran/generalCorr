#' Vector of generalized partial correlation coefficients (GPCC), 
#' always leaving out control variables, if any.
#'
#' This function calls  \code{parcor_ijk} function which
#' uses original data to compute
#' generalized partial correlations between \eqn{X_i}, the dependent variable,
#' and \eqn{X_j} which is the current regressor of interest. Note that
#' j can be any one of the remaining
#' variables in the input matrix \code{mtx}. Partial correlations remove the effect of
#' variables \eqn{X_k} other than \eqn{X_i} and \eqn{X_j}. 
#' Calculation merges control variable(s) (if any) into \eqn{X_k}.
#' Let the remainder effect
#' from kernel regressions of \eqn{X_i} on \eqn{X_k} equal the  
#' residuals u*(i,k). Analogously define u*(j,k). (asterisk for kernel regressions) 
#' Now partial correlation is generalized correlation
#' between  u*(i,k) and u*(j,k).
#' Calculation merges control variable(s) (if any) into \eqn{X_k}.
#' 
#'
#' @param mtx {Input data matrix with p (> or = 3) columns} 
#' @param ctrl {Input vector or matrix of data for control variable(s), 
#'     default is ctrl=0 when control variables are absent}
#' @param verbo Make this TRUE for detailed printing of computational steps
#' @param idep The column number of the dependent variable (=1, default)
#' @return A p by 1 `out' vector containing partials  r*(i,j | k).
#'
#' @note Generalized Partial Correlation Coefficients (GPCC) allow comparison of
#' the relative contribution of each \eqn{X_j} to the explanation of \eqn{X_i},
#' because GPCC are scale-free pure numbers 
#'
#' @note We want to get all partial
#'  correlation coefficient pairs removing other column effects. Vinod (2018) 
#'  shows why one needs more than one criterion to decide the causal paths or exogeneity.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See Also \code{\link{parcor_ijk}}.
#' @seealso See Also a hybrid version \code{\link{parcorVecH}}.
#' @concept  partial correlations
#' @references Vinod, H. D. 'Generalized Correlations and Instantaneous
#'  Causality for Data Pairs Benchmark,' (March 8, 2015)
#'  \url{https://www.ssrn.com/abstract=2574891}
#'
#' @references Vinod, H. D. 'Matrix Algebra Topics in Statistics and Economics
#'  Using R', Chapter 4 in Handbook of Statistics: Computational Statistics
#'  with R, Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#'  North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#'  
#' @references Vinod, H. D. 'New Exogeneity Tests and Causal Paths,' 
#' (June 30, 2018). Available at SSRN: 
#' \url{https://www.ssrn.com/abstract=3206096}    
#' 
#' @references Vinod, H. D. (2021) 'Generalized, Partial and Canonical Correlation
#' Coefficients' Computational Economics, 59(1), 1--28.
#'   
#' @examples
#' set.seed(234)
#' z=runif(10,2,11)# z is independently created
#' x=sample(1:10)+z/10  #x is partly indep and partly affected by z
#' y=1+2*x+3*z+rnorm(10)# y depends on x and z not vice versa
#' mtx=cbind(x,y,z)
#' parcorVec(mtx)
#'  
#'    
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c('V1', 'v2', 'V3')#some names needed
#' parcorVec(x)
#' }
#'
#' @export


parcorVec=
function (mtx, ctrl = 0, verbo = FALSE, idep=1) 
{
  n = NROW(mtx)
  p = NCOL(mtx)
  if (p < 3) 
    stop("stop: input matrix to parcorVec must have 3 or more columns")
  nam = colnames(mtx)
  if (length(nam) == 0) 
    nam = paste("V", 1:p, sep = "")
  if (verbo) 
    print(c("gen.zed partial Corr coef GPCC vec", nam[idep], "and one of the others"))
  if (verbo) 
    print(c("Removes effect of control variables in ctrl (if any)"))
  out = matrix(1, nrow = 1, ncol = (p-1))#single row output one less column
  namj=rep(NA,p-1)
  j.other = setdiff(1:p, idep)
  namj=nam[j.other]
  p.other = length(j.other)
  if (verbo) 
    print(c("j.other,p.other",j.other,p.other))
  xi=mtx[,idep]
  for (j in  1:p.other) {
    myj=j.other[j]
    xj = mtx[,myj]
    xk = mtx[, c(-idep, -myj)]
      if (length(ctrl) > 1) {
   p1 = parcor_ijk(xi = xi, xj = xj, xk = cbind(xk,  ctrl))
      }
      if (length(ctrl) == 1) {
        p1 = parcor_ijk(xi = xi, xj = xj, xk = xk)
      }
      out[1, j] = p1$ouij
  }
  colnames(out) = namj
  return(out)
}

