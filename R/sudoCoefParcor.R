#' Peudo regression coefficients from generalized partial correlation coefficients, 
#' (GPCC).
#'
#' This function gets GPCCs by calling  \code{parcorVec} function.
#' Pseudo regression coefficient of a kernel regression is obtained by
#' GPCC*(sd dep.var)/(sd regressor), that is
#' multiplying the GPCC by
#' the standard deviation (sd) of the dependent variable and dividing by the
#' sd of the regressor.  
#' 
#'
#' @param mtx {Input data matrix with p (> or = 3) columns},
#' @param ctrl {Input vector or matrix of data for control variable(s), 
#'     default is ctrl=0 when control variables are absent}
#' @param verbo Make this TRUE for detailed printing of computational steps
#' @param idep The column number of the dependent variable (=1, default)
#' @return A p by 1 `out' vector pseudo partial derivatives 
#'
#' @note Generalized Partial Correlation Coefficients (GPCC) allow comparison of
#' the relative contribution of each \eqn{X_j} to the explanation of \eqn{X_i},
#' because GPCC are scale-free pseudo regr coeff are GPCC*(sd dep.var)/(sd regressor)
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
#' sudoCoefParcor(mtx, idep=2)
#'  
#'    
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c('V1', 'v2', 'V3')#some names needed
#' sudoCoefParcor(x)
#' }
#'
#' @export


sudoCoefParcor=
function (mtx, ctrl = 0, verbo = FALSE, idep=1) 
{
mysd=apply(mtx,2,sd,na.rm=TRUE)
if (verbo) print(c("std dev",mysd))
p = NCOL(mtx)
if (p < 3) 
  stop("stop: input matrix to sudoCoefParcor must have 3 or more columns")
gpcc=parcorVec(mtx=mtx,ctrl=ctrl, verbo=verbo, idep=idep)
if (verbo) print(gpcc)
nam1=colnames(gpcc)
if(verbo) print(nam1)
pp=length(gpcc)
if(verbo) print(pp)
out=rep(NA,pp)
for (i in 1:pp){
out[i]=gpcc[i]*mysd[idep]/sd(mtx[,nam1[i]],na.rm = TRUE)
}#end for loop
  names(out) = nam1
  return(out)
}
