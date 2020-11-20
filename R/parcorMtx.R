#' Matrix of generalized partial correlation coefficients, 
#' always leaving out control variables, if any.
#'
#' This function calls  \code{parcor_ijk} function which
#' uses original data to compute
#' generalized partial correlations between \eqn{X_i} and \eqn{X_j}
#' where j can be any one of the remaining
#' variables in the input matrix \code{mtx}. Partial correlations remove the effect of
#' variables \eqn{x_k} other than \eqn{X_i} and \eqn{X_j}. Calculation further 
#' allows for the presence of control variable(s) (if any) to remain always outside
#' the input matrix and whose effect is also removed in computing partial correlations.
#' 
#'
#' @param mtx {Input data matrix with p columns. p is at least 3 columns.}
#' @param ctrl {Input vector or matrix of data for control variable(s), 
#'     default is ctrl=0 when control variables are absent}
#' @param dig The number of digits for reporting (=4, default)
#' @param verbo Make this TRUE for detailed printing of computational steps
#' @return A p by p `out' matrix containing partials  r*(i,j | k).
#'   and  r*(j,i | k).
#'
#' @note We want to get all partial
#'  correlation coefficient pairs removing other column effects. Vinod (2018) 
#'  shows why one needs more than one criterion to decide the causal paths or exogeneity.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See Also \code{\link{parcor_ijk}}.
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
#' @examples
#' set.seed(234)
#' z=runif(10,2,11)# z is independently created
#' x=sample(1:10)+z/10  #x is partly indep and partly affected by z
#' y=1+2*x+3*z+rnorm(10)# y depends on x and z not vice versa
#' mtx=cbind(x,y,z)
#' parcorMtx(mtx)
#'  
#'    
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c('V1', 'v2', 'V3')
#' parcorMtx(x)
#' }
#'
#' @export

parcorMtx <- function(mtx, ctrl=0, dig = 4, verbo = FALSE) {
    n = NROW(mtx)
    p = NCOL(mtx)
    if (p<3) stop("input matrix to parcorMtx must have 3 or more columns")
   nam = colnames(mtx)  #R makes nam=NULL of lenghth 0 if mtx column names Missing
    if (length(nam) == 0) 
        nam = paste("V", 1:p, sep = "")
  if(verbo) print(c("We want partial Corr Matrix"))
  if(verbo) print(c("The effect of control variables in ctrl is always removed"))
        out = matrix(1, nrow = p, ncol = p) #diag=1
        for (i in 1:(p-1)) {
        for (j in (i+1):(p)){
        xi=mtx[,i]
        xj=mtx[,j]
        xk=mtx[,c(-i,-j)]
            if (length(ctrl)>1){
            p1 = parcor_ijk(xi=xi, xj=xj, xk=cbind(xk,ctrl))}
            if (length(ctrl)==1){
              p1 = parcor_ijk(xi=xi, xj=xj, xk=xk)}
        out[i,j]= p1$ouij
        out[j,i]= p1$ouji
        }#end j loop
        }  #end i loop
    rownames(out) = nam
    colnames(out) =nam
    return(out)
}
