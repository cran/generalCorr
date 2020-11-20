#' Block version reports many generalized partial correlation coefficients 
#' allowing control variables.
#'
#' This function calls a block version \code{parcorBijk} of the function which
#' uses original data to compute
#' generalized partial correlations between \eqn{X_{idep}} and \eqn{X_j}
#' where j can be any one of the remaining
#' variables in the input matrix \code{mtx}. Partial correlations remove the effect of
#' variables \eqn{X_k} other than \eqn{X_i} and \eqn{X_j}. Calculation further 
#' allows for the presence of control variable(s) (if any) to remain always outside
#' the input matrix and whose effect is also removed in computing partial correlations.
#' 
#'
#' @param mtx {Input data matrix with at least 3 columns.}
#' @param ctrl {Input vector or matrix of data for control variable(s), 
#'     default is ctrl=0 when control variables are absent}
#' @param dig The number of digits for reporting (=4, default)
#' @param idep The column number of the first variable (=1, default)
#' @param verbo Make this TRUE for detailed printing of computational steps
#' @param blksiz {block size, default=10, if chosen blksiz >n, where n=rows in matrix
#'      then blksiz=n. That is, no blocking is done}
#' @return A five column `out' matrix containing partials. The first column
#'   has the name of the \code{idep} variable. The
#'    second column has the name of the j variable, while the third column 
#'    has partial correlation coefficients  r*(i,j | k).The last column
#'    reports the absolute difference between two partial correlations.
#'
#' @note This function reports all partial
#'  correlation coefficients, while avoiding ridge type adjustment.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See Also \code{\link{parcor_ijk}}, \code{\link{parcorMany}}.
#' @concept  partial correlations
#' @references Vinod, H. D. 'Generalized Correlations and Instantaneous
#'  Causality for Data Pairs Benchmark,' (March 8, 2015)
#'  \url{https://www.ssrn.com/abstract=2574891}
#'
#' @references Vinod, H. D. 'Matrix Algebra Topics in Statistics and Economics
#'  Using R', Chapter 4 in Handbook of Statistics: Computational Statistics
#'  with R, Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#'  North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#' @examples
#' set.seed(234)
#' z=runif(10,2,11)# z is independently created
#' x=sample(1:10)+z/10  #x is partly indep and partly affected by z
#' y=1+2*x+3*z+rnorm(10)# y depends on x and z not vice versa
#' mtx=cbind(x,y,z)
#' parcorBMany(mtx, blksiz=10)
#'  
#'    
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c('V1', 'v2', 'V3')
#' parcorBMany(x, idep=1)
#' }
#'
#' @export

parcorBMany <- function(mtx, ctrl=0, dig = 4, idep = 1, 
                       blksiz=10, verbo = FALSE) {
    n = NROW(mtx)
    p = NCOL(mtx)
    if (p<3) stop("input matrix to parcorBMany must have 3 or more columns")
   nam = colnames(mtx)  #R makes nam=NULL of lenghth 0 if mtx column names Missing
    if (length(nam) == 0) 
        nam = paste("V", 1:p, sep = "")
  if(verbo) print(c("We want partial Corr of", nam[idep], "w.r.t. others"))
    j.other = setdiff(1:p, idep)
    n.other = length(j.other)
#        p2 = length(nam)
#        print(c("nam length=", n))
        out1 = matrix(NA, nrow = p - 1, ncol = 4)
        partij = rep(NA, p - 1)  #place holders
        partji = rep(NA, p - 1)
        ii = 0
        for (i in 1:n.other) {
           myi=j.other[i]
            xk=mtx[,c(-idep,-myi)]
            if (length(ctrl)>1){
  p1 = parcorBijk(xi=mtx[,idep], xj=mtx[,myi], xk=cbind(xk,ctrl, blksiz=blksiz))}
            if (length(ctrl)==1){
  p1 = parcorBijk(xi=mtx[,idep], xj=mtx[,myi], xk=xk, blksiz=blksiz)}
            ii = ii + 1
            partij[ii] = p1$ouij
            partji[ii] = p1$ouji
        }  #end i loop
     rijMrji = (abs(partij) - abs(partji))
    cb1 = cbind(partij, partji, rijMrji)
    cb2 = apply(cb1, 2, round, dig)
    if (verbo) 
        print(cb2)
    m = length(partij)
    nami = rep(nam[idep], m)
    namj = nam[j.other]
    out = cbind(nami, namj, cb2)
    return(out)
}
