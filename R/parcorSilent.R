#' Silently compute generalized (ridge-adjusted) partial correlation coefficients from matrix R*.
#'
#' This function calls \code{parcor_ijkOLD} function which
#' uses a generalized correlation matrix R* as input to compute
#' generalized partial correlations between \eqn{X_i} and \eqn{X_j}
#' where j can be any one of the remaining
#' variables. Computation removes the effect of all other variables in the matrix.
#' It further adjusts the resulting partial correlation coefficients to be in the
#' appropriate [-1,1] range by using an additive constant in the fashion of ridge regression.
#' 
#'
#' @param gmc0 This must be a p by p matrix R* of generalized correlation coefficients.
#' @param dig The number of digits for reporting (=4, default)
#' @param idep The column number of the first variable (=1, default)
#' @param verbo Make this TRUE for detailed printing of computational steps
#' @param incr {incremental constant for iteratively adjusting `ridgek'
#'        where ridgek is the constant times the identity matrix used to
#'    make sure that the gmc0 matrix is positive definite. If not, this function iteratively
#'  increases the \code{incr} till relevant partial correlations are within the [-1,1] interval.}
#' @return A five column `out' matrix containing partials. The first column
#'   has the name of the \code{idep} variable. The
#'    second column has the name of the j variable, while the third column has  r*(i,j | k).
#'   The 4-th column has  r*(j,i | k) (denoted partji), and the 5-th column has rijMrji,
#'   that is the difference in absolute values (abs(partij) - abs(partji)).
#'
#' @note The ridgek constant created by the function during the first round
#'  may not be large enough to make sure that
#'  that other pairs of r*(i,j | k) are within the [-1,1] interval. The user may have to choose
#'  a suitably larger input \code{incr} to get all relevant partial
#'  correlation coefficients in  the correct [-1,1] interval. 
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See Also \code{\link{parcor_ijk}} for a better version using original data as input.
#' @keywords partial correlations, ridge biasing factor,
#' @references Vinod, H. D. 'Generalized Correlations and Instantaneous
#'  Causality for Data Pairs Benchmark,' (March 8, 2015)
#'  \url{http://ssrn.com/abstract=2574891}
#'
#' @references Vinod, H. D. 'Matrix Algebra Topics in Statistics and Economics
#'  Using R', Chapter 4 in Handbook of Statistics: Computational Statistics
#'  with R, Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#'  North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#' @references Vinod, H. D. "A Survey of Ridge Regression and Related Techniques 
#' for Improvements over Ordinary Least Squares," Review of Economics and Statistics, 
#' Vol. 60, February 1978, pp. 121-131.
#' @examples
#' set.seed(234)
#' z=runif(10,2,11)# z is independently created
#' x=sample(1:10)+z/10  #x is partly indep and partly affected by z
#' y=1+2*x+3*z+rnorm(10)# y depends on x and z not vice versa
#' mtx=cbind(x,y,z)
#' g1=gmcmtx0(mtx)
#' parcor_ijkOLD(g1,1,2) # ouji> ouij implies i=x is the cause of j=y
#' parcor_ridg(g1,idep=1)
#' parcorSilent(g1,idep=1)
#'  
#'    
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c('V1', 'v2', 'V3')
#' gm1=gmcmtx0(x)
#' parcorSilent(gm1, idep=1)
#' }
#'
#' @export

parcorSilent <- function(gmc0, dig = 4, idep = 1, verbo = FALSE, incr = 3) {
    n = NROW(gmc0)
    p = NCOL(gmc0)
    if (n != p) 
        stop("gmc0 is not a square matrix")
#    nam = colnames(gmc0)  #R makes nam=NULL of lenghth 0 if gmc0 column names Missing
#    if (length(nam) == 0) 
#        nam = paste("V", 1:p, sep = "")
#    print(c("We want partial Corr of", nam[idep], "w.r.t. others"))
    j.other = setdiff(1:n, idep)
    e0 = eigen(gmc0)
    sort.e0 = sort(e0$values)
    sort.abse0 = sort(abs(e0$values))
    min.e0 = sort.e0[1]
    min.e0
    diff.e0 = incr * (sort.abse0[2] - sort.abse0[1])

    sgn.e0 = sign(Re(min.e0))
    sgn.e0
    ridgek = 0  #initialize for positive definite case
    if (sgn.e0 == -1) 
        ridgek = sgn.e0 * round(Re(min.e0) + diff.e0, 3)
    for (jridge in 1:10) {
#        print(c("ridgek=", ridgek))
        
        if (jridge >= 2) 
            ridgek = ridgek + incr * diff.e0
        gmc1 = gmc0 + ridgek * diag(ncol(gmc0))
        nam = colnames(gmc0)
        if (jridge == 1) 
            {
 #               print(c("We want partial Corr of", nam[idep], "w.r.t. others"))
            }  #eindif printing
        n = length(nam)
#        print(c("nam length=", n))
        dig = 4  #digits of accuracy desired
        out1 = matrix(NA, nrow = n - 1, ncol = 4)
        partij = rep(NA, n - 1)  #place holders
        partji = rep(NA, n - 1)
        ii = 0
        for (i in 2:n) {
            p1 = parcor_ijkOLD(gmc1, 1, i)
            ii = ii + 1
            partij[ii] = p1$ouij
            partji[ii] = p1$ouji
        }  #end i lop
        if ((max(abs(partij)) < 1) & (max(abs(partji)) < 1)) 
            break
    }  #end jridge loop
    rijMrji = (abs(partij) - abs(partji))
#    print(c("final ridgek=", ridgek))
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
