#' Generalized canonical correlation, estimating alpha, beta, rho.
#'
#' What exactly is generalized? Canonical correlations start with Rij, a
#' symmetric matrix of Pearson correlation coefficients based on linear
#' relations. This function starts with a more general non-symmetric R*ij
#' produced by \code{gmcmtx0()} as an input. This is a superior measure
#' of dependence, allowing for nonlinear dependencies. It generalizes Hotelling's
#' derivation for the nonlinear case.
#' This function uses data on two sets of column vectors. LHS set [x1, x2 .. xr]
#' has r=nLHS number of columns 
#' with coefficients alpha, and 
#' the larger RHS set [xr+1, xr+2, .. xp] has nRHS=(p-r) columns and RHS
#' coefficients beta.  Must arrange the sets so that the larger set
#' in on RHS with coefficients beta estimated first from an eigenvector
#' of the problem [A* beta = rho^2 beta], where A* is a partitioning of our
#' generalized matrix of (non-symmetric) correlation coefficients.
#'
#' @param mtx {Input matrix of generalized correlation coefficients R*}
#' @param nLHS {number of columns in the LHS set, default=2}
#' @param sgn {preferred sign of coefficients default=1 for positive, 
#'    use sgn= -1 if prior knowledge suggests that negative signs of 
#'    coefficients are more realistic}
#' @param verbo {logical, verbo=FALSE default means do not print results}
#' @param ridg {two regularization constants added 
#'  before computing matrix inverses of S11 and S22, respectively, with
#'      default=c(0,0). Some suggest ridg=c(0.01,0.01) for stable results}
#' @references Vinod, H. D. 'Matrix Algebra Topics in Statistics and Economics
#' Using R', Chapter 4 in 'Handbook of Statistics: Computational Statistics
#' with R', Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#' @references Vinod, H. D. 'Canonical ridge and econometrics of joint production,'
#' Journal of Econometrics, vol. 4, 147--166.
#' @references Vinod, H. D. 'New exogeneity tests and causal paths,'
#'  Chapter 2 in 'Handbook of Statistics: Conceptual Econometrics 
#' Using R', Vol.32, co-editors: H. D. Vinod and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2019, pp. 33-64.
#'        
#' @return 
#' \item{A}{eigenvalue computing matrix for Generalized canonical correlations}
#' \item{rho}{Generalized canonical correlation coefficient}
#' \item{bet}{RHS coefficient vector}
#' \item{alp}{LHS coefficient vector}
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{gmcmtx0}}.
#' @note This function calls \code{\link{kern}}, 
#' @examples 
#' 
#' \dontrun{
#' set.seed(99)
#' mtx2=matrix(sample(1:25),nrow=5)
#' g1=gmcmtx0(mtx2)
#' canonRho(g1,verbo=TRUE)
#' }#' 
#' @export

canonRho=function(mtx, nLHS=2, sgn=1, verbo=FALSE, ridg=c(0,0)){
p=NCOL(mtx)
nRHS=p-nLHS
if(verbo){
print("p=no.of columns in mtx")
print("nLHS=no.of columns selected in the smaller set r=min(r, p-r)")
print("nRHS=no.of columns selected in the larger set")
print(c("p,nLHS,nRHS",p,nLHS,nRHS))  
}
if (p!=(nLHS+nRHS)) stop("wrong nLHS, nRHS or p",nLHS,nRHS,p)
#if (nLHS>nRHS) stop("nLHS exceeds nRHS")
r=nLHS
rp1=r+1


S11=mtx[1:r,1:r]
S12=mtx[1:r,rp1:p]
S21=mtx[rp1:p,1:r]
S22=mtx[rp1:p,rp1:p]


I1=diag(r)
I2=diag(nRHS)

##A=solve(S22)%*%(S21)%*%(solve(S11))%*%(S12)
SS11=S11+t(S11)+ridg[1]*I1
SS22=S22+t(S22)+ridg[2]*I2
SS21=S21+t(S12)
SS12=S12+t(S21)
A=(solve(SS22))%*%(SS21)%*%(solve(SS11)%*%(SS12))


if (verbo) print(A)
ei=eigen(A)
bet=ei$vector[,1]#asso.with largest e-valu
sbet=sum(sign(bet))
if(verbo) print(c("RHS coef=bet=",bet))
if(sbet!=0){if (sbet!=sgn) bet=-bet}
if(verbo) print(c("sign-revised RHS coef=bet=",bet))

rho=sqrt(ei$value[1])
if(verbo) { print(c("rho=",rho))
print(A%*%bet)
print("above A * beta, rho^2* beta below")
print((rho^2)*bet)
print(SS21)
print("SS21 above SS12 below")
print(SS12)
print(solve(SS11))
print("SS11 inverse above and bet below")
print(bet)
}
alp=(solve(SS11))%*%SS12 %*% bet/rho  
#alp=t((solve(SS11))%*%t(SS21)) %*%t(t(bet/rho)) #says Fred


salp=sum(sign(alp))
if(salp!=0){if (salp!=sgn) alp=-alp}

if(verbo) print(c("LHS coef=alp=",alp))
list(A=A, rho=rho, bet=bet, alp=alp)
} #end function canonRho
