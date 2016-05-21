#'  Compute the numerical integration by the trapezoidal rule.
#' 
#' See page 220 of Vinod's ``Hands-on Interemediate Econometrics Using R,'' cited below
#' for the trapezoidal integration formula
#' needed for stochastic dominance.  The book explains pre-multiplication by two
#' large sparse matrices denoted by \eqn{I_F,  I_f}.  Here we accomplish the 
#' same computation without actually creating the large sparse matrices. For example, the
#' \eqn{I_f} is replaced by \code{cumsum} in this code (unlike the R code in
#' my textbook).
#' 
#' @param d {A vector of distances from the smallest value from both data vectors}
#' @param p {Vector of probabilities of the type 1/2T, 2/2T, 3/2T, etc. to 1.}
#' @return Returns a result after pre-multiplication by \eqn{I_F,  I_f}
#' matrices, without actually creating the large sparse matrices. This is an internal function.
#' @note This is an internal function, called by the function \code{stochdom2}, for
#'  comparison of two portfolios in terms of stochastic dominance (SD) of orders
#'  1 to 4.
#' Typical usage is:
#' \code{sd1b=bigfp(d=dj, p=rhs)
#' sd2b=bigfp(d=dj, p=sd1b)
#' sd3b=bigfp(d=dj, p=sd2b)
#' sd4b=bigfp(d=dj, p=sd3b)}.
#' This produces numerical evaluation vectors for the four orders, SD1  to SD4.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @references Vinod, H. D.', 'Hands-On Intermediate Econometrics
#' Using R'  (2008) World Scientific Publishers: Hackensack, NJ.
#' \url{http://www.worldscibooks.com/economics/6895.html}
#' @keywords SD1 SD2 SD3 SD4
#' 
#' @export

bigfp <- function(d, p) {
    n = length(d)
    # Fp1 vector represents the first term
    Fp1 = 0.5 * (d * p)
    dplusdp = rep(0, n)  #[di+d(i+1)] prefix of second term
    ans = rep(NA, n)
    for (i in 2:n) {
        dplusdp[i] = (d[i - 1] + d[i]) * p[i - 1]
    }
    term2 = 0.5 * cumsum(dplusdp)
    ans[1] = Fp1[1]
    for (i in 2:n) {
        ans[i] = Fp1[i] + term2[i]
    }
    return(ans)
} 
