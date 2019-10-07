#' Compute vectors measuring stochastic dominance of four orders.
#' 
#' Stochastic dominance originated as a sophisticated comparison of two distributions of
#' stock market returns.  The dominating distribution is superior in terms of local
#' mean, variance, skewness and kurtosis respectively, representing dominance
#' orders 1 to 4, without simply computing the four moment summary measures for the entire
#' data.  Vinod (2008, sec. 4.3)
#' explains the details.  This function uses the output of `wtdpapb'.
#' 
#' 
#' @param dj {Vector of (unequal) distances of consecutive intervals defined on common support
#'    of two probability distributions being compared}
#' @param wpa Vector of the first set of (weighted) probabilities
#' @param wpb Vector of the second set of (weighted) probabilities
#' @return 
#' \item{sd1b}{Vector measuring stochastic dominance of order 1, SD1} 
#' \item{sd2b}{Vector measuring stochastic dominance of order 2, SD2} 
#' \item{sd3b}{Vector measuring stochastic dominance of order 3, SD3} 
#' \item{sd4b}{Vector measuring stochastic dominance of order 4, SD4} 
#' @note The input to this function is the output of the function \code{wtdpapb}.
## @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{wtdpapb}}
#' @references Vinod, H. D.', 'Hands-On Intermediate Econometrics 
#'  Using R'  (2008) World Scientific Publishers: Hackensack, NJ.
#'  \url{https://www.worldscientific.com/worldscibooks/10.1142/6895}
#'  
#'
#' @references Vinod, H. D. 'Ranking Mutual Funds Using 
#' Unconventional Utility Theory and Stochastic Dominance,'
#' Journal of Empirical Finance Vol. 11(3) 2004, pp. 353-377.
#' 
#' @concept stochastic dominance from local skewness 
#' @concept stochastic dominance from local kurtosis
#' @examples
#' 
#'  \dontrun{
#'  set.seed(234);x=sample(1:30);y=sample(5:34)
#'  w1=wtdpapb(x,y) #y should dominate x with mostly positive SDs
#'  stochdom2(w1$dj, w1$wpa, w1$wpb) }
#' 
#' @export

stochdom2 <- function(dj, wpa, wpb) {
    # input weighted pa and pb IsubF and I sub f matrices
    if (length(wpa) != length(dj)) 
        print("wrong rows dj")
    if (length(wpa) != length(wpb)) 
        print("wrong rows wpa wpb")
    rhs = cumsum(wpa - wpb)
    sd1b = bigfp(d = dj, p = rhs)
    # sd1=I.bigf %*% I.smallf %*% (wpa-wpb)
    sd2b = bigfp(d = dj, p = sd1b)
    sd3b = bigfp(d = dj, p = sd2b)
    sd4b = bigfp(d = dj, p = sd3b)
    list(sd1b = sd1b, sd2b = sd2b, sd3b = sd3b, sd4b = sd4b)
} 
