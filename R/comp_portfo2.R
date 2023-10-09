#' Compares two vectors (portfolios) using stochastic dominance of orders 1 to 4.
#'
#' Given two vectors of portfolio returns this function calls the internal function wtdpapb
#' to report the simple means of four sophisticated measures of stochastic dominance.
#' as explained in Vinod (2008).
#'
#' @param xa {Data on returns for portfolio A in the form of a T by 1 vector}
#' @param xb {Data on returns for portfolio B in the form of a T by 1 vector}
#' @return Returns four numbers which are averages of four sophisticated measures of stochastic
#' dominance measurements called SD1 to SD4.
#' @note It is possible to modify this function to report the median or standard
#' deviation or any other descriptive statistic by changing the line in the
#' code '\code{oumean = apply(outb, 2, mean)}' toward the end of this function.
#' A trimmed mean may be of interest when outliers are suspected.
#' @note require(np)
#' @note Make sure that functions wtdpapb, bigfp, stochdom2 are in the memory.
#' and options(np.messages=FALSE)
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso \code{\link{stochdom2}}
#' 
#' @references Vinod, H. D.", "Hands-On Intermediate Econometrics 
#' Using R"  (2008) World Scientific Publishers: Hackensack, NJ. (Chapter 4)
#' \url{https://www.worldscientific.com/worldscibooks/10.1142/12831}
#'
#' @concept stochastic dominance 
#' @concept financial portfolio choice
#' @examples
#'
#' set.seed(30)
#' xa=sample(20:30)#generally lower returns
#' xb=sample(32:40)# higher returns in xb
#' gp = comp_portfo2(xa, xb)#all Av(sdi) positive means xb dominates
#' ##positive SD1 to SD4 means xb dominates xa as it should
#'
#' @export

comp_portfo2 <- function(xa, xb) {
    # simplified:NAs already out
    gp = wtdpapb(xa, xb)  
    stdo2 = stochdom2(dj = gp$dj, wpa = gp$wpa, wpb = gp$wpb)
    outb = cbind(stdo2$sd1b, stdo2$sd2b, stdo2$sd3b, stdo2$sd4b)
    # print(outb) column sums for 4 orders of stochastic dominance
    oumean = apply(outb, 2, mean)
    return(oumean)
}
