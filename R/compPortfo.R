#' Compares two vectors (portfolios) using
#' momentVote, DecileVote and exactSdMtx functions.
#'
#' Given two vectors of portfolio returns this function summarizes their ranks
#' based on moments, deciles and exact measures of stochastic dominance.
#' as explained in Vinod (2021). This algorithm has model selection applications.
#'
#' @param xa {Data on returns for portfolio A in the form of a T by 1 vector}
#' @param xb {Data on returns for portfolio B in the form of a T by 1 vector}
#' @return Returns three numbers which represent signs based differences in
#' ranks (rank=1 for most desirable) measured by [rank(xa)-rank(xb)] using
#' momentVote, decileVote, and  exactSdMtx which are weighted
#' averages of four moments, nine deciles and exact measures of stochastic
#' dominance (from ECDFs of four orders, SD1 to SD4) respectively.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso \code{\link{exactSdMtx}}
#' @seealso \code{\link{momentVote}}
#' @seealso \code{\link{decileVote}}
#' @note There are model-selection applications where two models A and B are
#' compared and one wants to choose the model smaller absolute value of
#' residuals. This function when applied for model-selection will have
#' he inputs xa and xb as absolute residuals. We can compare the entire
#' probability distributions of absolute residuals by moments, deciles
#' or SD1 to SD4. Of course, care must be taken to choose xa or
#' xb depending on which model has smaller absolute residuals. This choice
#' is the exact opposite of portfolio choice application where
#' larger return is more desirable.  \code{silentPair2()}
#' and \code{siPair2Blk} call this
#' function for model selection application.
#' 
#' @references Vinod, H. D.", "Hands-On Intermediate Econometrics 
#' Using R"  (2008) World Scientific Publishers: Hackensack, NJ. (Chapter 4)
#' \url{https://www.worldscientific.com/worldscibooks/10.1142/6895}
#' @references Vinod, Hrishikesh D., R Package GeneralCorr 
#' Functions for Portfolio Choice 
#' (November 11, 2021). Available at SSRN: 
#' https://ssrn.com/abstract=3961683 
#'
#' @concept stochastic dominance 
#' @concept financial portfolio choice
#' @examples
#'
#' set.seed(30)
#' xa=sample(20:30)#generally lower returns
#' xb=sample(32:40)# higher returns in xb
#' gp = compPortfo(xa, xb)#all Av(sdi) positive means xb dominates
#' ##output (1,1,1) means xb dominates xa. xb are larger by consruction
#'
#' @export

compPortfo <- function(xa, xb) {
  # simplified:NAs already out
  #  we choose largest choice rank or smallest abs.resid
  xab=cbind(xa,xb)
  ouall=rep(NA,3)
  m1=momentVote(xab)
  mm1=NROW(m1)
  m2=m1[mm1,1:2] #evaluate last row
  ouall[1]=as.numeric(sign(m2[1]-m2[2]))
  d1=decileVote(xab)
  dd1=d1$out
  d2=dd1[NROW(dd1),1:2]
  ouall[2]=as.numeric(sign(d2[1]-d2[2]))
  s1=exactSdMtx(xab)
 # source("c:/r-hdvall/summaryRank.R")
  s11=summaryRank(s1$out)
  s2=s11[10,]  #tenth row has choice, larger is smaller 
  ouall[3]=as.numeric(sign(s2[1]-s2[2]))
  return(ouall)
}
