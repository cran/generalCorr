#' Probability of unambiguously correct (+ or -) sign from bootPairs output
#' 
#' If there are p columns of data, \code{bootSign} produces a p-1 by 1 vector
#' of probabilities of correct signs assuming that the mean of n999 values
#' has the correct sign and assuming that m of the 'sum' index values inside the
#' range [-tau, tau] are neither positive nor negative but 
#' indeterminate or ambiguous (being too close to zero). That is,
#' the denominator of P(+1) or P(-1) is (n999-m) if m signs are too close to zero. 
#' Thus it measures the bootstrap success rate in identifying the correct sign, when the sign
#' of the average of n999 bootstraps is assumed to be correct.  
#' 
#' @param out {output from bootPairs with p-1 columns and n999 rows}
#' @param tau {threshold to determine what value is too close to
#'  zero, default tau=0.476 is equivalent to 15 percent threshold for 
#'  the unanimity index ui}
#' @return  sgn {When \code{mtx} has p columns, \code{sgn}
#' reports pairwise p-1 signs  representing 
#' (fixing the first column in each pair)
#' the average sign after averaging the
#' output of of \code{bootPairs(mtx)} (a n999 by p-1 matrix)
#' each containing resampled `sum' values summarizing the weighted sums 
#' associated with all three  criteria from the 
#' function \code{silentPairs(mtx)}
#' applied to each bootstrap sample separately.} #' 
#' 
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{silentPairs}}, \code{\link{bootQuantile}},
#' \code{\link{bootSignPcent}}.
#' @references Vinod, H. D. `Generalized Correlation and Kernel Causality with 
#'  Applications in Development Economics' in Communications in 
#'  Statistics -Simulation and Computation, 2015, 
#'  \doi{10.1080/03610918.2015.1122048} 
#' @references Vinod, H. D. and Lopez-de-Lacalle, J. (2009). 'Maximum entropy bootstrap
#'  for time series: The meboot R package.' Journal of Statistical Software,
#'  Vol. 29(5), pp. 1-19. 
#' @references Vinod, H. D. Causal Paths and Exogeneity Tests 
#' in {Generalcorr} Package for Air Pollution and Monetary Policy 
#' (June 6, 2017). Available at SSRN: \url{https://www.ssrn.com/abstract=2982128}
#' @concept bootstrap
#' @concept  meboot
#' @concept  kernel regression
#' @concept  pairwise comparisons
#' @examples
#' \dontrun{
#' options(np.messages = FALSE)
#' set.seed(34);x=sample(1:10);y=sample(2:11)
#' bb=bootPairs(cbind(x,y),n999=29)
#' bootSign(bb,tau=0.476) #gives success rate in n999 bootstrap sum computations
#' 
#' bb=bootPairs(airquality,n999=999);options(np.messages=FALSE)
#' bootSign(bb,tau=0.476)#signs for n999 bootstrap sum computations
#' 
#' data('EuroCrime');options(np.messages=FALSE)
#' attach(EuroCrime)
#' bb=bootPairs(cbind(crim,off),n999=29) #col.1= crim causes off 
#' #hence positive signs are more intuitively meaningful.
#' #note that n999=29 is too small for real problems, chosen for quickness here.
#' bootSign(bb,tau=0.476)#gives success rate in n999 bootstrap sum computations
#' }
#' @export

bootSign=function(out,tau=0.476) {
  n999=NROW(out)
  signNull=sign(apply(out,2,mean,na.rm=TRUE))#average out columns over n999 rows
  pm1=NCOL(out)  #p m=minus 1
  bootSign=rep(NA,pm1)
for (j in 1:(pm1)) {
  zj=out[,j] #one column at a time
  m = length(zj[abs(zj)<tau])
  if (signNull[j]==1) bootSign[j]= length(zj[zj > tau]) /(n999-m)
  if (signNull[j]==-1) bootSign[j]= length(zj[zj < -tau]) /(n999-m)
}  #end j loop
return(bootSign)
}