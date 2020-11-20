#' Compute matrix of n999 rows and p-1 columns of bootstrap `sum' 
#' (strength from Cr1 to Cr3).
#' 
#' Maximum entropy bootstrap (meboot) package is used for statistical inference
#' using the sum of three signs sg1 to sg3 from the three criteria Cr1 to Cr3 to
#' assess preponderance of evidence in favor of a sign. (+1, 0, -1).
#' The bootstrap output can be analyzed to assess approximate
#' preponderance of a particular sign which determines
#' the causal direction.
#' 
#' @param mtx {data matrix with two or more columns}
#' @param ctrl {data matrix having control variable(s) if any}
#' @param n999 {Number of bootstrap replications (default=9)}
#' @importFrom meboot meboot
#' @importFrom stats complete.cases
#' @return  out {When \code{mtx} has p columns, \code{out}
#' of \code{bootPairs(mtx)} is a matrix of n999 rows and p-1 columns
#' each containing resampled `sum' values summarizing the weighted sums 
#' associated with all three  criteria from the function \code{silentPairs(mtx)}
#' applied to each bootstrap sample separately.} 
#' 
#' @note This computation is computer intensive and generally very slow. 
#'   It may be better to use
#'   it at a later stage in the investigation when a preliminary 
#'   causal determination 
#'   is already made.
#' A positive sign for j-th weighted sum reported in the column `sum' means
#' that the first variable listed in the argument matrix \code{mtx} is the 
#' `kernel cause' of the variable in the (j+1)-th column of \code{mtx}.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{silentPairs}}.
#' @references Vinod, H. D. `Generalized Correlation and Kernel Causality with 
#'  Applications in Development Economics' in Communications in 
#'  Statistics -Simulation and Computation, 2015, 
#'  \url{https://doi.org/gffn86} 
#' @references Zheng, S., Shi, N.-Z., and Zhang, Z. (2012). Generalized measures 
#'  of correlation for asymmetry, nonlinearity, and beyond. 
#'  Journal of the American Statistical Association, vol. 107, pp. 1239-1252.
#' @references Vinod, H. D. and Lopez-de-Lacalle, J. (2009). 'Maximum entropy bootstrap
#'  for time series: The meboot R package.' Journal of Statistical Software,
#'  Vol. 29(5), pp. 1-19. 
#' @references Vinod, H. D. Causal Paths and Exogeneity Tests 
#' in {Generalcorr} Package for Air Pollution and Monetary Policy 
#' (June 6, 2017). Available at SSRN: \url{https://www.ssrn.com/abstract=2982128}
#' @concept maximum entropy bootstrap
#' @examples
#' \dontrun{
#' options(np.messages = FALSE)
#' set.seed(34);x=sample(1:10);y=sample(2:11)
#' bb=bootPairs(cbind(x,y),n999=29)
#' apply(bb,2,summary) #gives summary stats for n999 bootstrap sum computations
#' 
#' bb=bootPairs(airquality,n999=999);options(np.messages=FALSE)
#' apply(bb,2,summary) #gives summary stats for n999 bootstrap sum computations
#' 
#' data('EuroCrime')
#' attach(EuroCrime)
#' bootPairs(cbind(crim,off),n999=29)#First col. crim causes officer deployment,
#' #hence positives signs are most sensible for such call to bootPairs
#' #note that n999=29 is too small for real problems, chosen for quickness here.
#' }
#' @export

bootPairs = function(mtx, ctrl = 0, n999 = 9) {
  ok= complete.cases(mtx) 
  p = NCOL(mtx[ok,])
  n = NROW(mtx[ok,])
  out = matrix(NA, nrow = n999, ncol = p - 1)
  Memtx <- array(NA, dim = c(n, n999, p))  #3 dimensional matrix
  for (i in 1:p) {
    Memtx[, , i] = meboot(x=mtx[ok, i], reps = n999)$ensem
  }
  for (k in 1:n999) {
    out[k, ] = silentPairs(mtx = Memtx[, k, 1:p], ctrl = ctrl)
    if (k%%50 ==1) print(c("k=",k)) #track the progress 
  }
  colnames(out) =colnames(mtx)[2:p]
  return(out)
}
