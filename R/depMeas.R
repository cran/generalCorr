#' depMeas Measure dependence between two vectors.
#'
#' An infant may depend on the mother for survival, but not vice versa.
#' Dependence relations need not be symmetric, yet correlation coefficients
#' are symmetric. One way to measure the extent of dependence is to find
#' the max of the absolute values of the two asymmetric correlations
#' using Vinod (2015) definition of generalized (asymmetric) correlation
#' coefficients.  It requires a kernel regression of x on y obtained by using 
#' the `np' package and its flipped version.  We use a block version of
#' `gmcmtx0'  called `gmcmtxBlk` to admit several bandwidths.
#' @param x {Vector of data on first variable}
#' @param y {Vector of data on second variable}
#' @param blksiz {block size, default=10, if chosen blksiz >n, where n=rows in matrix
#'      then blksiz=n. That is, no blocking is done}
#' 
#' @return A measure of dependence.
#' @note This function needs the gmcmtxBlk function which in turn needs the np package.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{gmcmtx0}} and \code{\link{gmcmtxBlk}}
#' @references Vinod, H. D. `Generalized Correlation and Kernel Causality with
#'  Applications in Development Economics' in Communications in
#'  Statistics -Simulation and Computation, 2015,
#'  \url{http://dx.doi.org/10.1080/03610918.2015.1122048}
#'
#' @references Vinod, H. D. 'Matrix Algebra Topics in Statistics and Economics
#' Using R', Chapter 4 in Handbook of Statistics: Computational Statistics
#' with R, Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#' @concept asymmetric  p-values
#' @examples
#' library(generalCorr)
#' options(np.messages = FALSE)
#' x=1:20;y=sin(x)
#' depMeas(x,y,blksiz=20)
#'
#' @export



depMeas = function(x, y, blksiz=length(x)) {
  g1 = gmcmtxBlk(cbind(x, y),blksiz=blksiz)
  sgn = sign(g1[1, 2])
  dep = sgn*max(abs(g1[1, 2]), abs(g1[2, 1]))
  return(dep)
}

