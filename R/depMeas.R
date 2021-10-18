#' depMeas Signed measure of nonlinear nonparametric dependence between two vectors.
#'
#' An infant may depend on the mother for survival, but not vice versa.
#' Dependence relations need not be symmetric, yet correlation coefficients
#' are symmetric. One way to measure the extent of dependence is to find
#' the max of the absolute values of the two asymmetric correlations
#' using Vinod (2015) definition of generalized (asymmetric) correlation
#' coefficients.  It requires a kernel regression of x on y obtained by using 
#' the `np' package and its flipped version.  We use a block version of
#' `gmcmtx0'  called `gmcmtxBlk` to admit several bandwidths.
#' @param x {Vector of data on the first variable}
#' @param y {Vector of data on the second variable}
#' @param blksiz {block size, default blksiz =n, where n=rows in the matrix
#'      or no blocking is done}
#' 
#' @return A measure of dependence having the same sign as Pearson correlation. Its
#' magnitude equals the larger of the two generalized correlation coefficients
#' @note This function needs the gmcmtxBlk function, which in turn needs the np package.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{gmcmtx0}} and \code{\link{gmcmtxBlk}}
#' @references Vinod, H. D. `Generalized Correlation and Kernel Causality with
#'  Applications in Development Economics' in Communications in
#'  Statistics -Simulation and Computation, 2015,
#'  \url{https://doi.org/gffn86}
#'
#' @references Vinod, H. D. 'Matrix Algebra Topics in Statistics and Economics
#' Using R', Chapter 4 in Handbook of Statistics: Computational Statistics
#' with R, Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#' 
#' @references Vinod, H. D. (2021) 'Generalized, Partial and Canonical Correlation
#' Coefficients' Computational Economics, 59(1), 1--28.
#' 
#' @concept asymmetric  p-values
#' @examples
#' library(generalCorr)
#' options(np.messages = FALSE)
#' x=1:20;y=sin(x)
#' depMeas(x,y,blksiz=20)
#'
#' @export



depMeas = function(x, y, blksiz=length(x)) {
  n=length(x)
  if (blksiz<=1){ blksiz=n
  print("bad blksiz reset to n in depMeas")}
  if (blksiz>n) {blksiz=n
  print("bad blksiz reset to n in depMeas")}
  g1 = gmcmtxBlk(cbind(x, y),blksiz=blksiz)
  sgn = sign(g1[1, 2])
  dep = sgn*max(abs(g1[1, 2]), abs(g1[2, 1]))
  return(dep)
}

