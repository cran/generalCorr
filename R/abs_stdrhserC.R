#' Absolute residuals kernel regressions of standardized x on y and control 
#' variables, Cr1 has abs(RHS*y) not gradients.
#'
#' 1) standardize the data to force mean zero and variance unity, 2) kernel
#' regress x on y and a matrix of control variables, 
#' with the option `residuals = TRUE' and finally 3) compute
#' the absolute values of residuals.
#'
#' The first argument is assumed to be the dependent variable.  If
#' \code{abs_stdrhserC(x,y)} is used, you are regressing x on y (not the usual y
#' on x). The regressors can be a matrix with 2 or more columns. The missing values
#' are suitably ignored by the standardization.
#'
#' @param x {vector of data on the dependent variable}
#' @param y {data on the regressors which can be a matrix}
#' @param ycolumn {if y has more than one column, the 
#' column number used when multiplying residuals times
#' this column of y, default=1 or first column of y matrix is used}
#' @param ctrl {Data matrix on the control variable(s) beyond causal path issues}
#' @importFrom stats sd
#' @return Absolute values of kernel regression residuals are returned after
#' standardizing the data on both sides so that the magnitudes of residuals are
#' comparable between regression of x on y on the one hand and regression of y
#' on x on the other.
### @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See  \code{\link{abs_stdres}}.
#' @references Vinod, H. D. 'Generalized Correlation and Kernel Causality with
#'  Applications in Development Economics' in Communications in
#'  Statistics -Simulation and Computation, 2015,
#'  \url{http://dx.doi.org/10.1080/03610918.2015.1122048}
#' @concept  kernel regression residuals
#' @examples
#'
#' \dontrun{
#' set.seed(330)
#' x=sample(20:50)
#' y=sample(20:50)
#' z=sample(21:51)
#' abs_stdrhserC(x,y,ctrl=z)
#' }
#'
#' @export


abs_stdrhserC=
  function (x, y, ctrl, ycolumn=1) 
  {
    stdx = function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
    stx = (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
    p = NCOL(y)
    q = NCOL(ctrl)#ctrl is a mtrix of control variables
    if (p == 1) 
      sty = (y - mean(y, na.rm = TRUE))/sd(y, na.rm = TRUE)
    if (p > 1) 
      sty = apply(y, 2, stdx)
    if (q == 1) 
      stz = (ctrl - mean(ctrl, na.rm = TRUE))/sd(ctrl, na.rm = TRUE)
    if (q > 1) 
      stz = apply(ctrl, 2, stdx)
    kk1 = kern_ctrl(dep.y = stx, reg.x = sty, ctrl = stz, residuals = TRUE)
    if (p==1)   ares = abs(sty*kk1$resid)
    if (p>1)   ares = abs(sty[,ycolumn]*kk1$resid)
    return(ares)
  }
