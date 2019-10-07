#' Absolute values of gradients (apd's) of kernel regressions of x on y when
#' both x and y are standardized and control variables are present.
#'
#' 1) standardize the data to force mean zero and variance unity, 2) kernel
#' regress x on y and a matrix of control variables, 
#' with the option `gradients = TRUE' and finally 3) compute
#' the absolute values of gradients
#'
#' The first argument is assumed to be the dependent variable.  If
#' \code{abs_stdapdC(x,y)} is used, you are regressing x on y (not the usual y
#' on x). The regressors can be a matrix with 2 or more columns. The missing values
#' are suitably ignored by the standardization.
#'
#' @param x {vector of data on the dependent variable}
#' @param y {data on the regressors which can be a matrix}
#' @param ctrl {Data matrix on the control variable(s) beyond causal path issues}
#' @importFrom stats sd
#' @return Absolute values of kernel regression gradients are returned after
#' standardizing the data on both sides so that the magnitudes of amorphous
#' partial derivatives (apd's) are comparable between regression of x on y on
#' the one hand and regression of y on x on the other.
## @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See  \code{\link{abs_stdapd}}.
#'
#' @concept  kernel regression gradients
#' @concept apd
#' @examples
#' \dontrun{
#' set.seed(330)
#' x=sample(20:50)
#' y=sample(20:50)
#' z=sample(20:50)
#' abs_stdapdC(x,y,ctrl=z)
#' }
#' @export


abs_stdapdC=
  function (x, y, ctrl) 
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
    kk1 = kern_ctrl(dep.y = as.vector(stx), reg.x = sty, ctrl=stz, gradients = TRUE)
    agrad = abs(kk1$grad[,1])
    return(agrad)
  }
