#' Absolute values of Hausman-Wu null in kernel regressions of x on y when
#' both x and y are standardized.
#'
#' 1) standardize the data to force mean zero and variance unity, 2) kernel
#' regress x on y, with the option `gradients = TRUE' and finally 3) compute
#' the absolute values of Hausman-Wu null hypothesis for testing exogeneity,
#' or E(RHS.regressor*error)=0 where error is approximated by kernel 
#' regression residuals
#'
#' The first argument is assumed to be the dependent variable.  If
#' \code{abs_stdrhserr(x,y)} is used, you are regressing x on y (not the usual y
#' on x). The regressors can be a matrix with 2 or more columns. The missing values
#' are suitably ignored by the standardization.
#'
#' @param x {vector of data on the dependent variable}
#' @param y {data on the regressors which can be a matrix}
#' @importFrom stats sd
#' @return Absolute values of kernel regression RHS*residuals are returned after
#' standardizing the data on both sides so that the magnitudes of 
#' Hausman-Wu null values are comparable between regression of x on y on
#' the one hand and flipped regression of y on x on the other.
## @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#'
#' @concept  kernel regression
#' @concept  Hausman-Wu statistic
#' @examples
#' \dontrun{
#' set.seed(330)
#' x=sample(20:50)
#' y=sample(20:50)
#' abs_stdrhserr(x,y)
#' }
#' @export


abs_stdrhserr = function (x, y) 
  #standardized data, absolute value of residuals* rhs variable y here
{
   kk1=stdres(x,y)
    astdrhserr = abs(kk1*y)
    return(astdrhserr)
}
