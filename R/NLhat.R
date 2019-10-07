#' Compute fitted values from kernel regression of x on y and y on x
#' 
#' This is an auxiliary function for `gmcmtxBlk.' It uses 
#' two numerical vectors (x, y) of same length to create two vectors
#' (xhat, yhat) of fitted values using nonlinear kernel regressions.
#' It
#' uses package `np' called by 
#' \code{kern} function to kernel regress x on y, and conversely y on x. 
#' It uses the option `residuals=TRUE' of `kern'
#' 
#' @param x {A column vector of x data}
#' @param y {A column vector of y data}
#' @return {two vectors named xhat and yhat for fitted values}
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also as \code{\link{gmcmtxBlk}}.

#' @examples
#' 
#' 
#' \dontrun{
#' set.seed(34);x=sample(1:15);y=sample(1:15)
#' NLhat(x,y)}
#' 
#' @export


NLhat=function(x,y){
  mod1=kern(dep.y=x,reg.x=y,residuals=TRUE)
  mod2=kern(dep.y=y,reg.x=x,residuals=TRUE)
  xhat=x-mod1$resid
  yhat=y-mod2$resid
  list(xhat=xhat,yhat=yhat)}


