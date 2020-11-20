#' Block version abs_stdrhser Absolute residuals kernel regressions of standardized x on y and control 
#' variables, Cr1 has abs(Resid*RHS).
#'
#' 1) standardize the data to force mean zero and variance unity, 2) kernel
#' regress x on y and a matrix of control variables, 
#' with the option `residuals = TRUE' and finally 3) compute
#' the absolute values of residuals.
#'
#' The first argument is assumed to be the dependent variable.  If
#' \code{absBstdrhserC(x,y)} is used, you are regressing x on y (not the usual y
#' on x). The regressors can be a matrix with 2 or more columns. The missing values
#' are suitably ignored by the standardization.
#'
#' @param x {vector of data on the dependent variable}
#' @param y {data on the regressors which can be a matrix}
#' @param ycolumn {if y has more than one column, the 
#' column number used when multiplying residuals times
#' this column of y, default=1 or first column of y matrix is used}
#' @param ctrl {Data matrix on the control variable(s) beyond causal path issues}
#' @param blksiz {block size, default=10, if chosen blksiz >n, where n=rows in matrix
#'      then blksiz=n. That is, no blocking is done}
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
#'  \url{https://doi.org/gffn86}
#' @concept  kernel regression residuals
#' @examples
#'
#' \dontrun{
#' set.seed(330)
#' x=sample(20:50)
#' y=sample(20:50)
#' z=sample(21:51)
#' absBstdrhserC(x,y,ctrl=z)
#' }
#'
#' @export


absBstdrhserC=
  function (x, y, ctrl, ycolumn=1, blksiz=10) 
  {
    stdx = function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
    stx = (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
    if (NCOL(x)>1) stop("too many columns of x in absBstdrhserC")
    p = NCOL(y)
    q=0
    if(length(ctrl)>1)  q = NCOL(ctrl)
    n = NROW(y)
    if (blksiz>n) blksiz=n
    ge=getSeq(n,blksiz=blksiz)
    ares=rep(NA,n) #absolute residuals vector
    LO=ge$sqLO
    UP=ge$sqUP
    k=length(LO)
    for (ik in 1:k){
      L1=LO[ik]  
      U1=UP[ik]
      stxx=stdx(x[L1:U1])     
      if (p == 1) {
        yy=y[L1:U1]
      sty = (yy - mean(yy, na.rm = TRUE))/sd(yy, na.rm = TRUE)}
    if (p > 1) {
      yy=y[L1:U1,] #pick all columns of y matrix
      sty = apply(yy, 2, stdx)} #name of the function is stdx
    if (q == 1) {
      ctrl2=ctrl[L1:U1]
      stz = (ctrl2 - mean(ctrl2, na.rm = TRUE))/sd(ctrl2, na.rm = TRUE)}
    if (q > 1) {
      ctrl2=ctrl[L1:U1,]#all columns of control variables
      stz = apply(ctrl2, 2, stdx)}
    if (q>=1){
      kk1 = kern_ctrl(dep.y = stxx, reg.x = sty, ctrl = stz, residuals = TRUE)}
    if (q==0)kk1=kern(dep.y = stxx, reg.x = sty, residuals = TRUE)
    if (p==1)   ares[L1:U1] = abs(sty*kk1$resid) #RHS var * residual here
    if (p>1)   ares[L1:U1] = abs(sty[,ycolumn]*kk1$resid)
    } #end ik loop over blocks
    return(ares)
  }
