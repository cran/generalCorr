#' Block version of Absolute values of residuals of kernel regressions of standardized x on 
#' standardized y and control variables.
#'
#' 1) standardize the data to force mean zero and variance unity, 2) kernel
#' regress x on y and a matrix of control variables, 
#' with the option `residuals = TRUE' and finally 3) compute
#' the absolute values of residuals.
#'
#' The first argument is assumed to be the dependent variable.  If
#' \code{abs_stdres(x,y)} is used, you are regressing x on y (not the usual y
#' on x). The regressors can be a matrix with two or more columns. The missing values
#' are suitably ignored by the standardization.
#'
#' @param x {vector of data on the dependent variable}
#' @param y {data on the regressors which can be a matrix}
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
#' @references Vinod, H. D.'Generalized Correlation and Kernel Causality with
#'  Applications in Development Economics' in Communications in
#'  Statistics -Simulation and Computation, 2015,
#'  \doi{10.1080/03610918.2015.1122048}
#' @concept  kernel regression residuals
#' @examples
#'
#' \dontrun{
#' set.seed(330)
#' x=sample(20:50)
#' y=sample(20:50)
#' z=sample(21:51)
#' absBstdresC(x,y,ctrl=z)
#' }
#'
#' @export


absBstdresC=
  function (x, y, ctrl, blksiz=10) 
  {
    stdx = function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
    stx = (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
    if (NCOL(x)>1) stop("too many columns of x in absBstdresC")
    p = NCOL(y)
    q=0 #ctrl=0, or when no control variables are present
    if(length(ctrl)>1)  q = NCOL(ctrl)#ctrl is a mtrix of control variables
    n = NROW(y)
    if (blksiz>n) blksiz=n
    ge=getSeq(n,blksiz=blksiz)
    ares=rep(NA,n) #storage for absolute residuals vector
    LO=ge$sqLO
    UP=ge$sqUP
    k=length(LO) #number of blocks is k
    for (ik in 1:k){
      L1=LO[ik]  
      U1=UP[ik]
      stxx=stdx(x[L1:U1])
    if (p == 1){
      yy=y[L1:U1]
      sty = (yy - mean(yy, na.rm = TRUE))/sd(yy, na.rm = TRUE)}
    if (p > 1) {
      yy=y[L1:U1,]
      sty = apply(yy, 2, stdx)}
       if (q == 1){ 
      ctrlz=ctrl[L1:U1]
      stz = (ctrlz - mean(ctrlz, na.rm = TRUE))/sd(ctrlz, na.rm = TRUE)}
    if (q > 1){ 
      ctrlz=ctrl[L1:U1,]
      stz = apply(ctrlz, 2, stdx)}
    if ( q>=1){
       kk1 = kern_ctrl(dep.y = stxx, reg.x = sty, ctrl = stz, residuals = TRUE)}
    if (q==0) kk1 = kern(dep.y = stxx, reg.x = sty, residuals = TRUE)
    ares[L1:U1] = abs(kk1$resid)
    } #end ik loop over blocks
    return(ares)
  }
