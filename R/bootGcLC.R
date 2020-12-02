#'  Compute vector of n999 nonlinear Granger causality paths
#' 
#' Maximum entropy bootstrap (meboot) package is used for statistical inference
#' The bootstrap output can be analyzed to estimate an approximate confidence
#' interval on sample-based direction of the causal path.
#' The LC in the function name stands for local constant.
#' Kernel regression np package options regtype="lc" for local constant, 
#' and bwmethod="cv.ls" for least squares-based bandwidth selection are fixed.
#'  
#' @param x1 {The data vector x1} 
#' @param px1 {number of lags of x1 in the data default px1=4} 
#' @param x2 {The data vector x2}
#' @param px2 {number of lags of x2 in the data, default px2=4} 
#' @param pwanted {number of lags of both x2 and x1 wanted for 
#'   Granger causal analysis, default =4} 
#' @param ctrl {data matrix having control variable(s) if any}
#' @param n999 {Number of bootstrap replications (default=9)}
#' @importFrom meboot meboot
#' @importFrom stats complete.cases
#' @return  out {is  n999 X 3 matrix for 3 outputs of GcauseX12 resampled} 
#' 
#' @note This computation is computer intensive and generally very slow. 
#'   It may be better to use this function
#'   it at a later stage in the investigation, after a preliminary 
#'   causal determination is already made. The 3 outputs of GauseX12 are
#'   two Rsquares and the difference between after subtracting the second
#'   from the first.  Col. 1 has (RsqX1onX2)
#'   Col.2 has (RsqX2onX1), and Col.3 has dif=(RsqX1onX2 -RsqX2onX1)
#'   Note that R-squares are always positive. 
#'   If dif>0, RsqX1onX2>RsqX2onX1, implying that x2 on RHS performs better
#'   that is, x2 --> x1 is the path, or x2 Granger-causes x1.
#'   If dif<0, x1 --> x2 holds. If dif is too close to zero,
#'   we may have bidirectional causality  x1 <--> x2. The proportion of
#'   resamples (out of n999) having dif<0 suggests level of confidence in
#'   the conclusion x1 --> x2.  The proportion of
#'   resamples (out of n999) having dif>0 suggests level of confidence in
#'   the conclusion x2 --> x1.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{GcRsqX12c}}.
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
#' library(Ecdat);options(np.messages=FALSE);attach(data.frame(MoneyUS))
#' bootGcLC(y,m,n999=9) 
#' }
#' \dontrun{
#' library(lmtest); data(ChickEgg);attach(data.frame(ChickEgg))
#' b2=bootGcLC(x1=chicken,x2=egg,pwanted=3,px1=3,px2=3,n999=99)
#' }
#' 
#' 
#' @export

bootGcLC = function(x1, x2, px2=4, px1=4, pwanted=4, 
              ctrl = 0, n999=9) {
  ok= complete.cases(x1,x2) 
  xx1=x1[ok]
  xx2=x2[ok]
  n=length(xx1)
  me1=meboot(xx1,reps=n999)$ensemble
  me2=meboot(xx2,reps=n999)$ensemble
  out=matrix(NA,nrow=n999,ncol=3)
  for (k in 1:n999) {
    outx = GcRsqX12c(x1=me1[,k], x2=me2[,k], 
            px2=px2, px1=px1, pwanted=pwanted, ctrl = ctrl)
    if(k==1) print(c("resample number k=",k))
    if(k==1) print(outx)
    out[k,1]=outx$RsqX1onX2
    out[k,2]=outx$RsqX2onX1
    out[k,3]=outx$dif
    if (k%%50 ==1) print(c("k=",k)) #track the progress 
  }
  print("dif>0 means Rsq-X1-on-X2 > Rsq-X2-on-X1 with causal path X2-->X1")
  colnames(out)=c("Rsq-X1-on-X2","Rsq-X2-on-X1","dif")
  Fn=function(x) quantile(x,probs=c(0.025,0.975),na.rm=TRUE)
  print("95 percent confidence intervals for each column")
  print(apply(out,2,Fn))
  print(apply(out,2,summary))
  bb=out[,3]
  px2to1=round(length(bb[bb>0])/n999,5)
  px1to2=round(length(bb[bb<0])/n999,5)
  print(c("boot prop. supporting x1-->x2",px1to2))
  print(c("boot prop. supporting x2-->x1",px2to1))
  return(out)
}
