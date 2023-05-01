#' Generalized Granger-Causality. If dif>0, x2 Granger-causes x1.
#' 
#' The usual Granger-causality assumes linear regressions. This function allows
#' nonlinear nonparametric kernel regressions using a local linear (ll) option.
#' Granger-causality (Gc) is generalized using nonlinear kernel 
#' regressions using local linear (ll) option. This
#' functionn computes two R^2 values. (i) R12 or kernel regression  R^2 of x1t on its
#' own lags and x2t and its lags. (ii) R21 or kernel regression R^2 of x2t on its
#' own lags and x1t and its lags. (iii) dif=R12-R21, the difference between the
#' two R^2 values. If dif>0 then x2 Granger-causes x1.
#' 
#'  Calls GcRsqYX for R-square from kernel regression (local linear version)
#'  R^2[x1=f(x1,x2)] choosing GcRsqYX(y=x1, x=x2). It predicts x1 from 
#'  both x1 and x2 using all information till time (t-1).
#'  It also calls GcRsqYX again after flipping x1 and x2.
#'  It returns RsqX1onX2, RsqX2onX1 and the difference dif=(RsqX1onX2-RsqX2onX1)
#'  If (dif>0) the regression x1=f(x1,x2) is better than the flipped
#'  version implying that x1 is more predictable or x2 Granger-causes x1,
#'  x2 --> x1, rather than vice versa. The kernel regressions use
#'  regtype="ll" for local linear, bwmethod="cv.aic" for AIC-based
#'  bandwidth selection.
#' 
#' @param x1 {The data vector x1} 
#' @param px1 {The number of lags of x1 in the data default px1=4} 
#' @param x2 {The data vector x2}
#' @param px2 {The number of lags of x2 in the data, default px2=4} 
#' @param pwanted {number of lags of both x2 and x1 wanted for Granger causal analysis, default =4} 
#' @param ctrl {data matrix for designated control variable(s) outside causal paths
#' default=0 means no control variables are present}
#' @return This function returns 3 numbers: RsqX1onX2, RsqX2onX1 and dif
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso \code{\link{bootGcRsq}}, 
#'  \code{\link{causeSummary}}, 
#'  \code{\link{GcRsqYX}}. 
#' @references Vinod, H. D. `Generalized Correlation and Kernel Causality with
#'    Applications in Development Economics' in Communications in
#'    Statistics -Simulation and Computation, 2015,
#'    \doi{10.1080/03610918.2015.1122048}
#' @references Vinod, H. D. 'New exogeneity tests and causal paths,'
#'  Chapter 2 in 'Handbook of Statistics: Conceptual Econometrics 
#' Using R', Vol.32, co-editors: H. D. Vinod and C.R. Rao. New York:
#' North-Holland, Elsevier Science Publishers, 2019, pp. 33-64.
#'  
#' @references Vinod, H. D. Causal Paths and Exogeneity Tests 
#' in {Generalcorr} Package for Air Pollution and Monetary Policy 
#' (June 6, 2017). Available at SSRN: 
#' \url{https://www.ssrn.com/abstract=2982128} 
#'    
#' @references Zheng, S., Shi, N.-Z., Zhang, Z., 2012. 
#' Generalized measures of correlation for
#' asymmetry, nonlinearity, and beyond. Journal of the American Statistical
#' Association 107, 1239-1252.
###' -at-note internal routine
#' @return 
#' returns a list of 3 numbers. RsqX1onX2=(Rsquare of
#' kernel regression of X1 on lags of X1 and X2 and its lags), 
#' RsqX2onX1= (Rsquare of kernel regression of x2 on own lags of X2 and X1), and 
#' the difference between the two Rquares (first minus second) called `dif.'
#' If dif>0 then x2 Granger-causes x1
#' @examples
#'
#'
#' \dontrun{
#' library(Ecdat);options(np.messages=FALSE);attach(data.frame(MoneyUS))
#' GcRsqX12(y,m)   
#' }
#' 
#' 
#' 
#' @export

GcRsqX12 = function(x1, x2, px1=4, px2=4, pwanted=4, ctrl = 0){
#print(head(x1,2))
#print(head(x2,2))
RsqX2onX1 =GcRsqYX(x1, x2, px=px2, py=px1, pwanted=pwanted, ctrl = ctrl)
R21=RsqX2onX1[1]
RsqX1onX2=GcRsqYX(x2, x1, px=px1, py=px2, pwanted=pwanted, ctrl = ctrl)     
R12=RsqX1onX2[1]
dif=R12-R21
list(RsqX1onX2=R12, RsqX2onX1=R21, dif=dif)
} #end of function
