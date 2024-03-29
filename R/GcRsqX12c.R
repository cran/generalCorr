#' Generalized Granger-Causality. If dif>0, x2 Granger-causes x1.
#' 
#' The usual Granger-causality assumes linear regressions. This allows
#' nonlinear nonparametric kernel regressions using a local constat (lc) option.
#'  Calls GcRsqYXc for R square from kernel regression.
#'  R^2[x1=f(x1,x2)] choosing GcRsqYXc(y=x1, x=x2). The name `c' in the function
#'  refers to local constant option of kernel regressions.`
#'  It predicts x1 from 
#'  both x1 and x2 using all information till time (t-1).
#'  It also calls GcRsqYXc again after flipping x1 and x2.
#'  It returns RsqX1onX2, RsqX2onX1 and the difference dif=(RsqX1onX2-RsqX2onX1)
#'  If (dif>0) the regression x1=f(x1,x2) is better than the flipped
#'  version implying that x1 is more predictable or x2 Granger-causes x1
#'  x2 --> x1, rather than vice versa. The kernel regressions use
#'  regtype="lc" for local constant, bwmethod="cv.ls" for least squares-based
#'  bandwidth selection.
#' 
#' @param x1 {The data vector x1} 
#' @param px1 {number of lags of x1 in the data default px1=4} 
#' @param x2 {The data vector x2}
#' @param px2 {number of lags of x2 in the data, default px2=4} 
#' @param pwanted {number of lags of both x2 and x1 wanted for Granger causal analysis, default =4} 
#' @param ctrl {data matrix for designated control variable(s) outside causal paths
#' default=0 means no control variables are present}
#' @return This function returns 3 numbers: RsqX1onX2, RsqX2onX1 and dif
###' @importFrom stats complete.cases
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso \code{\link{causeSummary}} 
#' @seealso \code{\link{GcRsqYXc}} 
#' @references Vinod, H. D. `Generalized Correlation and Kernel Causality with
#'    Applications in Development Economics' in Communications in
#'    Statistics -Simulation and Computation, 2015,
#'    \doi{10.1080/03610918.2015.1122048}
#' @references Vinod, H. D. 'New exogeneity tests and causal paths,'
#'  Chapter 2 in 'Handbook of Statistics: Conceptual Econometrics 
#' Using R', Vol.32, co-editors: H. D. Vinod and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2019, pp. 33-64.
#'  
#' @references Vinod, H. D. Causal Paths and Exogeneity Tests 
#' in {Generalcorr} Package for Air Pollution and Monetary Policy 
#' (June 6, 2017). Available at SSRN: 
#' \url{https://www.ssrn.com/abstract=2982128}    
#' @references Zheng, S., Shi, N.-Z., Zhang, Z., 2012. 
#' Generalized measures of correlation for
#' asymmetry, nonlinearity, and beyond. Journal of the American Statistical
#' Association 107, 1239-1252.
###' -at-note internal routine
#' @return 
#' returns a list of 3 numbers. RsqX1onX2=(Rsquare of
#' kernel regression of X1 on X1 and X2), 
#' RsqX2onX1= (Rsquare of kernel regression of x2 on X2 and X1), and 
#' the difference between the two Rquares called dif
#' 
#' @examples
#'
#'
#' \dontrun{
#' library(Ecdat);options(np.messages=FALSE);attach(data.frame(MoneyUS))
#' GcRsqX12c(y,m)   
#' }
#' 
#' 
#' 
#' @export

GcRsqX12c = function(x1, x2, px1=4, px2=4, pwanted=4, ctrl = 0){
#print(head(x1,2))
#print(head(x2,2))
RsqX2onX1 =GcRsqYXc(x1, x2, px=px1, py=px2, pwanted=pwanted, ctrl = ctrl)
R21=RsqX2onX1[1]
RsqX1onX2=GcRsqYXc(x2, x1, px=px2, py=px1, pwanted=pwanted, ctrl = ctrl)     
R12=RsqX1onX2[1]
dif=R12-R21
list(RsqX1onX2=R12, RsqX2onX1=R21, dif=dif)
} #end of function
