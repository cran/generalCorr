#' Nonlinear Granger causality between two time series workhorse function.(local 
#' constant version)
#' 
#' Function input is y=LHS=First time series and x=RHS=Second time series.
#' Kernel regression np package options regtype="lc" for local constant, 
#' and bwmethod="cv.ls" for least squares-based bandwidth selection are fixed.
#' Denote Rsq=Rsquare=R^2 in nonlinear kernel regression.
#' GcRsqYXc(.) computes the following two R^2 values.
#' out[1]=Rsqyyx = R^2 when we regress y on own lags of y and x.
#' out[2]=Rsqyy = R^2 when we regress y on own lags of y alone.
#'   
#' @param y {The data vector y for the Left side or dependent or first variable} 
#' @param py {number of lags of y in the data. px=4 for quarterly data} 
#' @param x {The data vector x for the right side or explanatory or second variable}
#' @param px {number of lags of x in the data} 
#' @param pwanted {number of lags of both x and y wanted for Granger causal analysis} 
#' @param ctrl {data matrix for designated control variable(s) outside causal paths
#' default=0 means no control variables are present}
#' @importFrom stats embed
#' @seealso \code{\link{GcRsqX12c}} 
#' @seealso \code{\link{kern_ctrl}} 
#' @return 
#' This function returns a set of 2 numbers measuring nonlinear Granger-causality 
#' for time series. out[1]=Rsqyyx, out[2]=Rsqyy.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
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
#' @note If data are annual or if no quarterly-type structure is present,
#' use this function with pwanted=px=py.  For example, the egg or chicken
#' data from lmtest package, Thurman W.N. and Fisher M.E. (1988)
#' 
#' @examples
#'
#'
#' \dontrun{
#' library(Ecdat);options(np.messages=FALSE);attach(data.frame(MoneyUS))
#' GcRsqYXc(y,m) 
#' }
#' 
#' 
#' 
#' @export

GcRsqYXc = function(y, x, px=4, py=4, pwanted=4, ctrl = 0)
   {
    # require(generalCorr); options(np.messages=FALSE)
    # if this is positive The c in YXc name is for local constant kernel
if (px != py) stop("px does not equal py, the number of lags")
if (pwanted > py) stop("pwanted lags exceeds lags in data")
    
    el=length(ctrl) # if el>1 control variables are present
#Task 1 compute corr and p-val with irwise deletion of NAs
    pxp1=px+1 #notation pxp1=px plus 1
    pxmpw=px-pwanted #notation m for minus, px minus pwanted
    newy=embed(y,pxp1) #last column has current values yt
    newx=embed(x,pxp1)#first col. has xt with px lags
# Task 2 compute kernel regression residual sum of squares of y on lagged y
 if(el>1) regyy=kern_ctrl(dep.y=newy[,pxp1], 
   reg.x=newy[,(pxmpw+1):px],ctrl=ctrl)
 if(el==1)   regyy=kern(dep.y=newy[,pxp1],#last column 
   reg.x=newy[,(pxmpw+1):px])
    Rsqyy=regyy$R2 #R Sq
# Task 3 compute kernel regression residualSumSq of y on lagged y & X
if(el>1)regyyx=kern_ctrl(dep.y=newy[,pxp1],#last column 
reg.x=cbind(newy[,(pxmpw+1):px],newx[,(pxmpw+1):pxp1]),ctrl=ctrl)
if(el==1)regyyx=kern(dep.y=newy[,pxp1],#last column 
 reg.x=cbind(newy[,(pxmpw+1):px],newx[,(pxmpw+1):pxp1]))
    Rsqyyx=regyyx$R2 # RSq
    out=rep(NA,2)
    out[1]=Rsqyyx
    out[2]=Rsqyy
    return(out)
}  #end function


