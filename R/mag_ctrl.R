#' After removing control variables, magnitude of effect of x on y, and of y on x.
#'
#' Uses Vinod (2015) and runs kernel regressions: \code{x~ y + ctrl}
#' and  \code{x~ ctrl} to evaluate the `incremental change' in R-squares.
#' Let (rxy;ctrl) denote the square root of that `incremental change' after its sign is made the
#' same as that of the Pearson correlation coefficient from
#' \code{cor(x,y)}). One can interpret (rxy;ctrl) as
#' a generalized partial correlation coefficient when x is regressed on y after removing
#' the effect of control variable(s) in \code{ctrl}.  It is more general than the usual partial 
#' correlation coefficient, since this one
#' allows for nonlinear relations among variables. 
#' Next, the function computes `dxdy' obtained by multiplying (rxy;ctrl) by the ratio of
#' standard deviations, \code{sd(x)/sd(y)}. Now our `dxdy' approximates the magnitude of the
#' partial derivative (dx/dy) in a causal model where y is the cause and x is the effect.
#' The function also reports entirely analogous `dydx' obtained by interchanging x and y.
#'  
#' @param x Vector of data on the dependent variable.
#' @param y Vector of data on the regressor.
#' @param ctrl {data matrix for designated control variable(s) outside causal paths.
#'  A constant vector is not allowed as a control variable.}
#' @importFrom stats sd cor
#' @return vector of two magnitudes `dxdy' (effect when x is regressed on y) and 
#' `dydx' for reverse regression.  Both regressions remove the effect of control variable(s).
#' @note This function is intended for use only after the causal path direction 
#' is already determined by various functions in this package (e.g. \code{someCPairs}).
#' That is, after the researcher knows whether x causes y or vice versa.
#' The output of this function is a vector of two numbers: (dxdy, dydx), in that order,
#' representing the magnitude of effect of one variable on the other.
#' We expect the researcher to use only `dxdy' if y is the 
#' known cause, or `dydx' if x is the cause. These approximate overall measures
#' may not be well-defined in some applications, because 
#' the real partial derivatives of nonlinear functions
#' are generally  distinct for each evaluation point.
#' 
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See  \code{\link{mag}}
#' @references Vinod, H. D. `Generalized Correlation and Kernel Causality with
#'  Applications in Development Economics' in Communications in
#'  Statistics -Simulation and Computation, 2015,
#'  \doi{10.1080/03610918.2015.1122048} 
#'
#' @references Vinod, H. D. 'Matrix Algebra Topics in Statistics and Economics
#' Using R', Chapter 4 in Handbook of Statistics: Computational Statistics
#' with R, Vol.32, co-editors: M. B. Rao and C. R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#' 
#' @concept apd amorphous partial derivatives
#' @examples
#'
#' set.seed(123);x=sample(1:10); z=runif(10); y=1+2*x+3*z+rnorm(10)
#' options(np.messages=FALSE)
#' mag_ctrl(x,y,z)#dx/dy=0.47 is approximately 0.5, but dy/dx=1.41 is not approx=2,
#'
#' @export

mag_ctrl = function(x, y, ctrl) {
    sgn=sign(cor(x,y))
    sdy = sd(y)
    if (sdy <= 0) 
        stop("sd of y is zero")
    sdx = sd(x)
    if (sdx <= 0) 
        stop("sd of x is zero")
    
    sdall = apply
    mod.1 = kern(dep.y = x, reg.x = cbind(y, ctrl))
    mod.2 = kern(dep.y = y, reg.x = cbind(x, ctrl))
    mod.3 = kern(dep.y = y, reg.x = ctrl)
    mod.4 = kern(dep.y = x, reg.x = ctrl)
    Gxy = mod.1$R2 - mod.4$R2
    dxdy = 0
    dydx = 0
    if (Gxy > 0) 
        dxdy = sgn*sqrt(Gxy) * sdx/sdy
    Gyx = mod.2$R2 - mod.3$R2
    if (Gyx > 0) 
        dydx = sgn*sqrt(Gyx) * sdy/sdx
    effxyyx = c(dxdy, dydx)
    return(effxyyx)
}



