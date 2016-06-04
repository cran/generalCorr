#' Approximate overall magnitudes of kernel regression partials dx/dy and dy/dx.
#'
#' Uses Vinod (2015) and runs kernel regression of x on y,  and also of y on x
#' by using the `np' package. The function goes on to compute a summary magnitude
#' of the overall approximate partial derivative dx/dy (and dy/dx), 
#' after adjusting for units by using
#' an appropriate  ratio of standard deviations.  Of course, 
#' the real partial derivatives of nonlinear functions
#' are generally  distinct for each observation.
#' 
#' @param x Vector of data on the dependent variable
#' @param y Vector of data on the regressor
#' @importFrom stats sd cor
#' @return vector of two magnitudes of kernel regression partials dx/dy and dy/dx.
#' @note This function is intended for use only after the direction of causal path
#' is already determined by various functions in this package (e.g. \code{somePairs}).
#' For example, if the researcher knows that x causes y, then only 
#' dy/dx denoted by dydx is relevant.
#' The other output of the function dxdy is to be ignored.
#' Similarly, only `dxdy' is relevant if y is known to be the cause of x.
#' 
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See  \code{\link{mag_ctrl}}.
#' @references Vinod, H. D.'Generalized Correlation and Kernel Causality with
#'  Applications in Development Economics' in Communications in
#'  Statistics -Simulation and Computation, 2015,
#'  \url{http://dx.doi.org/10.1080/03610918.2015.1122048}
#'
#' @references Vinod, H. D. 'Matrix Algebra Topics in Statistics and Economics
#' Using R', Chapter 4 in Handbook of Statistics: Computational Statistics
#' with R, Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#'
#' @keywords partial derivatives
#' @examples
#'
#' set.seed(123);x=sample(1:10);y=1+2*x+rnorm(10)
#' mag(x,y)#dxdy approx=.5 and dydx approx=2 will be nice.
#'
#' @export

mag = function(x, y) {
   sgn=sign(cor(x,y))
   sdy = sd(y)
    if (sdy <= 0) 
        stop("sd of y is zero")
    sdx = sd(x)
    if (sdx <= 0) 
        stop("sd of x is zero")
    mod.1 = kern(dep.y = x, reg.x = y)
    mod.2 = kern(dep.y = y, reg.x = x)
    R2xy = mod.1$R2
    dxdy = 0
    dydx = 0
    if (R2xy > 0) 
        dxdy = sgn*sqrt(R2xy) * sdx/sdy
    R2yx = mod.2$R2
    if (R2yx > 0) 
        dydx = sgn*sqrt(R2yx) * sdy/sdx
    effxyyx = c(dxdy, dydx)
    return(effxyyx)
}


