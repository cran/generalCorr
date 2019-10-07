#' Kernel regression with control variables and optional residuals and gradients.
#' 
#' Allowing matrix input of control variables, this function runs kernel regression 
#' with options for residuals and gradients.
#' 
#' @param dep.y {Data on the dependent (response) variable}
#' @param reg.x {Data on the regressor (stimulus) variable}
#' @param ctrl {Data matrix on the control variable(s) kept outside the 
#'  causal paths.
#'  A constant vector is not allowed as a control variable.}
#' @param tol {Tolerance on the position of located minima of the cross-validation 
#'  function (default=0.1)}
#' @param ftol {Fractional tolerance on the value of cross validation function
#'  evaluated at local minima  (default=0.1)}
#' @param gradients {Set to TRUE if gradients computations are desired}
#' @param residuals {Set to TRUE if residuals are desired}
#' @importFrom np npreg npregbw
#' @return Creates a model object `mod' containing the entire kernel regression output.
#' If this function is called as \code{mod=kern_ctrl(x,y,ctrl=z)}, the researcher can
#' simply type \code{names(mod)} to reveal the large variety of outputs produced by `npreg' 
#' of the `np' package.
#' The user can access all of them at will using the dollar notation of R.
#' @note This is a work horse for causal identification.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See \code{\link{kern}}.
#' @references Vinod, H. D. `Generalized Correlation and Kernel Causality with 
#'  Applications in Development Economics' in Communications in 
#'  Statistics -Simulation and Computation, 2015, 
#'  \url{http://dx.doi.org/10.1080/03610918.2015.1122048} 
#' @concept apd  amorphous partial derivative
#' @concept  kernel regression residuals
#' @concept  kernel regression gradients
#' @examples
#' 
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:50],ncol=5)
#' require(np)
#' k1=kern_ctrl(x[,1],x[,2],ctrl=x[,4:5])
#' print(k1$R2) #prints the R square of the kernel regression
#' }
#' 
#' @export

kern_ctrl=
  function (dep.y, reg.x, ctrl, tol = 0.1, ftol = 0.1, gradients = FALSE, 
            residuals = FALSE)  #ctrl is a matrix of control variables
  {
    gr = FALSE
    resz = FALSE
    if (gradients) 
      gr = TRUE
    if (residuals) 
      resz = TRUE
    ox=naTriplet(x=dep.y,y=reg.x,ctrl=ctrl)
        bw = npregbw(ydat = as.vector(ox$newx), 
                xdat = cbind(ox$newy,ox$newctrl), 
                 tol = tol, ftol = ftol)
    mod = npreg(bws = bw, gradients = gr, residuals = resz)
    return(mod)
  }
