#' Summary magnitudes after removing control variables in several pairs where dependent 
#' variable is fixed.
#'
#' This builds on the function \code{mag_ctrl}, where the input matrix \code{mtx}
#' has p columns. The first column is present in each of the (p-1) pairs. Its
#' output is a matrix with four columns containing the names of variables 
#' and approximate overall estimates of the magnitudes of
#' partial derivatives (dy/dx) and (dx/dy) for a distinct (x,y) pair in a row.  
#' The estimated overall derivatives are not always well-defined, because 
#' the real partial derivatives of nonlinear functions
#' are generally  distinct for each observation point. 
#' 
#' The function \code{mag_ctrl} has kernel regressions: \code{x~ y + ctrl}
#' and  \code{x~ ctrl} to evaluate the`incremental change' in R-squares.
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
#' \code{someMegPairs} function runs  the function \code{mag_ctrl} on several column
#' pairs in a matrix input \code{mtx} where the first column is held fixed and all others
#' are changed one by one, reporting two partial derivatives for each row. 
#' 
#' @param mtx {The data matrix with many columns where the first column is fixed and then 
#'  paired with all other columns, one by one.}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths.
#'  A constant vector is not allowed as a control variable.}
#' @param dig {Number of digits for reporting (default \code{dig}=6).}
#' @param verbo {Make \code{verbo= TRUE} for printing detailed steps.}
#' @return Table containing names of Xi and Xj and two magnitudes: (dXidXj, dXjdXi).
#' dXidXj is the magnitude of the effect on Xi when Xi is regressed on Xj 
#' (i.e., when Xj is the cause).  The analogous dXjdXi is the magnitude
#' when Xj is regressed on Xi.
#' @note This function is intended for use only after the causal path direction 
#' is already determined by various functions in this package (e.g. \code{someCPairs}).
#' That is, after the researcher knows whether Xi causes Xj or vice versa.
#' The output of this function is a matrix of 4 columns, where first columns list
#' the names of Xi and Xj and the next two numbers in each row are
#' dXidXj, dXjdXi, respectively,
#' representing the magnitude of effect of one variable on the other.
#' 
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See  \code{\link{mag_ctrl}}, \code{\link{someCPairs}}
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
#' @concept  partial derivatives 
#' @examples
#'
#'   set.seed(34);x=sample(1:10);y=1+2*x+rnorm(10);z=sample(2:11)
#'   w=runif(10)
#'   ss=someMagPairs(cbind(y,x,z),ctrl=w)
#'
#' @export

someMagPairs=
  function (mtx, ctrl, dig = 6, verbo = TRUE) 
  {
    n = NROW(mtx)
    p = NCOL(mtx)
    k = NCOL(ctrl) #set of column numbers representing control variables
    npair = p - 1 
    out = matrix(NA, nrow=npair, ncol=4)
    colnames(out) = c("Xi", "Xj", "dXi/dXj", "dXj/dXi")
    nam = colnames(mtx)
    if (length(nam) == 0) 
      nam = paste("V", 1:p, sep = "")
    ii = 0
    for (i in 2:p) {
      y0 = mtx[, i]
      x0 = mtx[, 1]  #x0 is first column
      z0 = ctrl
      na2 = naTriplet(x0, y0, z0)#triplet-wise delete missing data
      x = na2$newx
      y = na2$newy
      z = na2$newctrl
      if (verbo) {
        if (i > 2) 
          print(c("i=", i, "non-missing y=", length(y)), 
                quote = FALSE)
      }
      if (length(x) < 5) {
        print("available observations<5")
        break
      }
      ii = ii + 1
      if (verbo) 
        print(c("i=", i, "ii=", ii), quote = FALSE)
      mg1=round(mag_ctrl(x=x0,y=y0,ctrl=z0), dig)
      out[ii,1] = nam[1]
      out[ii,2] = nam[i]
      out[ii,3] = mg1[1]
      out[ii,4] = mg1[2]
    }
    print(out)
    return(out)
  }

