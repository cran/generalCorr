#' kernel causality (version 2) scores with control variables 
#' 
#' This function uses flipped kernel regressions to decide causal directions. 
#' This version 2 avoids Anderson's trapezoidal approximation used in
#' `silenPairs.' It calls functions: decileVote, momentVote, exactSdMtx, 
#' and summaryRank
#' after stochastic dominance is computed. It computes 
#' an average of ranks used. The column with the ``choice'' rank
#'  value helps in choosing the flip having the lowest 
#'  Hausman-Wu (residual times RHS regressor)
#'  and secondly the lowest absolute 
#'  residual. The chosen  flipped regression defines
#'  the ``cause" based on the variable on its right-hand side. In portfolio
#'  selection, choice rank 1 has the highest return. Here we want low residuals and
#'  low Hausman-Wu value, hence we choose choice=2 as the desirable flip.
#'  
#'  
#' The function develops a unanimity index regarding the particular
#' flip (y on xi) or (xi on y) is best. A summary of all relevant signs determines the
#' causal direction and unanimity index among three criteria.
#' The `2' in the name of the function suggests a second implementation
#' where exact stochastic dominance, decileVote, and momentVote algorithms are used.
#' @param mtx {The data matrix with p columns. Denote x1 as the first column, 
#' which is fixed in all rows of the output and then it is
#'  paired with all other columns, say: x2, x3, .., xp, one by one for the 
#'  purpose of flipping with x1. p must be 2 or more}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths,
#'   default is ctrl=0, which means that there are no control variables used.}
#' @param dig {Number of digits for reporting (default \code{dig}=6).}
#' @return A matrix with p columns in \code{mtx} argument to this function, x1 can be 
#' paired with a total of p-1 columns (x2, x3, .., xp). Note
#' we never flip any of the control variables with x1.  This function
#' produces i=1,2,..,p-1 numbers representing the summary sign, or `sum' from 
#' the signs sg1 to sg3 associated with the three criteria:
#' Cr1, Cr2, and Cr3.  Note that sg1 and sg2 themselves are weighted signs using 
#' a weighted sum of signs from four orders of stochastic dominance.
#' In general, a positive sign in the i-th location of the `sum' output of this function
#' means that x1 is the kernel cause while the variable in (i+1)-th column of \code{mtx} is the
#' `effect' or `response' or `endogenous.' The magnitude represents the strength (unanimity)
#' of the evidence for a particular sign. Conversely, a negative sign
#'  in the i-th location of the `sum' output of this function means 
#' that the first variable listed as the input to this function is the `effect,'
#' while the variable in (i+1)-th column of \code{mtx} is the exogenous kernel cause.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{summaryRank}}, \code{\link{decileVote}}
#' @seealso See  \code{\link{momentVote}}, \code{\link{exactSdMtx}}
#' @references H. D. Vinod  'Generalized Correlation and Kernel Causality with
#'    Applications in Development Economics' in Communications in
#'    Statistics -Simulation and Computation, 2015,
#'    \doi{10.1080/03610918.2015.1122048}
#' 
#' @references Vinod, H. D. Causal Paths and Exogeneity Tests 
#' in {Generalcorr} Package for Air Pollution and Monetary Policy 
#' (June 6, 2017). Available at SSRN: 
#' \url{https://www.ssrn.com/abstract=2982128}  
#' @concept Hausman-Wu exogeneity criteria
#' @concept stochastic dominance
#' @concept generalized correlations
#' @note The European Crime data has all three criteria correctly suggesting that a
#' high crime rate kernel causes the deployment of a large number of police officers.
#' The command \code{attach(EuroCrime); silentPairs(cbind(crim,off))}
#' returns only one number: 3.175, implying the highest unanimity strength index,
#' with the positive sign suggesting `crim' in the first column kernel causes
#' `off' in the second column of the argument \code{mtx} to this function.
#' 
#' @examples
#'
#'
#' \dontrun{
#' options(np.messages=FALSE)
#' colnames(mtcars[2:ncol(mtcars)])
#' silentPair2(mtcars[,1:3],ctrl=mtcars[,4:5]) # mpg paired with others
#' }
#' 
### \dontrun{
#'options(np.messages=FALSE)
#'set.seed(234)
#'z=runif(10,2,11)# z is independently created
#'x=sample(1:10)+z/10 #x is somewhat indep and affected by z
#'y=1+2*x+3*z+rnorm(10)
#'w=runif(10)
#'x2=x;x2[4]=NA;y2=y;y2[8]=NA;w2=w;w2[4]=NA
#'silentPair2(mtx=cbind(x2,y2), ctrl=cbind(z,w2))
### }
#' 
#' 
#' @export

silentPair2 = function(mtx, ctrl = 0, dig = 6) {
  len = length(ctrl)
  n = NROW(mtx)
  p = NCOL(mtx)
  if (p < 2) 
    stop("too few columns in mtx input to silentPairs")
  npair = p - 1
  cr1 = rep(NA, npair)
  cr2 = rep(NA, npair)
  cr3 = rep(NA, npair)
  crall = rep(NA, npair)
  for (typ in 1:3) {
   # outcause = matrix(NA, nrow = npair, ncol = 7) ii = 0
    for (i in 2:p) {
      x0 = mtx[, i]
      y0 = mtx[, 1]
      if (len > 1) {
        z0 = ctrl
        na2 = naTriplet(x0, y0, z0)
        x = na2$newx
        y = na2$newy
        z = na2$newctrl
      }
      if (len == 1) {
        na2 = napair(x0, y0)
        x = na2$newx
        y = na2$newy
      }
      
      
      # if (verbo) { if (i > 2) print(c('i=', i, 'non-missing y=',
      # length(y)), quote = FALSE) }
      if (length(x) < 5) {
        print("available observations<5")
        break
      }
      # ii = ii + 1
      im1 = i - 1
      if (len > 1) {
        if (typ == 1) 
          arxy = abs_stdrhserC(x, y, z)
        if (typ == 1) 
          aryx = abs_stdrhserC(y, x, z) #flipped model
        if (typ == 2) 
          arxy = abs_stdresC(x, y, z)
        if (typ == 2) 
          aryx = abs_stdresC(y, x, z)
        if (typ < 3) {
          av.crit4 =mean( compPortfo(arxy, aryx))
          
          if (typ == 1) {
            cr1[im1] = av.crit4
          }
          if (typ == 2) {
            cr2[im1] = av.crit4
          }
        }
        if (typ == 3) {
          par1 = parcor_ijk(x, y, z)
          rxy=par1$ouij
          ryx=par1$ouji
          del = rxy^2 - ryx^2
          #            print(c('delta',del),q=FALSE)
          cr3[im1] = as.numeric(sign(del))
        }
      }
      
      if (len == 1) {
        if (typ == 1) 
          arxy = abs_stdrhserr(x, y)
        if (typ == 1) 
          aryx = abs_stdrhserr(y, x)#flipped model here
        if (typ == 2) 
          arxy = abs_stdres(x, y)
        if (typ == 2) 
          aryx = abs_stdres(y, x)
        if (typ < 3) {
          av.crit4 =mean(compPortfo(arxy, aryx))
          
          if (typ == 1) {
            cr1[im1] = av.crit4
          }
          if (typ == 2) {
            cr2[im1] = av.crit4
          }
        }
        if (typ == 3) {
          gmc0 = gmcmtx0(cbind(x, y))
          rxy = gmc0[1, 2]
          ryx = gmc0[2, 1]
          del = rxy^2 - ryx^2
          #            print(c('delta',del),q=FALSE)
          cr3[im1] = as.numeric(sign(del))
        }
      }
      
    }  #end i loop
  }  #end typ loop
  
  
  for (j in 1:npair) {
    cr13 = c(cr1[j], cr2[j], cr3[j])
    crall[j] = round(sum(cr13, na.rm = TRUE), dig)
  }
  return(crall)
}



