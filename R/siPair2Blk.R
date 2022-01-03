#' Block Version of silentPair2 for causality scores with control variables
#' 
#'
#' Block version allows a new bandwidth (chosen by the np package)
#' while fitting kernel regressions for each block of data. This may
#' not be appropriate in all situations.  Block size is flexible. 
#' The function develops a unanimity index regarding the particular
#' flip (y on xi) or (xi on y) is best. Relevant signs determine the
#' causal direction and unanimity index among three criteria.
#' The `2' in the name of the function suggests a second implementation
#' where exact stochastic dominance, decileVote, and momentVote are used.
#' It avoids Anderson's trapezoidal approximation.
#' The summary results for all
#' three criteria are reported in a vector of numbers 
#' internally called \code{crall}.
#'   
#' @param mtx {The data matrix with p columns. Denote x1 as the first column, 
#' which is fixed and then 
#'  paired with all other columns, say: x2, x3, .., xp, one by one  
#'  flipping with x1.The number of columns, p, must be 2 or more}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths.
#'   The default ctrl=0 means that there are no control variables used.}
#' @param dig {Number of digits for reporting (default \code{dig}=6).}
#' @param blksiz {block size, default=10, if chosen blksiz >n, where n=rows in the
#'  matrix, then blksiz=n. That is, no blocking is done}
#' @return With p columns in \code{mtx} argument to this function, x1 can be 
#' paired with a total of p-1 columns (x2, x3, .., xp). Note
#' we never flip any of the control variables with x1.  This function
#' produces i=1,2,..,p-1 numbers representing the summary sign, or `sum' from 
#' the signs sg1 to sg3 associated with the three criteria:
#' Cr1, Cr2 and Cr3.  Note that sg1 and sg2 themselves are weighted signs using
#' the weighted sum of signs from four orders of stochastic dominance.
#' In general, a positive sign in the i-th location of the `sum' output of this function
#' means that x1 is the kernel cause while the variable in (i+1)-th column of \code{mtx} is the
#' `effect' or `response' or `endogenous.' The magnitude represents the strength (unanimity)
#' of the evidence for a particular sign. Conversely, a negative sign
#'  in the i-th location of the `sum' output of this function means that
#' that the first variable listed as the input to this function is the `effect,'
#' while the variable in (i+1)-th column of \code{mtx} is the exogenous kernel cause.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{bootPairs}}, \code{\link{silentMtx}}
#' @seealso See  \code{\link{someCPairs}}, \code{\link{compPortfo}}
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
#' @note The European Crime data has all three criteria correctly suggesting that
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
#' siPair2Blk(mtcars[,1:3],ctrl=mtcars[,4:5]) # mpg paired with others
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
#'siPair2Blk(mtx=cbind(x2,y2), ctrl=cbind(z,w2))
### }
#' 
#' 
#' @export

siPair2Blk = function(mtx, ctrl = 0, dig = 6, blksiz=10) {
  len = length(ctrl)
  n = NROW(mtx)
  p = NCOL(mtx)
  if (blksiz>n) blksiz=n
  if (p < 2){ 
    stop("too few columns in mtx input to siPairsBlk")}
  npair = p - 1
  cr1 = rep(NA, npair)
  cr2 = rep(NA, npair)
  cr3 = rep(NA, npair)
  crall = rep(NA, npair)
  for (typ in 1:3) {
    # outcause = matrix(NA, nrow = npair, ncol = 7) ii = 0
    for (i in 2:p) {  #i is for column numbers not rows
      x0 = mtx[, i]
      y0 = mtx[, 1] #first column has target y
      if (len > 1) { #ctrl is present case
        na2 = naTriplet(x0, y0, ctrl=ctrl)#matched delete NAs
        x = na2$newx
        y = na2$newy
        z = na2$newctrl
      } #endif len>1
      if (len == 1) { #ctrl is absent case
        na2 = napair(x0, y0)
        x = na2$newx
        y = na2$newy
        z = ctrl
      } #endif len==1 or absent ctrl case
      if (length(x) < 5) {
        print("available observations<5")
        break
      }
  
      im1 = i - 1  #im1 denotes i minus 1
      if (len > 1) { #std=standarzes, rhs=RHS variable, times er=error
        if (typ == 1) # C=control var present
          arxy = absBstdrhserC(x, y, ctrl=z, blksiz=blksiz)#block version 
        if (typ == 1) 
          aryx = absBstdrhserC(y, x, ctrl=z,blksiz=blksiz) #flipped model
        if (typ == 2) 
          arxy = absBstdresC(x, y, ctrl=z,blksiz=blksiz) #std=standardized
        if (typ == 2) #res=residuals C=control var present
          aryx = absBstdresC(y, x, ctrl=z,blksiz=blksiz) #block version
        if (typ < 3) {
          av.crit4 = mean(compPortfo(arxy, aryx))
          if (typ == 1) {
            cr1[im1] = av.crit4
          }
          if (typ == 2) {
            cr2[im1] = av.crit4
          }
        }
        if (typ == 3) {  #Cr3 here
          par1 = parcorBijk(x, y, z)
          rxy=par1$ouij
          ryx=par1$ouji
          del = rxy^2 - ryx^2
          #            print(c('delta',del),q=FALSE)
          cr3[im1] = as.numeric(sign(del))
        }
      }
      
      if (len == 1) {  #len=1 means no control variables
        if (typ == 1) #typ 1 is criterion 1
          arxy = absBstdrhserC(x, y,blksiz=blksiz,ctrl=ctrl)
        if (typ == 1) 
          aryx = absBstdrhserC(y, x,blksiz=blksiz,ctrl=ctrl)#flipped model here
        if (typ == 2) 
          arxy = absBstdresC(x, y,blksiz=blksiz,ctrl=ctrl)
        if (typ == 2) 
          aryx = absBstdresC(y, x,blksiz=blksiz,ctrl=ctrl)
        if (typ < 3) {
          av.crit4 = mean(compPortfo(arxy, aryx))
          
          if (typ == 1) {
            cr1[im1] = av.crit4
          }
          if (typ == 2) {
            cr2[im1] = av.crit4
          }
        }
        if (typ == 3) { #third criterion Cr3 here
          gmc0 = gmcmtxBlk(cbind(x, y),blksiz=blksiz)  #block version
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



