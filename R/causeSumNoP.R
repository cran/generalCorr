#' No print (NoP) version of causeSummBlk summary causal paths from three criteria 
#' 
#' Allowing input matrix of control variables, this function produces 
#' a 5 column matrix
#' summarizing the results where the estimated signs of
#' stochastic dominance order values, (+1, 0, -1), are weighted by 
#'  \code{wt=c(1.2,1.1, 1.05, 1)} to
#' compute an overall result for all orders of stochastic dominance by 
#' a weighted sum for
#' the criteria Cr1 and Cr2 and added to the Cr3 estimate as: (+1, 0, -1).
#' The final range for the unanimity of sign index is [--100, 100].
#' 
#' The reason for slightly declining weights on the signs from
#' SD1 to SD4 is simply that the higher order stochastic dominance
#' numbers are less reliable.
#' The summary results for all
#' three criteria are reported in one matrix called \code{out} but not printed: 
#'   
#' @param mtx {The data matrix with many columns, y the first column 
#' is a fixed target and then it is
#'  paired with all other columns, one by one, and still called x for the 
#'  purpose of flipping.}
#' @param blksiz {block size, default=10, if chosen blksiz >n, where n=rows in matrix
#'      then blksiz=n. That is, no blocking is done}
#' @param nam {vector of column names for \code{mtx}. Default: colnames(mtx)}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths}
#' @param dig {Number of digits for reporting (default \code{dig}=6).}
#' @param wt {Allows user to choose a vector of four alternative weights for SD1 to SD4.}
#' @param sumwt { Sum of weights can be changed here =4(default).}
#' @return If there are p columns in the input matrix, x1, x2, .., xp, say,
#' and if we keep x1 as a common member of all causal-direction-pairs
#' (x1, x(1+j)) for (j=1, 2, .., p-1) which can be flipped. That is, either x1 is
#' the cause or x(1+j) is the cause in a chosen pair.
#' The control
#' variables are not flipped. The printed output of this function
#' reports the results for p-1 pairs indicating which variable
#' (by name) causes which other variable (also by name).
#' It also prints strength or signed summary strength index in range [-100,100]. 
#' A positive sign of the strength index means x1 kernel causes x(1+j),
#' whereas negative strength index means x(1+j) kernel causes x1. The function 
#' also prints Pearson correlation and its p-value. This function also returns
#' a matrix of p-1 rows and 5 columns entitled: 
#' ``cause", ``response", ``strength", ``corr." and ``p-value", respectively
#' with self-explanatory titles. The first two columns have names of variables
#' x1 or x(1+j), depending on which is the cause. The `strength' column
#' has absolute value of summary index in range [0,100]  
#' providing summary of causal results
#' based on preponderance of evidence from Cr1 to Cr3 
#' from four orders of stochastic dominance, etc.  The order of input columns matters.
#' The fourth column `corr.' reports the Pearson correlation coefficient while
#' the fifth column has the p-value for testing the null of zero Pearson coeff.
#' This function calls  \code{siPairsBlk} allowing for control variables.
#' The output of this function can be sent to `xtable' for a nice Latex table. 
#' @importFrom stats complete.cases
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{bootPairs}},  \code{\link{causeSummary0}} has
#' an older version of this function.
#' @seealso See  \code{\link{causeAllPair}} 
#' @seealso \code{\link{siPairsBlk}}, \code{\link{causeSummary}} 
#' @references Vinod, H. D. `Generalized Correlation and Kernel Causality with
#'    Applications in Development Economics' in Communications in
#'    Statistics -Simulation and Computation, 2015,
#'    \doi{gffn86}
#' @references Vinod, H. D. 'New exogeneity tests and causal paths,'
#'  Chapter 2 in 'Handbook of Statistics: Conceptual Econometrics 
#' Using R', Vol.32, co-editors: H. D. Vinod and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2019, pp. 33-64.
#'  
#' @references Vinod, H. D. Causal Paths and Exogeneity Tests 
#' in {Generalcorr} Package for Air Pollution and Monetary Policy 
#' (June 6, 2017). Available at SSRN: 
#' \url{https://www.ssrn.com/abstract=2982128}    
#' @concept  causal path 
#' @concept stochastic dominance orders
#' @concept summary index
#' @note The European Crime data has all three criteria correctly suggesting that
#' high crime rate kernel causes the deployment of a large number of police officers.
#' Since Cr1 to Cr3 near unanimously suggest `crim' as the cause of `off', 
#' strength index 100 suggests unanimity. 
#' \code{attach(EuroCrime); causeSummary(cbind(crim,off))}
#' 
#' @examples
#'
#'
#' \dontrun{
#' mtx=data.frame(mtcars[,1:3])
#' ctrl=data.frame(mtcars[,4:5])
#'  causeSumNoP(mtx=mtx,ctrl=ctrl)
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
#'causeSumNoP(mtx=cbind(x2,y2), ctrl=cbind(z,w2))
### }
#' 
#' 
#' @export


causeSumNoP=   #noprint version
  function(mtx, nam = colnames(mtx), blksiz=10,
           ctrl = 0, dig = 6, wt = c(1.2, 1.1, 1.05, 1), sumwt = 4)
  {
# require(generalCorr); require(PerformanceAnalytics); options(np.messages=FALSE)
    p = NCOL(mtx)
    if (p < 2) 
      stop("stop:too few columns in input mtx to siPairsBlk")
    #Task 1 compute corr and p-val with irwise deletion of NAs
    pv=rep(NA,p-1)
    pearson=rep(NA,p-1)
    for ( i in 2:p){
      x=mtx[,1]
      y=mtx[,i]
      ok=complete.cases(x,y)   #non-missing data rows pairwise
      c1 = cor.test(x[ok], y[ok])
      pv[i] = c1$p.value
      pearson[i] = c1$estimate
    }
    si0 = siPairsBlk(mtx, ctrl = ctrl, dig = dig, wt = wt, 
                     blksiz=blksiz, sumwt = sumwt)
    si = round(100 * as.numeric(si0)/3.175, 3)
    out = matrix(NA, nrow = (p - 1), ncol = 5)
    for (i in 2:p) {
      if (si[i - 1] < 0) {
 #print(c(nam[i], "causes", nam[1], "strength=", si[i - 1]), quote = FALSE)
        out[i - 1, 1] = nam[i]
        out[i - 1, 2] = nam[1]
        out[i - 1, 3] = abs(si[i - 1]) #abs strength in out matrix
      }
      if (si[i - 1] > 0) {
  #  print(c(nam[1], "causes", nam[i], "strength=", si[i - 1]), quote = FALSE)
        out[i - 1, 1] = nam[1]
        out[i - 1, 2] = nam[i]
        out[i - 1, 3] = abs(si[i - 1]) #abs value in the out matrix
      }
   #print(c("corr=", round(pearson[i], 4), "p-val=", round(pv[i], 5)), quote = FALSE)
      out[i - 1, 4] = round(pearson[i], 4)
      out[i - 1, 5] = round(pv[i], 5)
    }  #end i loop 
    colnames(out) = c("cause","response","strength","corr.","p-value")
    return(out)
  }
