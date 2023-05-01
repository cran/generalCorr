#' All Pair Version Kernel (block) causality summary paths from three criteria 
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
#' SD1 to SD4 stochastic dominance orders is simply their
#' slightly increasing sampling
#' unreliability due to higher order trapezoidal approximations of
#' integrals of densities involved in definitions of SD1 to SD4.
#' The summary results for all
#' three criteria are reported in one matrix called \code{out}: 
#'   
#' @param mtx {The data matrix with many columns, We consider causal paths
#' among all
#' possible pairs of mtx columns.}
#' @param blksiz {block size, default=10, if chosen blksiz >n, where n=rows in matrix
#'      then blksiz=n. That is, no blocking is done}
#' @param nam {vector of column names for \code{mtx}. Default: colnames(mtx)}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths}
#' @param dig {Number of digits for reporting (default \code{dig}=6).}
#' @param wt {Allows user to choose a vector of four alternative weights for SD1 to SD4.}
#' @param sumwt { Sum of weights can be changed here =4(default).}
#' @return If there are p columns in the input matrix, x1, x2, .., xp, say,
#' there are choose(p,2) or  [p*(p-1)/2] possible pairs and as many causal paths.
#' This function returns
#' a matrix of p*(p-1)/2 rows and 5 columns entitled: 
#' ``cause", ``response", ``strength", ``corr." and ``p-value", respectively
#' with self-explanatory titles. The first two columns have names of variables
#' x1 or x(1+j), depending on which is the cause. The `strength' column
#' has absolute value of summary index in range [0,100]  
#' providing summary of causal results
#' based on preponderance of evidence from criteria  Cr1 to Cr3 
#' from four orders of stochastic dominance, etc.  
#' The fourth column `corr.' reports the Pearson correlation coefficient while
#' the fifth column has the p-value for testing the null of zero Pearson coeff.
#' This function merely calls \code{causeSumNoP} repeatedly to include all pairs.
#' The background function \code{siPairsBlk} allows for control variables.
#' The output of this function can be sent to `xtable' for a nice Latex table. 
#' @importFrom stats complete.cases
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{bootPairs}},  \code{\link{causeSummBlk}} 
#' @seealso See  \code{\link{someCPairs}} 
#' @seealso \code{\link{siPairsBlk}}, \code{\link{causeSummary}} 
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
#' mtx=data.frame(mtcars[,1:3]) #make sure columns of mtx have names
#' ctrl=data.frame(mtcars[,4:5])
#'  causeAllPair(mtx=mtx,ctrl=ctrl)
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
#'causeAllPair(mtx=cbind(x2,y2), ctrl=cbind(z,w2))
### }
#' 
#' 
#' @export



causeAllPair=
  function(mtx, nam = colnames(mtx), blksiz=10,
           ctrl = 0, dig = 6, wt = c(1.2, 1.1, 1.05, 1), sumwt = 4)  {
# require(generalCorr); require(PerformanceAnalytics); options(np.messages=FALSE)
    p = NCOL(mtx)
    if (p < 2)  stop("stop:too few columns for siPairsBlk mtx")
    pm1=p-1
    mtx1=mtx[,1:p]
    out=causeSumNoP(mtx=mtx1, nam =nam, blksiz=blksiz,
              ctrl = ctrl, dig = dig, wt = wt, sumwt = sumwt)
   if (p>2) {
    for(j in 2:pm1) {
     mtx1=mtx[,j:p]   
  ou=causeSumNoP(mtx=mtx1, nam =colnames(mtx1), blksiz=blksiz,
              ctrl = ctrl, dig = dig, wt = wt, sumwt = sumwt)
 #  print(ou)
    out=rbind(out,ou)
    }#end j loop
   }# end if p>2 condition
    return(out)
  }
