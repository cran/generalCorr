#' Block Version 2: Kernel causality summary of causal paths from three criteria 
#' 
#' The `2' in the name of the function suggests a second implementation
#' where exact stochastic dominance, `decileVote' and `momentVote' functions are used,
#' Block version allows a new bandwidth (chosen by the np package)
#' while fitting kernel regressions for each block of data. This may
#' not be appropriate in all situations.  Block size is flexible. 
#' The function develops a unanimity index regarding which regression
#' flip, (y on xi) or (xi on y) is the best. The ``cause'' is 
#' always on the right-hand side of a regression equation, and
#' the superior flip gives the correct sign. The summary of all signs determines the
#' causal direction and unanimity index among three criteria. This is
#' a block version of \code{causeSummary2()}.
#' While allowing the researcher to keep some variables as controls,
#' or outside the scope of causal path determination 
#' (e.g., age or latitude)  this function produces detailed causal path information 
#' in a 5 column matrix identifying the names of variables,
#' causal path directions, path strengths re-scaled to be in the 
#' range [--100, 100], (table reports absolute values of the strength)
#' plus Pearson correlation and its p-value.
#'  
#' The algorithm determines causal path directions from the sign
#' of the strength index and strength index values by comparing 
#' three aspects of flipped kernel regressions: 
#' [x1 on (x2, x3, .. xp)] and its flipped version [x2 on (x1, x3, .. xp)]
#' We compare (i) formal exogeneity test criterion, (ii) absolute residuals, and
#' (iii) R-squares of the flipped regressions implying three criteria Cr1, to Cr3.
#' The criteria are quantified by new methods using four orders
#' of stochastic dominance, SD1 to SD4. See Vinod (2021) two SSRN papers.
#'   
#' @param mtx {The data matrix with many columns, y the first column 
#' is a fixed target, and then it is
#'  paired with all other columns, one by one, and still called x for
#'  flipping.}
#' @param blksiz {block size, default=10, if chosen blksiz >n, where n=rows 
#' in the matrix then blksiz=n. That is, no blocking is done}
#' @param nam {vector of column names for \code{mtx}. Default: colnames(mtx)}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths}
#' @param dig {The number of digits for reporting (default \code{dig}=6).}
#' @return If there are p columns in the input matrix, x1, x2, .., xp, say,
#' and if we keep x1 as a common member of all causal-direction-pairs
#' (x1, x(1+j)) for (j=1, 2, .., p-1) which can be flipped. That is, either x1 is
#' the cause or x(1+j) is the cause in a chosen pair.
#' The control
#' variables are not flipped. The printed output of this function
#' reports the results for p-1 pairs indicating which variable
#' (by name) causes which other variable (also by name).
#' It also prints the strength or signed summary strength index in 
#' the range [-100,100]. 
#' A positive sign of the strength index means x1 kernel causes x(1+j),
#' whereas negative strength index means x(1+j) kernel causes x1. The function 
#' also prints Pearson correlation and its p-value. This function also returns
#' a matrix of p-1 rows and 5 columns entitled: 
#' ``cause", ``response", ``strength", ``corr." and ``p-value", respectively
#' with self-explanatory titles. The first two columns have names of variables
#' x1 or x(1+j), depending on which is the cause. The `strength' column
#' has an absolute value of the summary index in the range [0,100],  
#' providing a summary of causal results
#' based on the preponderance of evidence from Cr1 to Cr3 from deciles, moments,
#' from four orders of stochastic dominance.  
#' The order of input columns in "mtx" matters.
#' The fourth column, `corr.', reports the Pearson correlation coefficient, while
#' the fifth column has the p-value for testing the null of zero Pearson coefficient.
#' This function calls  \code{siPairsBlk}, allowing for control variables.
#' The output of this function can be sent to `xtable' for a nice Latex table. 
#' @importFrom stats complete.cases
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{bootPairs}},  \code{\link{causeSummary}} has
#' an older version of this function.
#' @seealso See  \code{\link{someCPairs}} 
#' @seealso \code{\link{siPair2Blk}}, \code{\link{causeSummary2}} 
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
#'  
#' @references Vinod, Hrishikesh D., R Package GeneralCorr 
#' Functions for Portfolio Choice 
#' (November 11, 2021). Available at SSRN: 
#' https://ssrn.com/abstract=3961683 
#' 
#' @references Vinod, Hrishikesh D., Stochastic Dominance 
#' Without Tears (January 26, 2021). Available at 
#' SSRN: https://ssrn.com/abstract=3773309 
#' 
#' @concept  causal path 
#' @concept stochastic dominance orders
#' @concept summary index
#' @note The European Crime data has all three criteria correctly suggesting that
#' high crime rate kernel causes the deployment of a large number of police officers.
#' If Cr1 to Cr3 near-unanimously suggest `crim' as the cause of `off', 
#' strength index would be near 100 suggesting unanimity. 
#' \code{attach(EuroCrime); causeSum2Blk(cbind(crim,off))}
#' 
#' @examples
#'
#' \dontrun{
#' mtx=as.matrix(mtcars[,1:3])
#' ctrl=as.matrix(mtcars[,4:5])
#' causeSum2Blk(mtx,ctrl,nam=colnames(mtx))
#' }
#' 
#' 
#' @export

causeSum2Blk = function(mtx, nam = colnames(mtx), blksiz=10,
     ctrl=0, dig = 6)
   {
    p = NCOL(mtx)
    n = NROW(mtx)
    if (blksiz>n) blksiz=n
    if (p < 2) 
        stop("too few columns in input to causeSum2Blk mtx")
#Task 1 compute corr and p-val with pairwise deletion of NAs
   pv=rep(NA,p-1) #place to hold p-values
   pearson=rep(NA,p-1)#place to hold corr.coeffs
   for ( i in 2:p){
   x=mtx[,1]
   y=mtx[,i]
   ok=complete.cases(x,y)   #non-missing data rows pairwise
   c1 = cor.test(x[ok], y[ok])
   pv[i] = c1$p.value
   pearson[i] = c1$estimate
   }
    si0 = siPair2Blk(mtx, ctrl=ctrl, dig = dig, 
     blksiz=blksiz) #block version
    si = round(100 * as.numeric(si0)/3.175, 3)
    out = matrix(NA, nrow = (p - 1), ncol = 5)
    for (i in 2:p) {
        if (si[i - 1] < 0) {
            print(c(nam[i], "causes", nam[1], "strength=", si[i - 1]), quote = FALSE)
            out[i - 1, 1] = nam[i]
            out[i - 1, 2] = nam[1]
            out[i - 1, 3] = abs(si[i - 1]) #abs strength in out matrix
        }
        if (si[i - 1] > 0) {
            print(c(nam[1], "causes", nam[i], "strength=", si[i - 1]), quote = FALSE)
            out[i - 1, 1] = nam[1]
            out[i - 1, 2] = nam[i]
            out[i - 1, 3] = abs(si[i - 1]) #abs value in the out matrix
        }
        print(c("corr=", round(pearson[i], 4), "p-val=", round(pv[i], 5)), quote = FALSE)
        out[i - 1, 4] = round(pearson[i], 4)
        out[i - 1, 5] = round(pv[i], 5)
    }  #end i loop 
    colnames(out) = c("cause","response","strength","corr.","p-value")
    return(out)
}  #end function
