#' Block Version 2: Kernel causality summary of causal paths from three criteria 
#' 
#' A block version of \code{causeSummary()} chooses new bandwidth for every
#' ten (blksiz=10) observations chosen by the `np' package injecting flexibility.
#' While allowing the researcher to keep some variables as controls,
#' or outside the scope of causal path determination 
#' (e.g., age or latitude), this function produces detailed causal path information. 
#' The output table is a 5-column matrix identifying the names of variables,
#' causal path directions, and path strengths re-scaled to be in the 
#' range [--100, 100], (table reports absolute values of the strength)
#' plus Pearson correlation coefficient and its p-value.
#' 
#' The algorithm determines causal path directions from the sign
#' of the strength index. The strength index magnitudes are computed by comparing 
#' three aspects of flipped kernel regressions: 
#' [x1 on (x2, x3, .. xp)] and its flipped version [x2 on (x1, x3, .. xp)].
#' The cause should be on the right-hand side of
#' the regression equation. The properties
#' of regression fit determine which flip is superior.
#' We compare (Cr1) formal exogeneity test criterion, (residuals times RHS
#' regressor, where smaller in absolute value is better) 
#' (Cr2) absolute values of residuals, where 
#' smaller in absolute value is better, and
#' (Cr3) R-squares of the flipped regressions implying three criteria Cr1, to Cr3.
#' The criteria are quantified by sophisticated methods using four orders
#' of stochastic dominance, SD1 to SD4. We assume slightly declining weights on 
#' the sign observed by Cr1 to Cr3. The user can change default weights.
#'   
#' @param mtx {The data matrix with many columns, y the first column 
#' is a fixed target, and then it is
#'  paired with all other columns, one by one, and still called x for the 
#'  purpose of flipping.}
#' @param blksiz {block size, default=10, if chosen blksiz >n, where n=rows in the matrix
#'      then blksiz=n. That is, no blocking is done}
#' @param nam {vector of column names for \code{mtx}. Default: colnames(mtx)}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths}
#' @param dig {The number of digits for reporting (default \code{dig}=6).}
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
#' It also prints a strength, or signed summary strength 
#' index forced to be in the range [-100,100] for easy interpretation. 
#' A positive sign of the strength index means x1 kernel causes x(1+j),
#' whereas negative strength index means x(1+j) kernel causes x1. The function 
#' also prints Pearson correlation and its p-value. This function also returns
#' a matrix of p-1 rows and 5 columns entitled: 
#' ``cause", ``response", ``strength", ``corr." and ``p-value", respectively
#' with self-explanatory titles. The first two columns have names of variables
#' x1 or x(1+j), depending on which is the cause. The `strength' column
#' has the absolute value of a summary index in the range [0,100],  
#' providing a summary of causal results
#' based on the preponderance of evidence from Cr1 to Cr3 
#' from four orders of stochastic dominance, etc.  The order of input columns matters.
#' The fourth column of the output matrix entitled `corr.' reports the Pearson 
#' correlation coefficient, while
#' the fifth column of the output matrix has the p-value for testing the 
#' null hypothesis of a zero Pearson coefficient.
#' This function calls  \code{siPairsBlk}, allowing for control variables.
#' The output of this function can be sent to `xtable' for a nice Latex table. 
#' @importFrom stats complete.cases
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{bootPairs}},  \code{\link{causeSummary}} has
#' an older version of this function.
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
#' a high crime rate kernel causes the deployment of a large number of police officers.
#' Since Cr1 to Cr3 near-unanimously suggest `crim' as the cause of `off', 
#' a strength index of 100 suggests unanimity. 
#' \code{attach(EuroCrime); causeSummBlk(cbind(crim,off))}
#' 
#' @examples
#'
#'
#' \dontrun{
#' mtx=as.matrix(mtcars[,1:3])
#' ctrl=as.matrix(mtcars[,4:5])
#'  causeSummBlk(mtx,ctrl,nam=colnames(mtx))
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
#'causeSummBlk(mtx=cbind(x2,y2), ctrl=cbind(z,w2))
### }
#' 
#' 
#' @export

causeSummBlk = function(mtx, nam = colnames(mtx), blksiz=10,
     ctrl=0, dig = 6, wt = c(1.2, 1.1, 1.05, 1), sumwt = 4)
   {
    p = NCOL(mtx)
    n = NROW(mtx)
    if (blksiz>n) blksiz=n
    if (p < 2) 
        stop("too few columns in input to summaryCause mtx")
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
    si0 = siPairsBlk(mtx, ctrl=ctrl, dig = dig, wt=wt, 
                     sumwt = 4, blksiz=blksiz) #block version
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
