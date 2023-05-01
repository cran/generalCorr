#' Kernel causality summary of evidence for causal paths from three criteria 
#' 
#' While allowing the researcher to keep some variables as controls,
#' or outside the scope of causal path determination 
#' (e.g., age or latitude)  this function produces detailed causal path information 
#' in a 5 column matrix identifying the names of variables,
#' causal path directions, path strengths re-scaled to be in the 
#' range [--100, 100], (table reports absolute values of the strength)
#' plus Pearson correlation and its p-value.
#' 
#' 
#' The algorithm determines causal path directions from the sign
#' of the strength index and strength index values by comparing 
#' three aspects of flipped kernel regressions: 
#' [x1 on (x2, x3, .. xp)] and its flipped version [x2 on (x1, x3, .. xp)]
#' We compare (i) formal exogeneity test criterion, (ii) absolute residuals, and
#' (iii) R-squares of the flipped regressions implying three criteria Cr1, to Cr3.
#' The criteria are quantified by sophisticated methods using four orders
#' of stochastic dominance, SD1 to SD4. We assume slightly declining weights on 
#' causal path signs because known reliability ranking. SD1 is better than SD2,
#' better than SD3, better than SD4. The user can optionally change our weights.
#'   
#' @param mtx {The data matrix with many columns, y the first column 
#' is fixed and then 
#'  paired with all columns, one by one, and still called x for the 
#'  purpose of flipping.}
#' @param nam {vector of column names for \code{mtx}. Default: colnames(mtx)}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths}
#' @param dig {Number of digits for reporting (default \code{dig}=6).}
#' @param wt {Allows user to choose a vector of four alternative weights for SD1 to SD4.}
#' @param sumwt { Sum of weights can be changed here =4(default).}
#' @return If there are p columns in the input matrix, x1, x2, .., xp, say,
#' and if we keep x1 as a common member of all causal direction pairs
#' (x1, x(1+j)) for (j=1, 2, .., p-1) which can be flipped. That is, either x1 is
#' the cause or x(1+j) is the cause in a chosen pair.
#' The control
#' variables are not flipped. The printed output of this function
#' reports the results for p-1 pairs indicating which variable
#' (by name) causes which another variable (also by name).
#' It also prints a signed summary strength index in the range [-100,100]. 
#' A positive sign of the strength index means x1 kernel causes x(1+j),
#' whereas negative strength index means x(1+j) kernel causes x1. The function 
#' also prints Pearson correlation and its p-value. In short, function returns
#' a matrix of p-1 rows and 5 columns entitled: 
#' ``cause", ``response", ``strength", ``corr." and ``p-value", respectively
#' with self-explanatory titles. The first two columns have names of variables
#' x1 or x(1+j), depending on which is the cause. The `strength' column
#' reports the absolute value of summary index, now in the range [0,100]  
#' providing summary of causal results
#' based on preponderance of evidence from Cr1 to Cr3 
#' from four orders of stochastic dominance, etc.  The order of input columns matters.
#' The fourth column `corr.' reports the Pearson correlation coefficient while
#' the fifth column has the p-value for testing the null of zero Pearson coeff.
#' This function calls  \code{silentPairs} allowing for control variables.
#' The output of this function can be sent to `xtable' for a nice Latex table. 
#' @importFrom stats complete.cases
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{bootPairs}},  \code{\link{causeSummary0}} has
#' an older version of this function.
#' @seealso See  \code{\link{someCPairs}} 
#' @seealso \code{\link{silentPairs}}
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
#' strength index 100 suggests unanimity. In portfolio
#' applications of stochastic dominance one wants higher returns.  Here we are
#' comparing two probability distributions of absolute residuals for two
#' flipped models. We choose that flip which has smaller absolute residuals
#' or better fit.
#' \code{attach(EuroCrime); causeSummary(cbind(crim,off))}
#' 
#' @examples
#'
#'
#' \dontrun{
#' mtx=as.matrix(mtcars[,1:3])
#' ctrl=as.matrix(mtcars[,4:5])
#'  causeSummary(mtx,ctrl,nam=colnames(mtx))
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
#'causeSummary(mtx=cbind(x2,y2), ctrl=cbind(z,w2))
### }
#' 
#' 
#' @export

causeSummary = function(mtx, nam = colnames(mtx), 
     ctrl = 0, dig = 6, wt = c(1.2, 1.1, 1.05, 1), sumwt = 4)
   {
    p = NCOL(mtx)
    if (p < 2) 
        stop("too few columns in input to summaryCause mtx")
#Task 1 compute corr and p-val with pairwise deletion of NAs
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
    si0 = silentPairs(mtx, ctrl = ctrl, dig = dig, wt = wt, sumwt = 4)
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
