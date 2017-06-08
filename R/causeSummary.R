#' Function for kernel causality summary of results 
#' 
#' Allowing input matrix of control variables, this function produces a 5 column matrix
#' summarizing the results where the estimated signs of
#' stochastic dominance order values (+1, 0, -1) are weighted by 
#'  \code{wt=c(1.2,1.1, 1.05, 1)} to
#' compute an overall result for all orders of stochastic dominance by a weighted sum for
#' the crietria Cr1 and Cr2 and added to the Cr3 estimate as: (+1, 0, -1).
#' 
#' The reason for slightly declining weights on the signs from
#' SD1 to SD4 is simply that the local mean comparisons 
#' implicit in SD1 are known to be
#' more reliable than local variance implicit in SD2, local skewness implicit in
#' SD3 and local kurtosis implicit in SD4. The reason for slightly declining sampling
#' unreliability of higher moments is simply that SD4 involves fourth power
#' of the deviations from the mean and SD3 involves 3rd power, etc.
#' The summary results for all
#' three criteria are reported in one matrix called \code{out}: 
#'   
#' @param mtx {The data matrix with many columns, y the first column is fixed and then 
#'  paired with all columns, one by one, and still called x for the purpose of flipping.}
#' @param nam {vector of column names for \code{mtx}. Default: colnames(mtx)}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths}
#' @param dig {Number of digits for reporting (default \code{dig}=6).}
#' @param wt {Allows user to choose a vector of four alternative weights for SD1 to SD4.}
#' @param sumwt { Sum of weights can be changed here =4(default).}
#' @return If there are p columns in the input matrix, there are p-1 possible
#' flipped X, Y pairs leaving the control variables alone.  
#' This function prints the results for p-1 pairs indicating which causes what
#' it also prints strength or signed summary index in range [-100,100] 
#' In general, a positive sign means first input variable causes the second, etc
#' that the first variable listed as the input to this function is the `kernel cause,' 
#' whereas negative strength index means second column variable causes the first.
#' Also, it prints Pearson correlation and its p-value. This function also returns
#' a matrix of p-1 rows and 5 columns providing summary of causal results.
#' the first column names the causal variable and second names the response.
#' the third column has absolute value of summary index in range [0,100] 
#' summarizing preponderance of evidence from Cr1 to Cr3 
#' from four orders of stochastic dominance, etc.  The order of input columns matters.
#' the fourth column of the output matrix has Pearson correlation coefficient
#' the fifth column has the p-value for testing the null of zero Pearson coeff.
#' Suggested column headings are Cause, Response, Strength, r, p-value
#' This function calls  \code{silentPairs} allowing for control variables. 
#' @importFrom xtable xtable
#' @importFrom stats complete.cases
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{bootPairs}}
#' @seealso See  \code{\link{someCPairs}} 
#' @seealso \code{\link{silentPairs}}
#' @references Vinod, H. D. 'Generalized Correlation and Kernel Causality with
#'    Applications in Development Economics' in Communications in
#'    Statistics -Simulation and Computation, 2015,
#'    \url{http://dx.doi.org/10.1080/03610918.2015.1122048}
#' @keywords causal path, SD1, SD2, SD3, SD4, summary index
#' @note The European Crime data has all three criteria correctly suggesting that
#' high crime rate kernel causes the deployment of a large number of police officers.
#' Cr1 to Cr3 unanimously suggest `crim' as the cause of `off', so index is 100.
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
    # require(generalCorr); require(PerformanceAnalytics); options(np.messages=FALSE)
    p = NCOL(mtx)
    if (p < 2) 
        stop("too few columns in input to summaryCause mtx")
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
