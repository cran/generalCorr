#' No-print kernel-causality unanimity score matrix with optional control variables
#' 
#' Allowing input matrix of control variables and missing data, this function produces a
#' p by p matrix summarizing the results, where the estimated signs of
#' stochastic dominance order values (+1, 0, --1) are weighted by 
#' \code{wt=c(1.2, 1.1, 1.05, 1)} to
#' compute an overall result for all orders of stochastic dominance by a weighted sum for
#' the criteria Cr1 and Cr2 and added to the Cr3 estimate as: (+1, 0, --1).
#' Final weighted index is always in the range [--3.175, 3.175]. It is converted
#' to the more intuitive range [--100, 100].
#' 
#' The reason for slightly declining weights on the signs from
#' SD1 to SD4 is simply that the local mean comparisons 
#' implicit in SD1 are known to be
#' more reliable than local variance implicit in SD2, local skewness implicit in
#' SD3 and local kurtosis implicit in SD4. Why are higher moment
#' estimates less reliable? The 
#' higher power of the deviations from the mean needed in their computations
#' lead to greater sampling variability.
#' The summary results for all
#' three criteria are reported in a vector of numbers internally called \code{crall}: 
#'   
#' @param mtx {The data matrix with p columns. Denote x1 as the first column 
#' which is fixed and then 
#'  paired with all other columns, say: x2, x3, .., xp, one by one for the 
#'  purpose of flipping with x1. p must be 2 or more}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths}
#' @param dig {Number of digits for reporting (default \code{dig}=6).}
#' @param wt {Allows user to choose a vector of four alternative weights for SD1 to SD4.}
#' @param sumwt { Sum of weights can be changed here =4(default).}
#' @return With p columns in \code{mtx} argument to this function, x1 can be 
#' paired with a total of p-1 columns (x2, x3, .., xp). Note
#' we never flip any of the control variables with x1.  This function
#' produces i=1,2,..,p-1 numbers representing the summary sign, or `sum' from 
#' the signs sg1 to sg3 associated with the three criteria:
#' Cr1, Cr2 and Cr3.  Note that sg1 and sg2 themselves are weighted signs using
#' weighted sum of signs from four orders of stochastic dominance.
#' In general, a positive sign in the i-th location of the `sum' output of this function
#' means that x1 is the kernel cause while the variable in (i+1)-th column of \code{mtx} is the
#' `effect' or `response' or `endogenous.' The magnitude represents the strength (unanimity)
#' of the evidence for a particular sign. Conversely a negative sign
#'  in the i-th location of the `sum' output of this function means that
#' that the first variable listed as the input to this function is the `effect,'
#' while the variable in (i+1)-th column of \code{mtx} is the exogenous kernel cause.
#' This function is a summary of \code{someCPairs} 
#' allowing for control variables. 
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{silentPairs}}.
#' @seealso See  \code{\link{someCPairs}}, \code{\link{some0Pairs}}
#' @references H. D. Vinod 'Generalized Correlation and Kernel Causality with
#'    Applications in Development Economics' in Communications in
#'    Statistics -Simulation and Computation, 2015,
#'    \url{http://dx.doi.org/10.1080/03610918.2015.1122048}
#' 
#' @references Vinod, H. D. Causal Paths and Exogeneity Tests 
#' in {Generalcorr} Package for Air Pollution and Monetary Policy 
#' (June 6, 2017). Available at SSRN: 
#' \url{https://www.ssrn.com/abstract=2982128}  
#' @concept causal criteria
#' @concept fourth order stochastic dominance
#' @concept  generalized correlations
#' @note The European Crime data has all three criteria correctly suggesting that
#' high crime rate kernel causes the deployment of a large number of police officers.
#' The command \code{attach(EuroCrime); silentPairs(cbind(crim,off))}
#' returns only one number: 3.175, implying a high unanimity strength.
#' The index 3.175 is the highest.
#' The positive sign of the index suggests that `crim' 
#' variable in the first column of the matrix input to this function kernel causes
#' `off' in the second column of the matrix argument \code{mtx} to this function.
#' @note Interpretation of the output matrix produced by this function is as follows.
#' A negative index means the variable named in the column kernel-causes 
#' the variable named in the row. A
#' positive index means the row name variable kernel-causes 
#' the column name variable. The
#' abs(index) measures unanimity by three criteria, Cr1 to Cr3 representing
#' the strength of evidence for the identified causal path.
#' 
#' @examples
#'
#'
#' \dontrun{
#' options(np.messages=FALSE)
#' colnames(mtcars[2:ncol(mtcars)])
#' silentMtx(mtcars[,1:3],ctrl=mtcars[,4:5]) # mpg paired with others
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
#'silentMtx(mtx=cbind(x2,y2), ctrl=cbind(z,w2))
### }
#' 
#' 
#' @export

silentMtx = function(mtx, ctrl = 0, dig = 6, 
    wt = c(1.2, 1.1, 1.05, 1), sumwt = 4) {
  len = length(ctrl)
  n = NROW(mtx)
  p = NCOL(mtx)
  if (p < 2) 
    stop("too few columns in mtx input to silentPairs")
  out=matrix(100,nrow=p,ncol=p)
  
  for (i in 1:(p-1)){
  ni=seq(i:p)
  mtxx=mtx[,ni]
  ss=silentPairs(mtxx, ctrl = ctrl, dig = dig, 
        wt = wt, sumwt = sumwt) 
  ssx=round(ss*100/3.175,3)
  ii=seq( (i+1),p)
  lenii=length(ii)
  if (lenii>=1){
#    print(c("i=",i,ii))
#    print(c("ss=",ss))
    out[i,ii]=ssx
  out[ii,i]=-ssx
  }}
  print("Negative index means the column named variable kernel-causes row named")
  print("Positive index means the row named variable kernel-causes column named")
  print ("abs(index)=sign unanimity by weighted sum of 3 signs from Cr1 to Cr3")
  colnames(out)=colnames(mtx)
  rownames(out)=colnames(mtx)
  return(out)
}



