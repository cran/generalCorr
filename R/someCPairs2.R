#' Kernel causality computations admitting control variables reporting 
#' a 7-column matrix, version 2.
#' 
#' Second version of \code{someCPairs} also allows input matrix of 
#' control variables, produce 7 column matrix
#' summarizing the results where the signs of
#' stochastic dominance order values (+1 or -1) are weighted by 
#' \code{wt=c(1.2,1.1, 1.05, 1)} to
#' compute an overall result for all orders of stochastic dominance by 
#' a weighted sum for the criteria Cr1 and Cr2. 
#' The weighting is obviously not needed for the third criterion Cr3.
#' 
#' The reason for slightly declining weights on the signs from
#' SD1 to SD4 is simply that the local mean comparisons 
#' implicit in SD1 are known to be
#' more reliable than local variance implicit in SD2, local skewness implicit in
#' SD3 and local kurtosis implicit in SD4. The source of slightly declining sampling
#' unreliability of higher moments is the
#' higher power of the deviations from the mean needed in their computations.
#' The summary results for all
#' three criteria are reported in one matrix called \code{outVote}: 
#'   
#' (typ=1) reports ('Y', 'X', 'Cause',
#' 'SD1.rhserr', 'SD2.rhserr', 'SD3.rhserr', 'SD4.rhserr') 
#' naming variables identifying the 'cause'
#' and measures of stochastic dominance using absolute values of kernel
#' regression abs(RHS first regressor*residual) values
#' comparing flipped regressions X on Y versus Y on X.
#' The letter C in the titles reminds presence of control variable(s).
#' 
#' 
#' typ=2 reports ('Y', 'X', 'Cause', 'SD1resC', 'SD2resC', 'SD3resC', 'SD4resC')
#' and measures of stochastic dominance using absolute values of kernel
#' regression residuals comparing regression of X on Y with that of Y on X.
#' 
#' 
#' typ=3 reports ('Y', 'X', 'Cause', 'r*x|yC', 'r*y|xC', 'r', 'p-val')
#' containing generalized correlation coefficients r*, 'r' refers to.
#' Pearson correlation coefficient p-val is the p-value for 
#' testing the significance of 'r'. 
#' The letter C in the titles reminds the presence of control variable(s).
#' @param mtx {The data matrix with many columns where the first column is fixed and then 
#'  paired with all other columns, one by one.}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths}
#' @param dig {Number of digits for reporting (default \code{dig}=6).}
#' @param verbo {Make \code{verbo= TRUE} for printing detailed steps.}
#' @param rnam {Make \code{rnam= TRUE} if cleverly created rownames are desired.}
#' @param wt {Allows user to choose a vector of four alternative weights for SD1 to SD4.}
#' @param sumwt { Sum of weights can be changed here =4(default).}
#' @return Prints three matrices detailing results for Cr1, Cr2 and Cr3.
#' It also returns a grand summary matrix called `outVote' which summarizes all three criteria.
#' In general, a positive sign for weighted sum reported in the column `sum' means
#' that the first variable listed as the input to this function is the `kernel cause.' 
#' This function is an extension of \code{some0Pairs} to allow for control variables. 
#' For example, crime `kernel causes' police officer deployment (not vice versa) is indicated by 
#' the positive sign of `sum' (=3.175) reported for that example included in this package.
#' @importFrom xtable xtable
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See Also \code{\link{somePairs}}, \code{\link{some0Pairs}}
#' @references Vinod, H. D. 'Generalized Correlation and Kernel Causality with
#'    Applications in Development Economics' in Communications in
#'    Statistics -Simulation and Computation, 2015,
#'    \url{http://dx.doi.org/10.1080/03610918.2015.1122048}
#' @concept causal criteria
#' @concept stochastic dominance
#' @concept  generalized correlations
#' @note The output matrix last column for `mtcars' example
#' has the sum of the scores by the three criteria
#' combined. If `sum' is positive, then variable X (mpg) is more likely to have been
#' engineered to kernel cause the response variable Y, rather than vice versa.
#' @note The European Crime data has all three criteria correctly suggesting that
#' high crime rate kernel causes the deployment of a large number of police officers.
#' 
#' @examples
#'
#'
#' \dontrun{
#' someCPairs2(mtcars[,1:3],ctrl=mtcars[4:5]) # first variable is mpg and effect on mpg is of interest
#' }
#' 
### \dontrun{
#' set.seed(234)
#' z=runif(10,2,11)# z is independently created
#' x=sample(1:10)+z/10  #x is somewhat indep and affected by z
#' y=1+2*x+3*z+rnorm(10)
#' w=runif(10)
#' x2=x;x2[4]=NA;y2=y;y2[8]=NA;w2=w;w2[4]=NA
#' someCPairs2(cbind(x2,y2), cbind(z,w2)) #yields x2 as correct cause
### }
#' 
#' 
#' @export

someCPairs2=
  function (mtx, ctrl, dig = 6, verbo = TRUE, rnam = FALSE, wt = c(1.2, 
                                                                   1.1, 1.05, 1), sumwt = 4) 
  {
    n = NROW(mtx)
    p = NCOL(mtx)
    if (p < 2) stop("too few columns in mtx input to someCPairs")
    k = NCOL(ctrl) #set of column numbers representing control variables
    npair = p - 1 
    cr1 = rep(NA, npair)
    cr2 = rep(NA, npair)
    cr3 = rep(NA, npair)
    crall = rep(NA, npair)
    rna = rep(NA, p)
    outVote = matrix(NA, nrow = npair, ncol = 7)
    colnames(outVote) = c("X", "Y", "Cause", "Cr1", "Cr2", "Cr3", 
                          "sum")
    nam = colnames(mtx)
    if (length(nam) == 0) 
      nam = paste("V", 1:p, sep = "")
    for (typ in 1:3) {
      outcause = matrix(NA, nrow = npair, ncol = 7)
      ii = 0
      for (i in 2:p) {
        x0 = mtx[, i]
        y0 = mtx[, 1]
        z0 = ctrl
        na2 = naTriplet(x0, y0, z0)
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
        rna[i] = paste("1", i, sep = ".")
        if (typ == 1) 
          outVote[ii, 1] = nam[1]
        if (typ == 1) 
          outVote[ii, 2] = nam[i]
        if (typ == 1) 
          arxy = abs_stdapdC(x, y, z)
        if (typ == 1) 
          aryx = abs_stdapdC(y, x, z)
        if (typ == 2) 
          arxy = abs_stdresC(x, y, z)
        if (typ == 2) 
          aryx = abs_stdresC(y, x, z)
        if (typ < 3) {
          crit4 = comp_portfo2(arxy, aryx)
          if (verbo) 
            print(crit4)
          round.crit4 = round(crit4, dig)
          wtdsign = wt * sign(round.crit4)
          av.crit4 = round((sum(wtdsign, na.rm = TRUE)/sumwt), 
                           dig)
          outcause[ii, 4:7] = round.crit4
          outcause[ii, 1] = nam[1]
          outcause[ii, 2] = nam[i]
          outcause[ii, 3] = nam[i]
          if (av.crit4 > 0) 
            outcause[ii, 3] = nam[1]
          if (typ == 1) {
            outVote[ii, 4] = av.crit4
            cr1[ii] = av.crit4
          }
          if (typ == 2) {
            outVote[ii, 5] = av.crit4
            cr2[ii] = av.crit4
          }
        }
        if (typ == 3) {
          cc=cor.test(x, y)
          par1 = parcor_ijk(x, y, z)
          rxy=par1$ouij
          ryx=par1$ouji
          del = rxy^2 - ryx^2
          cr3[ii] = as.numeric(sign(del))
          outcause[ii, 4] = round(rxy, dig)
          outcause[ii, 5] = round(ryx, dig)
          outcause[ii, 6] = round(cc$estimate, dig)
          outcause[ii, 7] = round(cc$p.value, dig)
          outcause[ii, 1] = nam[1]
          outcause[ii, 2] = nam[i]
          outcause[ii, 3] = nam[i]
          if (del > 0) 
            outcause[ii, 3] = nam[1]
          outVote[ii, 6] = sign(del)
        }
      }
      namout = c("Y", "X", "Cause", "SD1", "SD2", "SD3", "SD4")
      colnames(outcause) = namout
      if (typ == 1) 
        namout = c("Y", "X", "Cause", "SD1apdC", "SD2apdC", 
                   "SD3apdC", "SD4apdC")
      if (typ == 2) 
        namout = c("Y", "X", "Cause", "SD1resC", "SD2resC", 
                   "SD3resC", "SD4resC")
      if (typ == 3) 
        namout = c("Y", "X", "Cause", "r*x|yC", "r*y|xC", "r", 
                   "p-val")
      colnames(outcause) = namout
      if (rnam) 
        rownames(outcause) = rna[2:p]
      if (verbo) 
        print(outcause)
      if (verbo) 
        print(xtable(outcause))
      if (typ == 1) 
        outCr1 = outcause
      if (typ == 2) 
        outCr2 = outcause
      if (typ == 3) 
        outCr3 = outcause
    }
    for (j in 1:npair) {
      cr13 = c(cr1[j], cr2[j], cr3[j])
      crall[j] = round(sum(cr13, na.rm = TRUE), dig)
      if (verbo) 
        print(c("Scores", cr1[j], cr2[j], cr3[j], crall[j]), 
              quote = FALSE)
      outVote[j, 7] = crall[j]
      if (!is.na(crall[j])) {
        if (crall[j] > 0) 
          outVote[j, 3] = nam[1]
        if (crall[j] <= 0) 
          outVote[j, 3] = nam[j + 1]
      }
    }
    if (verbo) 
      print(outVote)
    if (verbo) 
      print(xtable(outVote))
    list(outCr1 = outCr1, outCr2 = outCr2, outCr3 = outCr3, outVote = outVote)
  }
