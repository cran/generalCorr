#' Report causal identification for all pairs of variables in a matrix 
#' (deprecated function). It is better to choose a target variable and pair
#' it with all others, instead of considering all possible targets.
#'
#' This studies all possible (perhaps too many) causal directions in a matrix.
#' It is deprecated because it uses older criterion 1 by caling \code{abs_stdapd}
#' I recommend using \code{causeSummary} or its block version \code{cuseSummBlk}.
#' This uses \code{abs_stdres}, \code{comp_portfo2}, etc. and returns
#' a matrix with 7 columns having detailed output. Criterion 1 has been revised
#' as described in Vinod (2019) and is known to work better.
#'
#' @param mtx {Input matrix with variable names}
#' @param dig {Digits of accuracy in reporting (=6, default)}
#' @param verbo {Logical variable, set to 'TRUE' if printing is desired}
#' @param typ {Causal direction criterion number (typ=1 is default)
#'  Criterion 1 (Cr1) compares kernel regression absolute values of gradients.
#'  Criterion 2 (Cr2) compares kernel regression absolute values of residuals.
#'  Criterion 3 (Cr3) compares kernel regression based r*(x|y) with r*(y|x)}.
#' @param rnam {Logical variable, default \code{rnam=FALSE} means the user does
#' not want the row names to be
#' (somewhat too cleverly) assigned by the function.}
#' @return A 7-column matrix called 'outcause' with names of variables
#'  X and Y in the first two columns and the name of the 'causal' variable in 3rd col.
#'  Remaining four columns report numerical computations of SD1 to SD4, r*(x|y),
#'  r*(y|x).  Pearson r and p-values for its traditional significance testing.
#' @note The cause reported in the third column
#'  is identified from the sign of the first SD1 only,
#'  ignoring SD2, SD3 and SD4  under both Cr1 and Cr2. It is
#'  a good idea to loop a call to this function with typ=1:3. One can print
#'  the resulting 'outcause' matrix with the 
#'  \code{xtable(outcause)} for the Latex output.
#'  A similar deprecated function included in this package,
#'   called \code{some0Pairs}, incorporates all SD1 to SD4 and all
#'  three criteria Cr1 rto Cr3 to report a `sum' of indexes representing the signed 
#'  number whose sign can more comprehensively help determine the causal direction(s).
#'  Since the Cr1 here is revised in later work, this is deprecated.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also  \code{somePairs}, \code{some0Pairs} \code{causeSummary}
#' @references Vinod, H. D.'Generalized Correlation and Kernel Causality with
#'  Applications in Development Economics' in Communications in
#'  Statistics -Simulation and Computation, 2015,
#'  \doi{10.1080/03610918.2015.1122048}
#' @references Vinod, H. D. 'New exogeneity tests and causal paths,'
#'  Chapter 2 in 'Handbook of Statistics: Conceptual Econometrics 
#' Using R', Vol.32, co-editors: H. D. Vinod and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2019, pp. 33-64.
#'  
#' @concept  amorphous partial derivative apd
#' @concept stochastic dominance
#' @concept  absolute residual values
#' @examples
#' 
#' data(mtcars)
#' options(np.messages=FALSE)
#' for(j in 1:3){
#' a1=allPairs(mtcars[,1:3], typ=j)
#' print(a1)}
#'
#' @export


allPairs <- function(mtx, dig = 6, verbo = FALSE, typ = 1, rnam = FALSE) {
    # author: H.D.Vinod, Prof of Economics, Fordham Univ. NY Sept2013 input matrix
    # of data with p columns with colnames input dig=digits for rounding typ==3
    # added for r*
    
    n = NROW(mtx)
    p = NCOL(mtx)
    
    print(c("n,p,digits", n, p, dig))
    if (typ == 1) 
        print("absolute apds compared")
    if (typ == 2) 
        print("absolute residuals compared")
    if (typ == 3) 
        print("r* compared")
    nam = colnames(mtx)
    npair = choose(p, 2)
    rna = rep(0, npair)  #for storing row names
    print(c("no. of pairs, typ ", npair, typ))
    outcause = matrix(NA, nrow = npair, ncol = 7)
    ii = 0
    for (i in 1:p) {
        x0 = mtx[, i]
        for (j in i:p) {
            if (j > i) 
                {
                  y0 = mtx[, j]
                  na2 = napair(x0, y0)
                  x = na2$newx
                  y = na2$newy
                  print(c("i,non-missing", i, length(x)))
                  ii = ii + 1
                  print(c("i,j,ii", i, j, ii))
                  rna[ii] = paste(i, j, sep = ".")
                  if (typ == 1) 
                    arxy = abs_stdapd(x, y)
                  if (typ == 1) 
                    aryx = abs_stdapd(y, x)
                  if (typ == 2) 
                    arxy = abs_stdres(x, y)
                  if (typ == 2) 
                    aryx = abs_stdres(y, x)
                  if (typ < 3) 
                    {
                      crit4 = comp_portfo2(arxy, aryx)
                      if (verbo) 
                        print(crit4)
                      round.crit4 = round(crit4, dig)
                      
                      outcause[ii, 4:7] = round.crit4
                      outcause[ii, 1] = nam[i]
                      outcause[ii, 2] = nam[j]
                      outcause[ii, 3] = nam[i]  #i has x and x is the implicit cause
                      if (crit4[1] > 0) 
                        outcause[ii, 3] = nam[j]  #SD1>0 then cause=y
                    }  #endif typ<3
                  if (typ == 3) 
                    {
                      rst = rstar(x, y)
                      rxy = rst$corxy
                      ryx = rst$coryx
                      # print(c(rxy,ryx,ii,i))
                      del = rxy^2 - ryx^2
                      # del>0 means rxy>ryx or x on y good or cause=y
                      outcause[ii, 4] = round(rst$corxy, dig)
                      outcause[ii, 5] = round(rst$coryx, dig)
                      outcause[ii, 6] = round(rst$pearson.r, dig)
                      outcause[ii, 7] = round(rst$pv, dig)
                      outcause[ii, 1] = nam[i]
                      outcause[ii, 2] = nam[j]
                      outcause[ii, 3] = nam[i]  #cause= x
                      if (del > 0) 
                        outcause[ii, 3] = nam[j]  #cause=y
                    }  #endif typ==3
                }  #end of if  j>i
        }  #end of j loop
    }  #end of i loop
    namout = c("X", "Y", "Cause", "SD1", "SD2", "SD3", "SD4")
    colnames(outcause) = namout
    if (typ == 1) 
        namout = c("X", "Y", "Cause", "SD1apd", "SD2apd", "SD3apd", "SD4apd")
    if (typ == 2) 
        namout = c("X", "Y", "Cause", "SD1res", "SD2res", "SD3res", "SD4res")
    if (typ == 3) 
        namout = c("X", "Y", "Cause", "r*x|y", "r*y|x", "r", "p-val")
    # rownames(outcause)=rownames(out2)
    colnames(outcause) = namout
    if (rnam) 
        rownames(outcause) = rna
    if (verbo) 
        print(outcause)
    if (verbo) 
        print(xtable(outcause))
    return(outcause)
}
