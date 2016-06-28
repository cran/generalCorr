#' Function for kernel causality into 3-column matrix admitting control variables 
#' 
#' Allowing input matrix of control variables, this function produce a 3 column matrix
#' summarizing the results where the estimated signs of
#' stochastic dominance order values (+1, 0, -1) are weighted by \code{wt=c(1.2,1.1, 1.05, 1)} to
#' compute an overall result for all orders of stochastic dominance by a weighted sum for
#' the crietria Cr1 and Cr2 and added to the Cr3 estimate as: (+1, 0, -1).
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
#' @param mtx {The data matrix with many columns, y the first column is fixed and then 
#'  paired with all other columns, one by one and called x.}
#' @param ctrl {data matrix for designated control variable(s) outside causal paths}
#' @param dig {Number of digits for reporting (default \code{dig}=6).}
#' @param wt {Allows user to choose a vector of four alternative weights for SD1 to SD4.}
#' @param sumwt { Sum of weights can be changed here =4(default).}
#' @return If there are p columns in the input matrix, there are p-1 possible
#' flipped X, Y pairs leaving the control variables alone.  This function
#' produces p-1 numbers representing the summary sign `sum' from 
#' the signs sg1 to sg3 associated with the three criteria:
#' Cr1, Cr2 and Cr3.  Note that sg1 and sg2 themselves are weithted signs using
#' weighted sum of signs from four orders of stochastic dominance.
#' In general, a positive sign for weighted sum reported in the column `sum' means
#' that the first variable listed as the input to this function is the `kernel cause.' 
#' This function is a summary of \code{someCPairs} allowing for control variables. 
#' For example, crime `kernel causes' police officer deployment (not vice versa) is indicated by 
#' the positive sign of `sum' (=3.175) reported for that example included in this package.
#' @importFrom xtable xtable
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{bootPairs}}.
#' @seealso See  \code{\link{someCPairs}}, \code{\link{some0Pairs}}
#' @references 'Generalized Correlation and Kernel Causality with
#'    Applications in Development Economics' in Communications in
#'    Statistics -Simulation and Computation, 2015,
#'    \url{http://dx.doi.org/10.1080/03610918.2015.1122048}
#' @keywords causal criteria, SD1, SD2, SD3, SD4, generalized correlations
#' @note The output matrix last column for `mtcars' example
#' has the sum of the scores by the three criteria
#' combined. If `sum' is positive, then variable X (mpg) is more likely to have been
#' engineerd to kernel cause the response variable Y, rather than vice versa.
#' @note The European Crime data has all three criteria correctly suggesting that
#' high crime rate kernel causes the deployment of a large number of police officers.
#' 
#' @examples
#'
#'
#' \dontrun{
#' silentPairs(mtcars[,1:3],ctrl=mtcars[4:5]) # mpg paired with others
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
#'silentPairs(cbind(x2,y2), cbind(z,w2))
### }
#' 
#' 
#' @export

silentPairs = function(mtx, ctrl = 0, dig = 6, wt = c(1.2, 1.1, 1.05, 
    1), sumwt = 4) {
    len = length(ctrl)
    n = NROW(mtx)
    p = NCOL(mtx)
    if (p < 2) 
        stop("too few columns in mtx input to silentPairs")
    npair = p - 1
    cr1 = rep(NA, npair)
    cr2 = rep(NA, npair)
    cr3 = rep(NA, npair)
    crall = rep(NA, npair)
    for (typ in 1:3) {
        # outcause = matrix(NA, nrow = npair, ncol = 7) ii = 0
        for (i in 2:p) {
            x0 = mtx[, i]
            y0 = mtx[, 1]
            if (len > 1) {
                z0 = ctrl
                na2 = naTriplet(x0, y0, z0)
                x = na2$newx
                y = na2$newy
                z = na2$newctrl
            }
            if (len == 1) {
                na2 = napair(x0, y0)
                x = na2$newx
                y = na2$newy
            }
            
            
            # if (verbo) { if (i > 2) print(c('i=', i, 'non-missing y=',
            # length(y)), quote = FALSE) }
            if (length(x) < 5) {
                print("available observations<5")
                break
            }
            # ii = ii + 1
            im1 = i - 1
            if (len > 1) {
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
                  round.crit4 = round(crit4, dig)
                  wtdsign = wt * sign(round.crit4)
                  av.crit4 = round((sum(wtdsign, na.rm = TRUE)/sumwt), 
                    dig)
                  if (typ == 1) {
                    cr1[im1] = av.crit4
                  }
                  if (typ == 2) {
                    cr2[im1] = av.crit4
                  }
                }
                if (typ == 3) {
                  gmc0 = gmcmtx0(cbind(x, y, z))
                  out2 = parcorSilent(gmc0, idep = 2)
                  rxy = as.numeric(out2[1, 3])   #third column has r*xy
                  ryx = as.numeric(out2[1, 4]) #col. 4 has r*yx
            del = rxy^2 - ryx^2
 #            print(c('delta',del),q=FALSE)
            cr3[im1] = as.numeric(sign(del))
                }
            }
            
            if (len == 1) {
                if (typ == 1) 
                  arxy = abs_stdapd(x, y)
                if (typ == 1) 
                  aryx = abs_stdapd(y, x)
                if (typ == 2) 
                  arxy = abs_stdres(x, y)
                if (typ == 2) 
                  aryx = abs_stdres(y, x)
                if (typ < 3) {
                  crit4 = comp_portfo2(arxy, aryx)
                  round.crit4 = round(crit4, dig)
                  wtdsign = wt * sign(round.crit4)
                  av.crit4 = round((sum(wtdsign, na.rm = TRUE)/sumwt), 
                    dig)
                  if (typ == 1) {
                    cr1[im1] = av.crit4
                  }
                  if (typ == 2) {
                    cr2[im1] = av.crit4
                  }
                }
                if (typ == 3) {
                  gmc0 = gmcmtx0(cbind(x, y))
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
 


