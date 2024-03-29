#' Matrix R* of generalized correlation coefficients captures nonlinearities.
#' 
#' This function checks for missing data for each pair individually. It then uses the
#' \code{kern} function to kernel regress x on y, and conversely y on x. It
#'  needs the R package `np', which reports the R-squares of each regression. 
#' \code{gmcmtx0()} function
#' reports their square roots after assigning them the observed sign of the Pearson 
#' correlation coefficient. Its threefold advantages are: (i)
#' It is asymmetric, yielding causal direction information
#' by relaxing the assumption of linearity implicit in usual correlation coefficients.
#' (ii) The r* correlation coefficients are generally larger upon admitting 
#' arbitrary nonlinearities.  (iii) max(|R*ij|, |R*ji|) measures (nonlinear) 
#' dependence.
#' For example, let x=1:20 and y=sin(x). This y has a perfect (100 percent)
#' nonlinear dependence on x, and yet Pearson correlation coefficient r(xy)
#' -0.0948372 is near zero, and the 95\% confidence interval (-0.516, 0.363)
#' includes zero, implying that r(xy) is not significantly different from zero.  
#' This shows a miserable failure of traditional r(x,y) to measure dependence
#' when nonlinearities are present. 
#' \code{gmcmtx0(cbind(x,y))} will correctly reveal
#'  perfect (nonlinear) dependence with generalized correlation coefficient =-1.
#' 
#' @param mym {A matrix of data on variables in columns}
#' @param nam {Column names of the variables in the data matrix}
#' @importFrom stats cor
#' @return A non-symmetric R* matrix of generalized correlation coefficients
### @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also as \code{\link{gmcmtxBlk}} for a more general version using
#' blocking allowing several bandwidths.
#' @references Vinod, H. D.'Generalized Correlation and Kernel Causality with 
#'  Applications in Development Economics' in Communications in 
#'  Statistics -Simulation and Computation, 2015, 
#'  \doi{10.1080/03610918.2015.1122048} 
#' @references Vinod, H. D. 'Matrix Algebra Topics in Statistics and Economics
#' Using R', Chapter 4 in 'Handbook of Statistics: Computational Statistics
#' with R', Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#' @references Vinod, H. D. 'New exogeneity tests and causal paths,'
#'  Chapter 2 in 'Handbook of Statistics: Conceptual Econometrics 
#' Using R', Vol.32, co-editors: H. D. Vinod and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2019, pp. 33-64.
#' @references Zheng, S., Shi, N.-Z., and Zhang, Z. (2012). 'Generalized measures 
#'  of correlation for asymmetry, nonlinearity, and beyond,' 
#'  Journal of the American Statistical Association, vol. 107, pp. 1239-1252.
#' @concept  kernel regression
#' @concept R* asymmetric matrix of generalized correlation coefficients
#' @examples
#' 
#' gmcmtx0(mtcars[,1:3])
#' 
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c('V1', 'v2', 'V3')
#' gmcmtx0(x)}
#' 
#' @export

gmcmtx0 <- function(mym, nam = colnames(mym)) {
    # mym is a data matrix with n rows and p columns some NAs may be present in the
    # matrix
    p = NCOL(mym)
    # print(c('p=',p))
    out1 = matrix(1, p, p)  # out1 stores asymmetric correlations
    for (i in 1:p) {
        x = mym[, i]
        for (j in 1:p) {
            if (j > i) 
                {
                  y = mym[, j]
                  ava.x = which(!is.na(x))  #ava means available
                  ava.y = which(!is.na(y))
                  ava.both = intersect(ava.x, ava.y)
                  newx = x[ava.both]
                  newy = y[ava.both]
                  c1 = cor(newx, newy)
                  sig = sign(c1)
                  # begin non parametric regressions
                  bw = npregbw(formula = newx ~ newy, tol = 0.1, ftol = 0.1)
                  mod.1 = npreg(bws = bw, gradients = FALSE, residuals = TRUE)
                  corxy = sqrt(mod.1$R2) * sig
                  out1[i, j] = corxy
                  bw2 = npregbw(formula = newy ~ newx, tol = 0.1, ftol = 0.1)
                  mod.2 = npreg(bws = bw2, gradients = FALSE, residuals = TRUE)
                  coryx = sqrt(mod.2$R2) * sig
                  out1[j, i] = coryx
                }  #end i loop
        }  #end j loop
    }  #endif
    colnames(out1) = nam
    rownames(out1) = nam
    return(out1)
}
