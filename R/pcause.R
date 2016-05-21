#' Compute the bootstrap probability of correct causal direction.
#' 
#' Maximum entropy bootstrap (meboot) package is used for statistical inference
#' regarding \eqn{\delta} which equals GMC(X|Y)-GMC(Y|X) defined by Zheng et al (2012).
#' The bootstrap provides an approximation to chances of correct determination of
#' the causal direction.
#' 
#' @param x {Vector of x data}
#' @param y {Vector of y data}
#' @param n999 {Number of bootstrap replications (default=999)}
#' @importFrom meboot meboot
#' @return  P(cause) the bootstrap proportion of correct causal determinations.
#' @note 'pcause' is computer intensive and generally slow. It may be better to use
#'   it at a later stage in the investigation when a preliminary causal determination 
#'   is already made.  Its use may slow the exploratory phase. In my experience, if
#'   P(cause) is less than 0.55, there is a cause for concern.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @references Vinod, H. D.'Generalized Correlation and Kernel Causality with 
#'  Applications in Development Economics' in Communications in 
#'  Statistics -Simulation and Computation, 2015, 
#'  \url{http://dx.doi.org/10.1080/03610918.2015.1122048} 
#' @references Zheng, S., Shi, N.-Z., and Zhang, Z. (2012). Generalized measures 
#'  of correlation for asymmetry, nonlinearity, and beyond. 
#'  Journal of the American Statistical Association, vol. 107, pp. 1239-1252.
#' @references Vinod, H. D. and Lopez-de-Lacalle, J. (2009). 'Maximum entropy bootstrap
#'  for time series: The meboot R package.' Journal of Statistical Software,
#'  Vol. 29(5), pp. 1-19. 
#' @keywords kernel regression, asymmetric
#' @examples
#'
#' \dontrun{
#' set.seed(34);x=sample(1:10);y=sample(2:11)
#' pcause(x,y,n999=29)
#' 
#' data('EuroCrime')
#' attach(EuroCrime)
#' pcause(crim,off,n999=29)
#' }
#' @export

pcause = function(x, y, n999 = 999) {
    if (n999 <= 1) {
        p.cause = NA
        return(p.cause)
    } else {
        xboot = meboot(x = x, reps = n999)$ensemble
        xb = x
        yboot = meboot(x = y, reps = n999)$ensemble
        yb = y
        out.diff = rep(NA, n999)
        out.corxy = rep(NA, n999)
        out.coryx = rep(NA, n999)
        for (i in 1:n999) {
            xb = xboot[, i]
            yb = yboot[, i]
            gm = gmcxy_np(xb, yb)
            out.corxy[i] = gm$corxy
            out.coryx[i] = gm$coryx
            out.diff[i] = out.corxy[i] - out.coryx[i]
        }  #end of i loop
        ou.nega = length(out.diff[out.diff < 0])
        ou.posi = length(out.diff[out.diff > 0])
        p.cause = max(ou.nega, ou.posi)/n999
    }  #end of else
    return(p.cause)
} 
