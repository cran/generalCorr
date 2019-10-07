#' Compute usual summary stats of 'sum' indexes from bootPairs output
#' 
#' Begin with the output of bootPairs function, a (n999 by p-1) matrix when
#' there are p columns of data, \code{bootSummary} produces a (6 by p-1) mtx
#' of summary of bootstrap ouput (Min, 1st Qu,Median, Mean, 3rd Qi.,Max)
#' 
#' @param out {output from bootPairs with p-1 columns and n999 rows in input here}
#' @param per100 {logical (default per100=TRUE) to change the range of
#' 'sum' to [-100, 100] values which are easier to interpret}
#' @return  summ {summary output from the (n999 by p-1) matrix
#' output of \code{bootPairs(mtx)} 
#' each containing resampled `sum' values summarizing the weighted sums 
#' associated with all three  criteria from the 
#' function \code{silentPairs(mtx)}
#' applied to each bootstrap sample separately.} 
#' 
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{silentPairs}}.
#' @references Vinod, H. D. `Generalized Correlation and Kernel Causality with 
#'  Applications in Development Economics' in Communications in 
#'  Statistics -Simulation and Computation, 2015, 
#'  \url{http://dx.doi.org/10.1080/03610918.2015.1122048} 
#' @references Vinod, H. D. and Lopez-de-Lacalle, J. (2009). 'Maximum entropy bootstrap
#'  for time series: The meboot R package.' Journal of Statistical Software,
#'  Vol. 29(5), pp. 1-19. 
#' @references Vinod, H. D. Causal Paths and Exogeneity Tests 
#' in {Generalcorr} Package for Air Pollution and Monetary Policy 
#' (June 6, 2017). Available at SSRN: \url{https://ssrn.com/abstract=2982128}
#' @concept bootstrap
#' @concept  meboot
#' @concept  kernel regression
#' @concept  pairwise comparisons
#' @examples
#' \dontrun{
#' options(np.messages = FALSE)
#' set.seed(34);x=sample(1:10);y=sample(2:11)
#' bb=bootPairs(cbind(x,y),n999=29)
#' bootSummary(bb) #gives summary stats for n999 bootstrap sum computations
#' 
#' bb=bootPairs(airquality,n999=999);options(np.messages=FALSE)
#' bootSummary(bb)#signs for n999 bootstrap sum computations
#' 
#' data('EuroCrime')
#' attach(EuroCrime)
#' bb=bootPairs(cbind(crim,off),n999=29) #col.1= crim causes off 
#' #hence positive signs are more intuitively meaningful.
#' #note that n999=29 is too small for real problems, chosen for quickness here.
#' bootSummary(bb)#signs for n999 bootstrap sum computations
#' }
#' @export

bootSummary=function(out, per100=TRUE) {
  out2=as.matrix(out)
  if(per100) {out2=as.matrix(out)*(100/3.175)}
  summ=apply(out2,2,summary)
  return(summ)
}

