#' Compute confidence intervals [quantile(s)] of indexes from bootPairs output
#' 
#' Begin with the output of bootPairs function, a (n999 by p-1) matrix when
#' there are p columns of data, \code{bootQuantile} produces a (k by p-1) mtx
#' of quantile(s) of bootstrap ouput assuming that there are k quantiles needed.
#' 
#' @param out {output from bootPairs with p-1 columns and n999 rows}
#' @param probs {quantile evaluation probabilities. The default is k=2,
#' probs=c(.025,0.975) for a 95 percent confidence interval. Note
#' that there are k=2 quantiles desired for each column with this specification}
#' @param per100 {logical (default per100=TRUE) to change the range of
#' 'sum' to [-100, 100] values which are easier to interpret}
#' @return  CI {k quantiles evaluated at probs as a matrix with k rows 
#' and quantile of pairwise p-1 indexes representing p-1 column pairs
#' (fixing the first column in each pair)
#' This function summarizes the 
#' output of of \code{bootPairs(mtx)} (a n999 by p-1 matrix)
#' each containing resampled `sum' values summarizing the weighted sums 
#' associated with all three  criteria from the 
#' function \code{silentPairs(mtx)}
#' applied to each bootstrap sample separately.} #' 
#' 
#' @importFrom stats  quantile
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{silentPairs}}.
#' @references Vinod, H. D. `Generalized Correlation and Kernel Causality with 
#'  Applications in Development Economics' in Communications in 
#'  Statistics -Simulation and Computation, 2015, 
#'  \doi{10.1080/03610918.2015.1122048} 
#' @references Vinod, H. D. and Lopez-de-Lacalle, J. (2009). 'Maximum entropy bootstrap
#'  for time series: The meboot R package.' Journal of Statistical Software,
#'  Vol. 29(5), pp. 1-19. 
#' @references Vinod, H. D. Causal Paths and Exogeneity Tests 
#' in {Generalcorr} Package for Air Pollution and Monetary Policy 
#' (June 6, 2017). Available at SSRN: \url{https://www.ssrn.com/abstract=2982128}
#' @concept bootstrap confidence intervals
#' @concept meboot
#' @concept  kernel regression
#' @concept pairwise comparisons
#' @examples
#' \dontrun{
#' options(np.messages = FALSE)
#' set.seed(34);x=sample(1:10);y=sample(2:11)
#' bb=bootPairs(cbind(x,y),n999=29)
#' bootQuantile(bb) #gives summary stats for n999 bootstrap sum computations
#' 
#' bb=bootPairs(airquality,n999=999);options(np.messages=FALSE)
#' bootQuantile(bb,tau=0.476)#signs for n999 bootstrap sum computations
#' 
#' data('EuroCrime')
#' attach(EuroCrime)
#' bb=bootPairs(cbind(crim,off),n999=29) #col.1= crim causes off 
#' #hence positive signs are more intuitively meaningful.
#' #note that n999=29 is too small for real problems, chosen for quickness here.
#' bootQuantile(bb)# quantile matrix for n999 bootstrap sum computations
#' }
#' @export

bootQuantile=function(out,probs=c(0.025, 0.975),per100=TRUE) {
  out2=as.matrix(out)
  if(per100) {out2=as.matrix(out)*(100/3.175)}  
  pm1=NCOL(out2)  #p m=minus 1
  k=length(probs)
  CI=matrix(NA,nrow=k,ncol=pm1) #place to keep quantiles
  if (pm1==1){
    zj=as.numeric(out2[,1])
    qu=quantile(zj,probs = probs,na.rm = TRUE)
    CI[,1]=qu    
  }
  if (pm1>1) {
    for (j in 1:(pm1)) {
      zj=as.numeric(out2[,j]) #one column at a time
      qu=quantile(zj,probs = probs,na.rm = TRUE)
      CI[,j]=qu
    }  #end j loop
  }# end pm1 loop
  colnames(CI)=colnames(out)
  rownames(CI)=names(qu)
  return(CI)
}


