#' Function to do pairwise deletion of missing rows. 
#' 
#' The aim in pair-wise deletions is to retain the largest 
#' number of available data pairs with all non-missing data.
#' 
#' @param x Vector of x data
#' @param y Vector of y data
#' @return 
#' \item{newx}{A new vector x after removing pairwise missing data} 
#' \item{newy}{A new vector y after removing pairwise missing data} 
## @note %% ~~further notes~~
#' @importFrom stats complete.cases
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @examples
#' 
#' \dontrun{
#' x=sample(1:10);y=sample(1:10);x[2]=NA; y[3]=NA
#' napair(x,y)}
#' 
#' @export


napair <- function(x, y) {
    # author: H D Vinod, Fordham University, 2013
    ok=complete.cases(x,y)
    list(newx = x[ok], newy = y[ok])  #delete NAs from x and y
} 
