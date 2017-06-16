#' Intermediate weighting function giving Non-Expected Utility theory weights.
#'
#' Computes cumulative probabilities and difference between consecutive
#' cumulative probabilities described in Vinod (2008) textbook.  This is a simpler version
#' of the version in the book without mapping to non-expected utility theory weights.
#'
#' @param n {A (usually small) integer.}
#' @return
#' \item{x}{sequence 1:n}
#' \item{p}{probabilities p= x[i]/n}
#' \item{pdif}{consecutive differences p[i] - p[i - 1]}
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @references Vinod, H. D. `Hands-On Intermediate Econometrics
#' Using R'  (2008) World Scientific Publishers: Hackensack, NJ.
#' \url{http://www.worldscibooks.com/economics/6895.html}
#' @keywords Prelec
#' @examples
#'
#'  \dontrun{prelec2(10)}
#'
#' @export

prelec2 <- function(n) 
  {
    x <- seq_len(n)
    p <- x/n
    pdif <- c(p[1], diff(p))
    
    list(x=x, p=p, pdif=pdif) 
  }

 
