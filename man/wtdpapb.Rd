% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wtdpapb.R
\name{wtdpapb}
\alias{wtdpapb}
\title{Creates input for the stochastic dominance function stochdom2}
\usage{
wtdpapb(xa, xb)
}
\arguments{
\item{xa}{{Vector of (excess) returns for the first investment option A or
values of any random variable being compared to another.}}

\item{xb}{Vector of returns for the second option B}
}
\value{
\item{wpa}{Weighted vector of probabilities for option A} 
\item{wpb}{Weighted vector of probabilities for option B} 
\item{dj}{Vector of interval widths (distances) when both sets of data are forced on a common support}
}
\description{
Stochastic dominance is a sophisticated comparison of two distributions of
stock market returns.  The dominating distribution is superior in terms of
mean, variance, skewness and kurtosis respectively, representing dominance
orders 1 to 4, without directly computing four moments.  Vinod(2008) sec. 4.3
explains the details.  The `wtdpapb' function creates the input
for stochdom2 which in turn computes the stochastic dominance.
See Vinod (2004) for details about quantitative stochastic dominance.
}
\note{
Function is needed before using stochastic dominance

In Vinod (2008) where the purpose of \code{wtdpapb} is to map from standard
`expected utility theory' weights to more sophisticated `non-expected utility
theory' weights using Prelec's (1998, Econometrica, p. 497) method.  These
 weights are not needed here. Hence we provide the function \code{prelec2}
 which does not use Prelec weights at all, thereby simplifying and speeding up
 the R code provided in Vinod (2008). This function avoids sophisticated `non-expected'
 utility theory which incorporates commonly observed human behavior favoring
 loss aversion and other anomalies inconsistent with precepts of the
 expected utility theory. Such weighting is not needed for our application.
}
\examples{

 \dontrun{
 set.seed(234);x=sample(1:30);y=sample(5:34)
 wtdpapb(x,y)}

}
\references{
Vinod, H. D.', 'Hands-On Intermediate Econometrics 
Using R'  (2008) World Scientific Publishers: Hackensack, NJ.
\url{https://www.worldscientific.com/worldscibooks/10.1142/12831}

Vinod, H. D. 'Ranking Mutual Funds Using 
Unconventional Utility Theory and Stochastic Dominance,'
Journal of Empirical Finance Vol. 11(3) 2004, pp. 353-377.
}
\seealso{
See Also \code{\link{stochdom2}}
}
\author{
Prof. H. D. Vinod, Economics Dept., Fordham University, NY
}
\concept{stochastic dominance}
