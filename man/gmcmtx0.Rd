% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmcmtx0.R
\name{gmcmtx0}
\alias{gmcmtx0}
\title{Matrix R* of generalized correlation coefficients captures nonlinearities.}
\usage{
gmcmtx0(mym, nam = colnames(mym))
}
\arguments{
\item{mym}{{A matrix of data on variables in columns}}

\item{nam}{{Column names of the variables in the data matrix}}
}
\value{
A non-symmetric R* matrix of generalized correlation coefficients
}
\description{
This function checks for missing data for each pair individually. It then uses the
\code{kern} function to kernel regress x on y, and conversely y on x. It
 needs the R package `np', which reports the R-squares of each regression. 
\code{gmcmtx0()} function
reports their square roots after assigning them the observed sign of the Pearson 
correlation coefficient. Its threefold advantages are: (i)
It is asymmetric, yielding causal direction information
by relaxing the assumption of linearity implicit in usual correlation coefficients.
(ii) The r* correlation coefficients are generally larger upon admitting 
arbitrary nonlinearities.  (iii) max(|R*ij|, |R*ji|) measures (nonlinear) 
dependence.
For example, let x=1:20 and y=sin(x). This y has a perfect (100 percent)
nonlinear dependence on x, and yet Pearson correlation coefficient r(xy)
-0.0948372 is near zero, and the 95\% confidence interval (-0.516, 0.363)
includes zero, implying that r(xy) is not significantly different from zero.  
This shows a miserable failure of traditional r(x,y) to measure dependence
when nonlinearities are present. 
\code{gmcmtx0(cbind(x,y))} will correctly reveal
 perfect (nonlinear) dependence with generalized correlation coefficient =-1.
}
\examples{

gmcmtx0(mtcars[,1:3])

\dontrun{
set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
colnames(x)=c('V1', 'v2', 'V3')
gmcmtx0(x)}

}
\references{
Vinod, H. D.'Generalized Correlation and Kernel Causality with 
 Applications in Development Economics' in Communications in 
 Statistics -Simulation and Computation, 2015, 
 \doi{10.1080/03610918.2015.1122048}

Vinod, H. D. 'Matrix Algebra Topics in Statistics and Economics
Using R', Chapter 4 in 'Handbook of Statistics: Computational Statistics
with R', Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
North Holland, Elsevier Science Publishers, 2014, pp. 143-176.

Vinod, H. D. 'New exogeneity tests and causal paths,'
 Chapter 2 in 'Handbook of Statistics: Conceptual Econometrics 
Using R', Vol.32, co-editors: H. D. Vinod and C.R. Rao. New York:
North Holland, Elsevier Science Publishers, 2019, pp. 33-64.

Zheng, S., Shi, N.-Z., and Zhang, Z. (2012). 'Generalized measures 
 of correlation for asymmetry, nonlinearity, and beyond,' 
 Journal of the American Statistical Association, vol. 107, pp. 1239-1252.
}
\seealso{
See Also as \code{\link{gmcmtxBlk}} for a more general version using
blocking allowing several bandwidths.
}
\author{
Prof. H. D. Vinod, Economics Dept., Fordham University, NY
}
\concept{R* asymmetric matrix of generalized correlation coefficients}
\concept{kernel regression}
