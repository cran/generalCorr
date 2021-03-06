% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abs_stdrhserr.R
\name{abs_stdrhserr}
\alias{abs_stdrhserr}
\title{Absolute values of Hausman-Wu null in kernel regressions of x on y when
both x and y are standardized.}
\usage{
abs_stdrhserr(x, y)
}
\arguments{
\item{x}{{vector of data on the dependent variable}}

\item{y}{{data on the regressors which can be a matrix}}
}
\value{
Absolute values of kernel regression RHS*residuals are returned after
standardizing the data on both sides so that the magnitudes of 
Hausman-Wu null values are comparable between regression of x on y on
the one hand and flipped regression of y on x on the other.
}
\description{
1) standardize the data to force mean zero and variance unity, 2) kernel
regress x on y, with the option `gradients = TRUE' and finally 3) compute
the absolute values of Hausman-Wu null hypothesis for testing exogeneity,
or E(RHS.regressor*error)=0 where error is approximated by kernel 
regression residuals
}
\details{
The first argument is assumed to be the dependent variable.  If
\code{abs_stdrhserr(x,y)} is used, you are regressing x on y (not the usual y
on x). The regressors can be a matrix with 2 or more columns. The missing values
are suitably ignored by the standardization.
}
\examples{
\dontrun{
set.seed(330)
x=sample(20:50)
y=sample(20:50)
abs_stdrhserr(x,y)
}
}
\author{
Prof. H. D. Vinod, Economics Dept., Fordham University, NY
}
\concept{Hausman-Wu statistic}
\concept{kernel regression}
