% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parcorBijk.R
\name{parcorBijk}
\alias{parcorBijk}
\title{Block version of generalized partial correlation coefficients between Xi 
and Xj, after removing the
effect of xk, via nonparametric regression residuals.}
\usage{
parcorBijk(xi, xj, xk, blksiz = 10)
}
\arguments{
\item{xi}{{Input vector of data for variable xi}}

\item{xj}{{Input vector of data for variable xj}}

\item{xk}{{Input data for variables in xk, usually control variables}}

\item{blksiz}{{block size, default=10, if chosen blksiz >n, where n=rows in matrix
then blksiz=n. That is, no blocking is done}}
}
\value{
\item{ouij}{Generalized partial correlation Xi with Xj (=cause) after removing xk}
\item{ouji}{Generalized partial correlation Xj with Xi (=cause) after removing xk}
allowing for control variables.
}
\description{
This function uses data on two column vectors, xi, xj and a third
xk which can be a vector or a matrix, usually of the remaining 
variables in the model, including control variables, if any.
It first removes missing data from all input variables. Then,
it computes residuals of kernel regression (xi on xk) and (xj on xk). 
This is a block version of parcor_ijk.
}
\note{
This function calls \code{\link{kern}},
}
\examples{

\dontrun{
set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
options(np.messages=FALSE)
parcorBijk(x[,1], x[,2], x[,3], blksi=10)
}#' 
}
\seealso{
See  \code{\link{parcor_ijk}}.
}
\author{
Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
}
