% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parcorHijk2.R
\name{parcorHijk2}
\alias{parcorHijk2}
\title{Generalized partial correlation coefficients between 
Xi and Xj,}
\usage{
parcorHijk2(xi, xj, xk)
}
\arguments{
\item{xi}{{Input vector of data for variable xi}}

\item{xj}{{Input vector of data for variable xj}}

\item{xk}{{Input data for variables in xk, usually control variables}}
}
\value{
\item{ouij}{Generalized partial correlation Xi with Xj (=cause) after removing xk}
\item{ouji}{Generalized partial correlation Xj with Xi (=cause) after removing xk}
allowing for control variables.
}
\description{
The 2 in the name of the function means second version.
The H in the function name means hybrid. This removes the
effect of Xk, via OLS regression residuals.
This function uses data on two column vectors, xi, xj, and a third set
xk, which can be a vector or a matrix, usually of the remaining 
variables in the model, including control variables, if any.
It first removes missing data from all input variables. Then,
it computes residuals of OLS regression (xi on xk) and (xj on xk).
The function reports the generalized correlation between two OLS residuals.
This hybrid version uses both OLS and then generalized correlation among
OLS residuals. This second version works when `parcorVecH` fails.
It is called by the function `parcorVecH2'.
}
\note{
This function calls \code{\link{kern}},
}
\examples{

\dontrun{
set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
options(np.messages=FALSE)
parcorHijk2(x[,1], x[,2], x[,3])
}#' 
}
\seealso{
See  \code{\link{parcor_ijk}}.
}
\author{
Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
}
