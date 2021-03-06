% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cofactor.R
\name{cofactor}
\alias{cofactor}
\title{Compute cofactor of a matrix based on row r and column c.}
\usage{
cofactor(x, r, c)
}
\arguments{
\item{x}{{matrix whose cofactor is desired to be computed}}

\item{r}{{row number}}

\item{c}{{column number}}
}
\value{
cofactor of x,  w.r.t. row r and column c.
}
\description{
Compute cofactor of a matrix based on row r and column c.
}
\note{
needs the function `minor'' in memory. attaches sign (-1)^(r+c) to the minor.
}
\examples{

## The function is currently defined as
function (x, r, c) 
{
    out = minor(x, r, c) * ((-1)^(r + c))
    return(out)
  }
}
\seealso{
\code{minor(x,r,c)}
}
\author{
Prof. H. D. Vinod, Economics Dept., Fordham University, NY
}
\concept{cofactor of a matrix}
