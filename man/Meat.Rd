% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Meat-data.R
\docType{data}
\name{Meat}
\alias{Meat}
\title{Meat Data}
\format{A list with two components:
\describe{
\item{x}{Homogenized raw meat spectra. A matrix with 231 rows and 1050 columns.}
\item{y}{A vector containing the true class memberships.}}}
\source{
McElhinney, Downey and Fearn (1999) <doi:10.1255/jnirs.245>
}
\usage{
data(Meat)
}
\description{
This is the near-infrared spectroscopic meat data used in Murphy, Dean and Raftery (2009) <doi:10.1214/09-AOAS279> and originally collected by McElhinney, Downey and Fearn (1999) <doi:10.1255/jnirs.245>.
}
\examples{
data(Meat)
Meat$x[1:5,1:5]
Meat$y
}
\references{
Murphy, Dean and Raftery (2010) <doi:10.1214/09-AOAS279>
}
\keyword{datasets}
