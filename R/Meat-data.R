#' Meat Data
#'
#' This is the near-infrared spectroscopic meat data used in Murphy, Dean and Raftery (2009) <doi:10.1214/09-AOAS279> and originally collected by McElhinney, Downey and Fearn (1999) <doi:10.1255/jnirs.245>.
#'
#' @docType data
#'
#' @usage data(Meat)
#'
#' @format A list with two components:
#' \describe{
#' \item{x}{Homogenized raw meat spectra. A matrix with 231 rows and 1050 columns.}
#' \item{y}{A vector containing the true class memberships.}}
#'
#' @keywords datasets
#'
#' @references Murphy, Dean and Raftery (2010) <doi:10.1214/09-AOAS279>
#'
#' @source McElhinney, Downey and Fearn (1999) <doi:10.1255/jnirs.245>
#'
#' @examples
#' data(Meat)
#' Meat$x[1:5,1:5]
#' Meat$y
"Meat"
