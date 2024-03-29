% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ability.R
\name{ability.mml}
\alias{ability.mml}
\title{Abilities for a 1 and 2pl}
\usage{
ability.mml(dataSrc, parms, predicate = NULL, method = c("MLE", "WLE"))
}
\arguments{
\item{dataSrc}{a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score}

\item{parms}{object produced by function fit_1pl or fit_2pl or possibly a data.frame of parameters}

\item{predicate}{An optional expression to subset data, if NULL all data is used}

\item{method}{Maximum Likelihood (MLE), Weighted Likelihood (WLE)}
}
\value{
data.frame with variables person_id, theta and se
}
\description{
Note that in a 2PL it is possible for a person with fewer items correct to get a higher ability estimate than
a person with more items correct, on the same test. Use of a 2PL for scoring a summative test should therefore be deemed unethical.
}
\details{
When using a data.frame of parameters, be sure that you use the correct parametrisation. See \code{\link{fit_1pl}} for details.
}
