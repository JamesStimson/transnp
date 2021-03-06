% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{filterInfectedDate}
\alias{filterInfectedDate}
\title{Filter record on specified date of infection}
\usage{
filterInfectedDate(record, case, infDate, keepAfter = TRUE)
}
\arguments{
\item{record}{MCMC output produced by inferTTree}

\item{case}{Case index, either an integer or a string matching one of the case names in the data}

\item{infDate}{Date for filtering}

\item{keepAfter}{Flag, TRUE to keep cases infected after infDate, FALSE to keep cases infected on or before}
}
\value{
Filtered record
}
\description{
Filter record on specified date of infection
}
