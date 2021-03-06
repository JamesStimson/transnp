% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visuals.R
\name{networkTPlot}
\alias{networkTPlot}
\title{Dynamic plot of the transmission network with edges weighted according to likelihood of transmission}
\usage{
networkTPlot(thisRecord, mcmcIndex = 1, burnin = 0.5,
  missLabel = "Unsampled", colours = c("lightblue", "orange"),
  fontSize = 24)
}
\arguments{
\item{thisRecord}{Posterior MCMC sample set from TransPhylo}

\item{mcmcIndex}{Sample set index of network to be displayed}

\item{burnin}{Proportion of the MCMC output for edge weights to be discarded as burnin}

\item{missLabel}{Label to be used for missing/unsampled cases}

\item{colours}{Colours for sampled and unsampled nodes, vector of 2 colours}
}
\description{
Dynamic plot of the transmission network with edges weighted according to likelihood of transmission
}
\examples{
networkTPlot(record[[5]])
}
