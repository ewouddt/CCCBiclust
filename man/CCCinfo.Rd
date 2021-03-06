% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{CCCinfo}
\alias{CCCinfo}
\title{Filter CCCinfo in info slot of CCC Result}
\usage{
CCCinfo(resbic, method = "none", alpha = 0.01, filter_overlap = FALSE)
}
\arguments{
\item{resbic}{Result from CCC.}

\item{method}{Adjust p-values for multiplicity. Can be one of the following: \code{c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")}.}

\item{alpha}{Significance level.}

\item{filter_overlap}{Filter Bicluster results based on overlap set in CCC algorithm.}
}
\value{
Returns a filtered and ordered CCCinfo (info slot of CCC result)
}
\description{
Compute adjusted p-values and order biclusters based on them.
}
\examples{
\dontrun{
data(heatdata)
out <- CCC(heatdata)
CCCinfo(out,method="bonferroni")
}
}
\author{
Ewoud De Troyer
}

