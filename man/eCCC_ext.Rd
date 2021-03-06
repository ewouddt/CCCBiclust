% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{eCCC_ext}
\alias{eCCC_ext}
\title{The Extended e-CCC Algorithm}
\usage{
eCCC_ext(matrix, minr = 1, minc = 1, maxErrors = 1, overlap = 0.1,
  missings = "allow", anticorrelation = FALSE, restrictedErrors = FALSE)
}
\arguments{
\item{matrix}{The input matrix in which the columns are ordered by time.}

\item{minr}{Integer containing the row quorum (= minimum number of genes allowed in e-CCC-Biclusters).}

\item{minc}{Integer containing the column quorum (= minimum number of contiguous time points allowed in e-CCC-Biclusters).}

\item{maxErrors}{Integer containing the amount of errors allowed, per gene, in the e-CCC-Biclustering algorithm (value of e).}

\item{overlap}{Numeric in ]0,1[ containing the maximum percentage of overlapping allowed. The \code{info} slot will contain a logical vector which can filter the biclusters which overlap more than \code{overlap*100} percent.}

\item{missings}{How to handle missing values:
\itemize{
\item \code{"remove"}: remove genes with missing values
\item \code{"allow"}: allow missing values as valid errors
\item \code{"jump"}: 'jump over' missing values
}}

\item{anticorrelation}{Logical value to allow or not allow anticorrelation. If allowed (\code{TRUE}), the algorithm will look for e-CCC-Biclusters with Sign-Changes.}

\item{restrictedErrors}{Logical value restricting or not restricting errors. 
If restricting (\code{TRUE}), errors are restricted to the symbols in the 1-neighbourhood of the symbols in 
the alphabet. Since the alphabet {D,N,U} is used in the predefined discretization step provided in this version 
of the prototype, the number of neighbours used in the restricted errors extension can only be equal to 1.}
}
\value{
A Biclust S4 Class object containing extra information of the e-CCC algorithm result in the \code{info} slot.
}
\description{
A R-wrapper which directly calls the original Java code for the Extended e-CCC algorithm (\url{http://kdbio.inesc-id.pt/software/e-ccc-biclustering/}) and transforms it to the output format of the \code{Biclust} R package.
}
\details{
PLACEHOLDER
}
\examples{
\dontrun{
data(heatdata)
out <- eCCC_ext(heatdata,minr=3,minc=2,missings="jump",anticorrelation = TRUE)
eCCCinfo(out,filter="Bonf0.01")
out@info$eCCCpatterns$BC53[1:10,]
}
}
\author{
Ewoud De Troyer
}
\references{
Sara C. Madeira and Arlindo L. Oliveira, "A polynomial time biclustering algorithm for finding genes with approximate expression patterns in gene expression time series", Algorithms for Molecular Biology 2009, 4:8 (4 June 2009)
}

