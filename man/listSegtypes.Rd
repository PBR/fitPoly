% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitPolyTools.R
\name{listSegtypes}
\alias{listSegtypes}
\title{conversion of segtype code to F1 segregation ratios}
\usage{
listSegtypes(ploidy, ploidy2=NULL)
}
\arguments{
\item{ploidy}{the ploidy of parent 1 of the F1, or of the F1 itself if
ploidy2 is NULL}

\item{ploidy2}{the ploidy of parent 2 of the F1}
}
\value{
a matrix with one row for each segregation type and one column for
each possible F1 dosage, with the integer ratios of the dosages for each
segregation type.
}
\description{
Produce a matrix with the F1 segregation ratios (as integers)
for all segregation types for the given ploidy
}
\details{
If ploidy2 is not NULL, ploidy and ploidy2 are the ploidy levels of
the two parents, and both must be even. If ploidy2 is NULL, ploidy is
the ploidy of the F1; if even, both parents are assumed to have the same
ploidy; if odd, parent 1 and parent 2 are asumed to have ploidy-1 and ploidy+1.
This function calls calcSegtypeInfo and is a convenience function; it is not
used by any other functions. For more information, including parental
dosages for each segregation, use calcSegtypeInfo and segtypeInfoSummary.
}
