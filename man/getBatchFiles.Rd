% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitPolyTools.R
\name{getBatchFiles}
\alias{getBatchFiles}
\title{get the names of all batch files present}
\usage{
getBatchFiles(filePrefix, set)
}
\arguments{
\item{filePrefix}{the same filePrefix as used in the fitMarkers call}

\item{set}{one of "models" or "scores"}
}
\value{
a vector of batch file names, either the scores set or the models set
}
\description{
This function gets the names of all batch files already produced
by fitMarkers; used to recover data after a crash.
}
