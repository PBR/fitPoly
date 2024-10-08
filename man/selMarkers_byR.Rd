% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitPolyTools.R
\name{selMarkers_byR}
\alias{selMarkers_byR}
\title{Select markers at specified R levels}
\usage{
selMarkers_byR(Rstats, Rlevels, mrkperlevel=1, stat="q95")
}
\arguments{
\item{Rstats}{a data frame as returned by calcRstats}

\item{Rlevels}{a vector of R levels in increasing order (preferably from
the minimum to the maximum of stat over all markers in Rstats)}

\item{mrkperlevel}{number, default 1: the number of markers to select at
each of the levels in Rlevels}

\item{stat}{the name of one of the statistics columns in Rstats}
}
\value{
a selection from Rstats, in ascending order of column stat.
The first <mrkperlevel> markers with stat above each Rlevel are returned;
duplicated selections are removed (so the number of returned markers may be
less than length(Rlevels) * mrkperlevel)
}
\description{
Select markers based on their R statistics, for studying
the relation between R level and marker quality
}
\details{
The return value of this function is intended to be used for
studying the relation between the marker quality and the value of the
chosen R statistic, e.g. by drawing XY-plots of each of the selected
markers. By finding a suitable threshold for the R statistic bas markers
could be excluded from evaluation by fitPoly, saving time because bad
markers take the longest to be scored and are then often rejected anyway.
}
