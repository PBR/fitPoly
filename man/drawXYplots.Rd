% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitPolyTools.R
\name{drawXYplots}
\alias{drawXYplots}
\title{Draws a series of pages, each with 6 XY-plots showing allele signals
and assigned dosages}
\usage{
drawXYplots(dat, markers=NA, out, genocol="grey", pch=1, cex=1,
sel.samples=as.character(unique(XYdat$SampleName)), omit.pch=".",
omit.col=c(rep("black", length(genocol)-1), "grey"), sample.groups=list(),
groups.col="black", groups.pch=1, groups.cex=1, groups.rnd=FALSE,
R.col="black", R.lty=1, drawRthresholds=FALSE,
Rthreshold.param=c(0.95, 0.5, 0), ploidy)
}
\arguments{
\item{dat}{a data.frame with at least columns MarkerName, SampleName, X and Y;
column geno (if present) is also used}

\item{markers}{either a vector with names of markers to plot or a
data.frame with at least column MarkerName containing the names of the markers
to plot; default NA means all markers in dat. If markers is a data.frame that
also has a column shift, the geno values in dat will be shifted accordingly.
This allows to use (a selection from) the output of checkF1 (with or without
parameter shiftmarkers) as input for drawXYplots.}

\item{out}{base for filenames of output, will be extended with
_pagenumber.png; may include a path, but the directory where the plot files
are to be saved must already exist}

\item{genocol}{a vector of color values to be used for plotting
the sel.samples according on their geno (dosage) value; if only one value is
given (default) all samples are plotted in that color}

\item{pch}{the plot character to plot the sel.samples; default 1 is an
open circle}

\item{cex}{the relative size of the sample symbols}

\item{sel.samples}{character vector (not a factor) with the names of the
important samples: those that must be plotted in colors genocol and symbol pch}

\item{omit.pch}{the plot character to use for the other samples, default a dot}

\item{omit.col}{vector of two colors to use for the other samples; the first
color is used for sample with a geno (dosage) value, the second color for
unscored samples}

\item{sample.groups}{a list specifying samples to be highlighted in a
different color and/or symbol and/or size. For each group of sample the list
has one vector of sample names; the list may also be empty}

\item{groups.col}{a vector or color values, one for each item (vector of
sample names) in sample.groups; recycled if shorter than sample.groups}

\item{groups.pch}{a vector of plot symbols, one for each item in
sample.groups; recycled if shorter than sample.groups}

\item{groups.cex}{a vector of relative symbol sizes, one for each item in
sample.groups; recycled if shorter than sample.groups}

\item{groups.rnd}{FALSE (default) or TRUE. If FALSE, all samples in
sample.groups are drawn in the order in which they appear in XYdat; if TRUE
they are drawn in a random order. Note that the samples are never drawns in
group order (except if the samples are already in group order in XYdat and
groups.rnd=FALSE).}

\item{R.col}{a vector of color values for drawing R thresholds, one for
each value of Rthreshold.param; recycled if needed}

\item{R.lty}{a vector of line types for drawing the R thresholds, one for
each value of Rthreshold.param; recycled if needed}

\item{drawRthresholds}{FALSE (default) or TRUE: whether R thresholds should be
drawn}

\item{Rthreshold.param}{either a list of vectors each of length 3 (for
multiple R thresholds) or one vector of length 3 (for one R threshold). Each
vector defines one R threshold and is based on a specified quantile of the
distribution of R values for the current marker. The first number in each
vector is the R quantile, the second is a number to multiply that R quantile
with, and the third is the minimum value of the result. The default of
c(0.95, 0.5, 0) means that 0.5 * the 95% quantile of R is used (which is
always higher than 0, the minimum result). This is often a good cut-off
value to discard samples, or to signal markers with many samples below that
value.}

\item{ploidy}{a single integer specifying the ploidy, only needed if
dat contains a column geno}
}
\value{
The function produces a series of pages with plots and returns NULL
}
\description{
Draws 6 XY-plots per page for a series of markers; each XY-plot is
drawn by function XY_plot
}
