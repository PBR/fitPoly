% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitPolyTools.R
\name{correctDosages}
\alias{correctDosages}
\title{Check if dosage scores may have to be shifted}
\usage{
correctDosages(chk, scores, parent1, parent2, ploidy, ploidy2,
polysomic=TRUE, disomic=FALSE, mixed=FALSE, parentsScoredWithF1,
absent.threshold=0.04, outfile=NA)
}
\arguments{
\item{chk}{data frame returned by function checkF1 when called without
shiftmarkers}

\item{scores}{data.frame with scores as produced by the fitMarkers
function of package fitPoly; at least columns MarkerName, SampleName and geno
must be present, any other columns are ignored}

\item{parent1}{character vector with names of the samples of parent 1}

\item{parent2}{character vector with names of the samples of parent 2}

\item{ploidy}{The ploidy of parent 1 (must be even, 2 (diploid) or larger).}

\item{ploidy2}{The ploidy of parent 2. If omitted it is
assumed to be equal to ploidy.}

\item{polysomic}{if TRUE at least all polysomic segtypes are considered;
if FALSE these are not specifically selected (but if e.g. disomic is TRUE,
any polysomic segtypes that are also disomic will still be considered);
same as used in the call to checkF1 that generated data.frame chk}

\item{disomic}{if TRUE at least all disomic segtypes are considered (see
param polysomic); same as used in the call to checkF1 that generated
data.frame chk}

\item{mixed}{if TRUE at least all mixed segtypes are considered (see
param polysomic). A mixed segtype occurs when inheritance in one parent is
polysomic (random chromosome pairing) and in the other parent disomic (fully
preferential chromosome pairing); same as used in the call to checkF1 that
generated data.frame chk}

\item{parentsScoredWithF1}{single logical. TRUE means that parents are scored
in the same experiment and the same fitPoly run as the F1, else FALSE.
If missing and ploidy2==ploidy, TRUE is assumed.
If FALSE, parental scores will not be shifted along with the F1 scores.}

\item{absent.threshold}{the threshold for the fraction of ALL samples
that has the dosage that is assumed to be absent due to mis-fitting of
fitPoly; should be at least the assumed error rate of the fitPoly scoring
assuming the fitted model is correct}

\item{outfile}{file name to which the result is written. If NA no file
is written.}
}
\value{
a data frame with columns
\itemize{
\item{markername}
\item{segtype: the bestfit (not bestParentfit!) segtype from chk}
\item{parent1, parent2: the consensus parental dosages; possibly
low-confidence, so may be different from those reported in chk}
\item{shift: -1, 0 or 1: the amount by which this marker should be shifted}
}
The next fields are only calculated if shift is not 0:
\itemize{
\item{fracNotOk: the fraction of ALL samples that are in the dosage
(0 or ploidy) that should be empty if the marker is indeed shifted.}
\item{parNA: the number of parental dosages that is missing (0, 1 or 2)}
}
}
\description{
fitPoly sometimes uses a "shifted" model to assign dosage
scores (e.g. all samples are assigned a dosage one higher than the true
dosage). This happens mostly when there are only few dosages present
among the samples. This function checks if a shift of +/-1 is possible.
}
\details{
A shift of -1 (or +1) is proposed when (1) the fraction of all
samples with dosage 0 (or ploidy) is below absent.threshold, (2) the
bestfit (not bestParentfit!) segtype in chk has one empty dosage on the
low (or high) side and more than one empty dosage at the high (or low) side,
and (3) the shifted consensus parental dosages do not conflict with the
shifted segregation type.\cr
The returned data.frame (or a subset, e.g. based on the values in the
fracNotOk and parNA columns) can serve as parameter shiftmarkers in a
new call to checkF1.\cr
Based on the quality scores assigned by checkF1 to
the original and shifted versions of each marker the user can decide if
either or both should be kept. A data.frame combining selected rows
of the original and shifted versions of the checkF1 output (which may
contain both a shifted and an unshifted version of some markers) can then be
used as input to compareProbes or writeDosagefile.
}
