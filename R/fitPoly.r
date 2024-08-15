# package fitPoly ***************************************************************
#
#    fitPoly is the successor of fitTetra. It is an R package for assigning
#    polyploid genotype scores for bi-allelic marker assays.
#    (C) 2010-2024 Roeland E. Voorrips, Gerrit Gort, Alejandro Therese Navarro
#    Wageningen University and Research
#    Version 4.0, August 2024
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
#    or obtain one at http://www.gnu.org
#
#    contact: Alejandro Therese Navarro
#    e-mail: alejandro.theresenavarro@wur.nl
#    address: A. Therese Navarro
#             Wageningen University and Research - Plant Breeding
#             P.O. Box 16
#             6700 AA Wageningen
#             The Netherlands
#
# ******************************************************************************

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# There is a problem in RStudio (Version 0.99.902) when using
# Build - Check Package:
# When testing the examples it generates an error saying
# "There is no package called ...".
# A workaround is to first Build a source package, install it, and then
# run Build - Check package
# Another similar issue is found if the address of the package is too long,
# int that case devtools::check() does not correctly identify the package. 
# The solution is to just move the package to a directory closer to the base 
# for example C:/temp in windows, and then check the package.
# Packages required for building, checking and testing:
# roxygen2, devtools, testthat
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#'fitPoly: a package for assigning dosage scores based on SNP array data
#'
#'fitPoly (an evolved version of package fitTetra) fits mixture models to
#'the distribution of intensity ratios Y/(X+Y) (where X and Y are the
#'intensities of the signals produced by the A and B alleles of bi-allelic
#'markers) and uses these to assign genotypes (dosages).
#'The main differences compared with fitTetra are that it can handle
#'any ploidy level, and multiple populations that can be either F1 populations
#'(and their parents) or panels of accessions. There are also improvements
#'in accuracy, speed and the possibility to use prior dosage information.
#'
#'@docType package
#'@name fitPoly
NULL

#'@importFrom foreach foreach %do% %dopar%
#'@importFrom grDevices dev.cur dev.list dev.off pdf png svg
#'@importFrom graphics axis axTicks barplot box close.screen hist layout lines
#'            par plot points screen split.screen text
#'@importFrom stats dbinom dhyper dnorm integrate kmeans lm
#'            nls predict sd weighted.mean
#'@importFrom utils read.table write.table
#'

#the following line serves to stop devtools::check complaining about
#some identifiers that do exist but it doesn't find because they are
#loaded from RData files:
utils::globalVariables(names=c("batch", "batdat", "batchlog",
                               "batchallmodeldata", "batchmodeldata",
                               "batchscores", "batchdiploscores"),
                       add=TRUE)


#'@title Small fitPoly input datasets for testing and examples
#'@description A list with small datasets of four different
#'ploidy levels for testing and examples
#'@details list fitPoly_data contains the following items:
#'\itemize{
#'  \item ploidy2: a diploid dataset with only the SNP array signal ratios
#'  \item ploidy3: a triploid dataset with in addition to the signal ratios
#'                 two data.frames specifying the population structure
#'  \item ploidy4: a tetraploid dataset similar to the triploid dataset and
#'                 additionally prior dosage information of the F1 population parents
#'                 and of a few other samples
#'  \item ploidy6: a hexaploid dataset similar to the tetraploid dataset and
#'                 additionally the 7 starting means for some of the markers
#'}
#'Each of the items contains one or more elements, postfixed by 2x, 3x, 4x or 6x
#'depending on the ploidy:
#'\itemize{
#'  \item dat: a data.frame with at least columns MarkerName, SampleName and ratio
#'  with the signal ratios to be analyzed
#'  \item pop: a data.frame with columns SampleName and Population, specifying the
#'  population to which each sample belongs
#'  \item pop.par: a matrix specifying what are the parents of each population (if any)
#'  \item parPriors: a data.frame specifying prior known allele dosages for the F1 parents
#'  \item sampPriors: a data.frame specifying the prior known dosages for other samples
#'  \item startmeans: a data.frame with prior known means for the (ploidy+1) mixture model
#'  components
#'}
#'In addition the ploidy6 component has elements pop and pop.parents (no suffix)
#'which are equivalent to pop6x and pop.par6x, in the format required
#'by function codomMarker.
#'@name fitPoly_data
#'@usage data(fitPoly_data)
#'@docType data
NULL


# fitMarkers **************************************************
# This is a user function that calls fitOneMarker for a series of markers
# and saves the tabular and log output to files.
# *******************************************************************

getBatchsize <- function(mrkcount, ncores) {
  #an efficient target batchsize is around 500, but we vary it
  #to avoid processors standing idle too much with small data sets.
  #Also very small batchsizes are unwanted as the time to select a batch from
  #a large data.frame plus the combination of the batchresults takes more time
  #than processing a single marker.
  batchsize <- 500
  if (ncores > 1) {
    #cyc: how many cycles per core are needed
    cyc <- ceiling(mrkcount / (ncores * batchsize * 1.1))
    batchsize <- ceiling(mrkcount / (ncores * cyc))
    batchsize[batchsize < 3] <- 3
  } else batchsize <- min(batchsize, mrkcount) #only for the sMM message
  batchsize
} #getBatchsize

#'@title Function to fit mixture models for series of markers and save the
#'results to files
#'
#'@description This is the main function that calls fitOneMarker
#'for a series of markers and saves the tabular, graphical and log output to
#'files. Most of the arguments are identical to those of fitOneMarker and
#'are directly passed through.
#'
#'@usage fitMarkers(ploidy, markers=NA, data, diplo=NULL, select=TRUE,
#'diploselect=TRUE, pop.parents=NULL, population=NULL, parentalPriors=NULL,
#'samplePriors=NULL, startmeans=NULL, maxiter=40, maxn.bin=200, nbin=200,
#'sd.threshold=0.1, p.threshold=0.9, call.threshold=0.6, peak.threshold=0.85,
#'try.HW=TRUE, dip.filter=1, sd.target=NA,
#'filePrefix, rdaFiles=FALSE, allModelsFile=FALSE,
#'plot="none", plot.type="png", ncores=1)
#'
#'@param ploidy The ploidy level, 2 or higher: 2 for diploids, 3 for triploids
#'etc.
#'@param markers NA or a character or numeric vector specifying the markers to be
#'fitted. If a character vector, names should match the MarkerName column of
#'data; if numeric, the numbers index the markers based on the alphabetic order
#'of the MarkerNames in data.
#'@param data A data frame with the polyploid samples, with (at least) columns
#'MarkerName, SampleName and ratio, where ratio is the Y-allele signal
#'divided by the sum of the X- and Y-allele signals: ratio == Y/(X+Y)
#'@param diplo NULL or a data frame like data, with the diploid samples and (a subset
#'of) the same markers as in data. Genotypic scores for diploid samples are
#'calculated according to the best-fitting model calculated for the polyploid
#'samples and therefore may range from 0 (nulliplex) to <ploidy>, with the
#'expected dosages 0 and <ploidy> for the homozygotes and <ploidy/2> for the
#'heterozygotes.\cr
#'diplo can also be used for any other samples that need to be
#'scored, but that should not affect the fitted models.
#'@param select A logical vector, recycled if shorter than nrow(data):
#'indicates which rows of data are to be used (default TRUE, i.e. keep all rows)
#'@param diploselect A logical vector like select, matching diplo instead of data
#'@param pop.parents NULL or a data.frame specifying the population structure. The
#'data frame has 3 columns: the first containing population ID's, the 2nd and 3rd
#'with the population ID's of the parents of these populations (if F1's) or NA
#'(if not). The population ID's should match those in parameter population. If
#'pop.parents is NULL all samples are considered to be in one population, and
#'parameter population should be NULL (default).
#'@param population NULL or a data.frame specifying to which population each
#'sample belongs. The data frame has two columns, the first containing
#'the SampleName (containing all SampleNames occurring in data),
#'the second column containing population ID's that match pop.parents. In both
#'columns NA values are not allowed. Parameters pop.parents and population
#'should both be NULL (default) or both be specified.
#'@param parentalPriors NULL or a data frame specifying the prior dosages for
#'the parental populations. The data frame has one column MarkerName
#'followed by one column for each F1 parental population. Column names (except
#'first) are population ID's matching the parental populations in pop.parents.
#'In case there is just one F1 population in pop.parents, it is possible to
#'have two columns for both parental populations instead of one (allowing two
#'specify two different prior dosages); in that case both columns for each
#'parent have the same caption. Each row specifies the priors for
#'one marker. The contents of the data frame are dosages, as integers from 0
#'to <ploidy>; NA values are allowed.\cr
#'Note: when reading the data frame with read.table or read.csv, set
#'check.names=FALSE so column names (population ID's) are not changed.
#'@param samplePriors NULL or a data.frame specifying prior dosages for individual
#'samples. The first column called MarkerName is followed by one column per
#'sample; not all samples in data need to have a column here, only
#'those samples for which prior dosages for one or more markers are available.
#'Each row specifies the priors for one marker. The contents of the data frame
#'are dosages, as integers from 0 to <ploidy>; NA values are allowed.\cr
#'Note: when reading the data frame with read.table or read.csv, set
#'check.names=FALSE so column names (population ID's) are not changed.
#'@param startmeans NULL or a data.frame specifying the prior means of
#'the mixture distributions. The data frame has one column MarkerName,
#'followed by <ploidy+1> columns with the prior ratio means on the original
#'(untransformed) scale. Each row specifies the
#'means for one marker in strictly ascending order (all means NA is allowed, but
#'markers without start means can also be omitted).
#'@param maxiter A single integer, passed to CodomMarker, see there for explanation
#'@param maxn.bin A single integer, passed to CodomMarker, see there for explanation
#'@param nbin A single integer, passed to CodomMarker, see there for explanation
#'@param sd.threshold The maximum value allowed for the (constant) standard
#'deviation of each peak  on the arcsine - square root transformed scale,
#'default 0.1. If the optimal model has a larger standard deviation the marker
#'is rejected. Set to a large value (e.g. 1) to disable this filter.
#'@param p.threshold The minimum P-value required to assign a genotype (dosage)
#'to a sample; default 0.9. If the P-value for all possible genotypes is less
#'than p.threshold the sample is assigned genotype NA. Set to 1 to disable
#'this filter.
#'@param call.threshold The minimum fraction of samples to have genotypes
#'assigned ("called"); default 0.6. If under the optimal model the fraction of
#'"called" samples is less than call.threshold the marker is rejected. Set to 0
#'to disable this filter.
#'@param peak.threshold The maximum allowed fraction of the scored samples that
#'are in one peak; default 0.85. If any of the possible genotypes (peaks in the
#'ratio histogram) contains more than peak.threshold of the samples the marker
#'is rejected (because the remaining samples offers too little information for
#'reliable model fitting).
#'@param try.HW Logical: if TRUE (default), try models with and without a
#'constraint on the mixing proportions according to Hardy-Weinberg equilibrium
#'ratios. If FALSE, only try models without this constraint. Even when the HW
#'assumption is not applicable, setting try.HW to TRUE often still leads to
#'a better model. For more details on how try.HW is used see the Details
#'section of function fitOneMarker.
#'@param dip.filter if 1 (default), select best model only from models
#'that do not have a dip (a lower peak surrounded by higher peaks: these are not
#'expected under Hardy-Weinberg equilibrium or in cross progenies). If all
#'fitted models have a dip still the best of these is selected. If 2, similar,
#'but if all fitted models have a dip the marker is rejected. If 0, select best
#'model among all fitted models, including those with a dip.
#'@param sd.target If the fitted standard deviation of the peaks on the
#'transformed scale is larger than sd.target a penalty is given (see Details
#'section of function fitOneMarker);
#'default NA i.e. no penalty is given.
#'@param filePrefix partial file name, possibly including an absolute or
#'relative file path. filePrefix must always be specified.
#'All output files will have filePrefix prefixed to their name so it is clear
#'they are all derived from the same call to fitMarkers.
#'If filePrefix includes a file path all output files
#'will be saved there; if a filePrefix is specified that does not include a
#'a path the output will be saved in the working directory.
#'@param rdaFiles logical, default FALSE. The tabular output (scorefile,
#'diploscorefile, modelfile, allmodelsfile) is saved as tab-separated text files
#'with extension .dat or as an .RData file if this parameter is FALSE or TRUE
#'respectively.
#'@param allModelsFile logical, default FALSE. If TRUE an allmodelsfile is saved
#'with all models that have been tried for each marker; also the log file will
#'contain a few lines for each marker. This information is mostly useful
#'for debugging and locating problems.
#'@param plot String, "none" (default), "fitted" or "all". If "fitted" a plot
#'of the best fitting model and the assigned genotypes is saved with filename
#'<marker number><marker name>.<plot.type>, preceded by "rejected_" if the
#'marker was rejected. If "all", small plots of all models are saved to files
#'(8 per file) with filename
#'<"plots"><marker number><marker name><pagenr>.<plot.type> in addition to the
#'plot of the best fitting model.
#'@param plot.type String, "png" (default), "emf", "svg" or "pdf". Indicates
#'format for saving the plots.
#'@param ncores The number of processor cores to use for parallel processing,
#'default 1. Specifying more cores than available may cause problems.
#'Note that the implementation under Windows involves duplicating the input data
#'(under Linux that does not happen, nor under Windows if ncores=1), so if
#'under Windows memory size is a problem it would be better to run several
#'R instances simultaneously, each with ncores=1, each processing part of the
#'data.
#'@return NULL. The result of fitMarkers is a set of output files.
#'@details fitMarkers calls fitOneMarker for all markers specified by
#'parameter markers. The markers are processed in batches; the number of markers
#'per batch is printed to the console when fitMarkers is started. If
#'multiple cores are used the batches are processed in parallel.\cr
#'During the processing a series of RData files (2 for each batch) is saved in the
#'directory specified in filePrefix. At the end these are combined into the required
#'output files and then deleted.
#'If something goes wrong at any stage, the files for the completed batches are
#'still available and can be combined manually, avoiding the need to re-run the
#'process for the completed batches.
#'The output files consist of:
#'\itemize{
#'\item{<filePrefix>.log: a logfile containing several lines listing the
#'input parameters. If parameter allModelsFile is TRUE the logfile also
#'contains several text lines per marker, corresponding to component "log"
#'in the result of fitOneMarker}
#'\item{<filePrefix>_scores.dat (or .RData) a file containing one line per
#'polyploid sample for every marker that could be fitted, corresponding to
#'component "scores" in the result of fitOneMarker}
#'\item{<filePrefix>_diploscores.dat (or .RData) a file containing one line per
#'diploid sample for every marker that could be fitted, corresponding to
#'component "diploscores" in the result of fitOneMarker. This file is only produced
#'if parameter diplo is not missing}
#'\item{<filePrefix>_models.dat (or .RData) a file containing one line per
#'marker, corresponding to component "modeldata" in the result of fitOneMarker: the
#'selected model for each marker, with several statistics}
#'\item{<filePrefix>_allmodels.dat (or .RData) as the models file, but
#'containing all models fitted for each marker, not only the selected model,
#'marker, corresponding to component "allmodeldata" in the result of fitOneMarker.
#'This file is only produced if parameter allModelsFile is TRUE}
#'}
#'Additionally, if plot != "none", plot files are generated in directory
#'<filePrefix>_plots
# NOTE that in the examples we should not use ncores > 1!
# This results in hard-to-understand errors in devtools::check
#'@examples
#' \donttest{
#'  # These examples run for a total of about 55 sec.
#'  # All output files are saved in tempdir() and subdirectories of it.
#'
#'  data(fitPoly_data)
#'
#'  # tetraploid, with no populations and with sample prior dosages
#'  fitMarkers(ploidy=4, data=fitPoly_data$ploidy4$dat4x,
#'                   samplePriors=fitPoly_data$ploidy4$sampPriors4x,
#'                   filePrefix=paste0(tempdir(),"/4xA"),
#'                   allModelsFile=TRUE,
#'                   plot="fitted")
#'
#'  # tetraploid, with specified populations and parental and sample prior dosages
#'  fitMarkers(ploidy=4, data=fitPoly_data$ploidy4$dat4x,
#'                   population=fitPoly_data$ploidy4$pop4x,
#'                   pop.parents=fitPoly_data$ploidy4$pop.par4x,
#'                   parentalPriors=fitPoly_data$ploidy4$parPriors4x,
#'                   samplePriors=fitPoly_data$ploidy4$sampPriors4x,
#'                   filePrefix=paste0(tempdir(),"/4xB"),
#'                   allModelsFile=TRUE,
#'                   plot="fitted")
#'
#'  # hexaploid, no populations or prior information
#'  fitMarkers(ploidy=6, data=fitPoly_data$ploidy6$dat6x,
#'                   filePrefix=paste0(tempdir(),"/6xA"),
#'                   allModelsFile=TRUE,
#'                   plot="fitted")
#'
#'  # hexaploid, with specified populations, prior dosages of parents and other samples
#'  # and prior means of the mixture components
#'  fitMarkers(ploidy=6, data=fitPoly_data$ploidy6$dat6x,
#'                   population=fitPoly_data$ploidy6$pop6x,
#'                   pop.parents=fitPoly_data$ploidy6$pop.par6x,
#'                   startmeans=fitPoly_data$ploidy6$startmeans6x,
#'                   parentalPriors=fitPoly_data$ploidy6$parPriors6x,
#'                   samplePriors=fitPoly_data$ploidy6$sampPriors6x,
#'                   filePrefix=paste0(tempdir(),"/6xB"),
#'                   plot="fitted")
#' }
#'
#'@export
fitMarkers <- function(
  ploidy,
  markers=NA, #not NULL
  data, diplo=NULL,
  select=TRUE, diploselect=TRUE, # by default all samples of polyploids and diploids selected
  pop.parents=NULL,
  population=NULL,
  parentalPriors=NULL,
  samplePriors=NULL,
  startmeans=NULL,
  maxiter=40, maxn.bin=200, nbin=200,
  sd.threshold=0.1, p.threshold=0.9,
  call.threshold=0.6, peak.threshold=0.85,
  try.HW=TRUE, dip.filter=1,
  sd.target=NA, #not NULL
  filePrefix, rdaFiles=FALSE, allModelsFile=FALSE,
  plot="none", plot.type="png",
  ncores=1) {
  sMMnow <- format(Sys.time(), format='%Y%m%d-%H%M%S')
  pars <- as.list(match.call()) #parameters of call to this function
  if (missing(filePrefix)) stop("filePrefix must be a valid (path +) start of filename")
  fpf <- checkFilePrefix(filePrefix)
  fext <- if (rdaFiles) ".RData" else ".dat"
  logfile <- paste(fpf, ".log", sep="")
  modelfile <- paste(fpf, "_models", fext, sep="")
  if (allModelsFile) {
    allmodelsfile <- paste(fpf, "_allmodels", fext, sep="")
  } else allmodelsfile <- ""
  scorefile <- paste(fpf, "_scores", fext, sep="")
  diploscorefile <- paste(fpf, "_diploscores", fext, sep="")
  #plot.dir will be calculated based on fpf if needed

  # check ploidy:
  if (!(class(ploidy)[1] %in% c("numeric", "integer")) ||
      length(ploidy) != 1 ||
      ploidy < 2) {
    stop("fitMarkers: ploidy must be  >= 2")
  }

  # check (diplo)data and select
  # and get all markernames, samplenames and diplonames:
  # marker numbers refer to markers actually present in data,
  # in the alphabetic order according to the current OS and R version
  checkInputData(data)
  select <- checkSelect(select, data)
  markernames <- getSortedLevels(data$MarkerName)
  samplenames <- getSortedLevels(data$SampleName)
  IntendedSamples <- getIntendedSamples(select, data)
  if (!is.null(diplo) && !is1NA(diplo)) {
    checkInputData(diplo)
    diploselect <- checkSelect(diploselect, diplo)
    diplonames <- getSortedLevels(diplo$SampleName)
  } else diploscorefile <- ""

  # check population structures:
  origpopstruct <- getOrigPopstruct(samplenames=samplenames,
                                    pop.parents=pop.parents,
                                    population=population,
                                    parentalPriors=parentalPriors,
                                    try.HW=try.HW,
                                    ploidy=ploidy)
  if (is.character(origpopstruct))
    stop(paste("fitMarkers: invalid population structure:", origpopstruct))

  if (!is.null(startmeans) && !is1NA(startmeans)) checkStartmeans(ploidy, startmeans)
  if (!is.null(samplePriors) && !is1NA(samplePriors))
    checkSamplePriors(ploidy, samplePriors, origpopstruct)

  if (!(dip.filter %in% 0:2))
    stop("fitMarkers: dip.filter must be in 0:2")

  # arrange plotting output:
  plot <- checkPlot(plot, plot.type, plot.dir=NA, fpf) #plot.dir not NULL here!
  if (plot[3] != "") message(plot[3]) #cat(paste(plot.type[3], "\n", sep=""))
  plot.type <- plot[[2]]
  plot.dir <- plot[[4]]
  plot <- plot[[1]]

  # check markers:
  markers <- markers[!is.na(markers)]
  if (length(markers) == 0) markers <- 1:length(markernames) else {
    if (is.character(markers)) {
      mrknrs <- match(markers, markernames)
      if (anyNA(mrknrs))
        stop("fitMarkers: some marker names in markers don't occur in data")
      markers <- mrknrs
    } else if (!all(markers %in% 1:length(markernames))) {
      stop("fitMarkers: invalid markers")
    }
  }

  sMMinfo <-
    list(markernames=markernames,
         samplenames=samplenames,
         IntendedSamples=IntendedSamples,
         origpopstruct=origpopstruct)
  if (exists("diplonames")) {
    sMMinfo$diplonames <- diplonames
    #sMMinfo$IntendedDiplosamples <- IntendedDiplosamples
  }

  suppressWarnings( {
    file.remove(scorefile)
    file.remove(diploscorefile)
    file.remove(allmodelsfile)
    file.remove(modelfile)
    file.remove(logfile)
    }
  )

  suppressPackageStartupMessages({ #suppress message from foreach

  #NOTE: if fitPoly is sourced instead of used as a package,
  #library(foreach) is needed before calling fitMarkers.
  #The reason is that %do% and %dopar% cannot be prefixed with 'foreach::'

  #There is a problem with foreach + doParallel under Windows, see
  #http://stackoverflow.com/questions/34653567/doparallel-foreach-inconsistently-inherits-objects-from-parent-environment-e?rq=1
  #date: 9-2-2016, R version: 3.2.3
  #The following adresses this, but apparently at the cost of duplicating the
  #input data for each core:
  #(duplication does not occur even under Windows if ncores==1)
  foreachExports <- NULL
  if (.Platform$OS.type == "windows" && ncores > 1)
    foreachExports <- setdiff(ls(envir=globalenv()),
                              c("parentalPriors", "samplePriors", "startmeans"))
  #  If the three last are not removed, under Windows with ncores>1 we get:
  #  Warning message:
  #  In e$fun(obj, substitute(ex), parent.frame(), e$data) :
  #    already exporting variable(s): parentalPriors, samplePriors, startmeans
  #  But if we don't export anything (in the same situation) we get:
  #  Error in { : task 1 failed - "could not find function "processBatch""

  if (length(ncores) != 1 || is.na(ncores) || ncores < 2) ncores <- 1
  if (ncores > 1 && !requireNamespace("doParallel", quietly=TRUE)) {
    ncores <- 1
    message("fitMarkers: package doParallel not available, using only 1 core")
  }
  batchsize <- getBatchsize(length(markers), ncores)
  message(paste("fitMarkers: batchsize =", batchsize))
  batchcount <- ceiling(length(markers) / batchsize)
  # batches give a speed increase because the selection of a subset of rows
  # from a large data frame takes very much longer than from a small data frame.
  # Also batches are efficiently parallellized.
  if (logfile != "") {
    write (paste("fitMarkers started at", sMMnow), file=logfile)
    write (paste("ploidy =", ploidy), file=logfile, append=TRUE)
    if ("data" %in% names(pars))
      write (paste("data frame passed as data =", pars$data), file=logfile, append=TRUE)
    write (paste("sample count in data =", length(samplenames)), file=logfile, append=TRUE)
    if (exists("diplonames")) {
      write (paste("data frame passed as diplo =", pars$diplo), file=logfile, append=TRUE)
      write (paste("sample count in diplo =", length(diplonames)), file=logfile, append=TRUE)
    }
    write (paste("marker count in data =", length(markernames)), file=logfile, append=TRUE)
    write (paste("marker count to score =", length(markers)), file=logfile, append=TRUE)
    write (paste("maxiter =", maxiter), file=logfile, append=TRUE)
    write (paste("maxn.bin =", maxn.bin), file=logfile, append=TRUE)
    write (paste("nbin =", nbin), file=logfile, append=TRUE)
    write (paste("sd.threshold =", sd.threshold), file=logfile, append=TRUE)
    write (paste("p.threshold =", p.threshold), file=logfile, append=TRUE)
    write (paste("call.threshold =", call.threshold), file=logfile, append=TRUE)
    write (paste("peak.threshold =", peak.threshold), file=logfile, append=TRUE)
    write (paste("try.HW =", try.HW), file=logfile, append=TRUE)
    write (paste("dip.filter =", dip.filter), file=logfile, append=TRUE)
    write (paste("sd.target =", sd.target), file=logfile, append=TRUE)
    write (paste("ncores =", ncores), file=logfile, append=TRUE)
    write (paste("batchsize =", batchsize), file=logfile, append=TRUE)
  }

  if (ncores > 1) {
    doParallel::registerDoParallel(cores=ncores)
    stats <-
      foreach::foreach(batch = 1:batchcount, .combine='rbind',
                       .multicombine=TRUE, .maxcombine=batchcount+1,
                       .export=foreachExports) %dopar% {
      processBatch(
        batch=batch, batchsize=batchsize, batchcount=batchcount,
        markers=markers,
        data=data, diplo=diplo,
        select=select, diploselect=diploselect,
        parentalPriors=parentalPriors,
        samplePriors=samplePriors,
        startmeans=startmeans,
        maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
        sd.threshold=sd.threshold,
        p.threshold=p.threshold,
        call.threshold=call.threshold,
        peak.threshold=peak.threshold,
        try.HW=try.HW,
        dip.filter=dip.filter,
        sd.target=sd.target,
        plot=plot, plot.type=plot.type, plot.dir=plot.dir,
        ploidy=ploidy,
        savelog=TRUE,
        savediplo=is.data.frame(diplo),
        saveallmodels=allmodelsfile!="",
        filePrefix=fpf,
        sMMnow=sMMnow,
        sMMinfo=sMMinfo)
    } # foreach batch
  } else {
    #ncores == 1
    #no parallel backend, and %do% instead of %dopar% and no .export
    stats <-
      foreach::foreach(batch = 1:batchcount, .combine='rbind',
                       .multicombine=TRUE, .maxcombine=batchcount+1) %do% {
      prbat <- processBatch(
        batch=batch, batchsize=batchsize, batchcount=batchcount,
        markers=markers,
        data=data, diplo=diplo,
        select=select, diploselect=diploselect,
        parentalPriors=parentalPriors,
        samplePriors=samplePriors,
        startmeans=startmeans,
        maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
        sd.threshold=sd.threshold,
        p.threshold=p.threshold,
        call.threshold=call.threshold,
        peak.threshold=peak.threshold,
        try.HW=try.HW,
        dip.filter=dip.filter,
        sd.target=sd.target,
        plot=plot, plot.type=plot.type, plot.dir=plot.dir,
        ploidy=ploidy,
        savelog=logfile!="",
        savediplo=diploscorefile!="",
        saveallmodels=allmodelsfile!="",
        filePrefix=fpf,
        sMMnow=sMMnow,
        sMMinfo=sMMinfo)
        #print(paste("batch", batch, "finished"))
        prbat # the line to be cbind-ed by foreach
    } # foreach batch
    #print("batch-loop finished")
  }
  }) #suppressPackageStartupMessages
  if (is.vector(stats)) stats <- matrix(stats, nrow=1) #happens if batchcount==1
  #concatenate and sort data of all batches:
  suppressWarnings(rm(batchscores, batchdiploscores,
                      batchallmodeldata, batchmodeldata, batchlog,
                      batdat))
  #we process the scores and models in two successive loops to
  #reduce memory requirements
  samplevels <- list()
  if (is.factor(data$SampleName)) {
    samplevels[[1]] <- levels(data$SampleName)
  } else samplevels[[1]] <- sort(unique(data$SampleName))
  if (is.data.frame(diplo)) {
    if (is.factor(diplo$SampleName)) {
      samplevels[[2]] <- levels(diplo$SampleName)
    } else samplevels[[2]] <- sort(unique(diplo$SampleName))
  }
  if (is.factor(data$MarkerName)) {
    mrklevels <- levels(data$MarkerName)
  } else mrklevels <- sort(unique(data$MarkerName))
  combineScoreBatches(ploidy+1, sMMnow, rowcounts=stats[,2:3, drop=FALSE],
                      outfnames=c(scorefile, diploscorefile), ncores,
                      mrklevels=mrklevels, samplevels,
                      filePrefix=fpf)
  combineModelBatches(ploidy+1, sMMnow, rowcounts=stats[,4:6, drop=FALSE],
                      outfnames=c(modelfile, allmodelsfile,
                                  ifelse(allModelsFile, logfile, "")),
                      ncores,
                      mrklevels=mrklevels,
                      pop.parents=origpopstruct$orig_pop.parents,
                      filePrefix=fpf)
  if (logfile != "")
    write (paste("fitMarkers finished at",
                 format(Sys.time(), format='%Y%m%d-%H%M%S')), file=logfile,
           append=TRUE)
    invisible(NULL)
} # fitMarkers


#'@title DEPRECATED: Function to fit mixture models for series of markers and save the
#'results to files
#'
#' @description This is the old name of the function `fitMarkers()`, the main function
#' of `fitPoly`. It is kept only for backwards-compatibility. The only difference between
#' the two is the default parameter of the genotyping probability threshold `p.threshold`, 
#' it is 0.99 in `saveMarkerModels()` (as was originally set) and it is 0.9 in the current
#' `fitMarkers()`. 
#' 
#' @param ... All parameters allowed in the function `fitMarkers()`. For a full description
#' see the help of that function.
#' @param p.threshold The minimum P-value required to assign a genotype (dosage)
#'to a sample; default 0.99 (very stringent). If the P-value for all possible genotypes is less
#'than p.threshold the sample is assigned genotype NA. Set to 1 to disable
#'this filter.
#'
#' @return See `fitMarkers()` documentation for a full description.
#' @export
#'
saveMarkerModels <- function(..., p.threshold = 0.99){
  
  warning("saveMarkerModels() has been deprecated in favour of fitMarkers().\nBe aware that default p.threshold of 0.99 is considered too stringent.")
  
  #Just here for backwards compatibility
  fitMarkers(p.threshold = p.threshold, ...)
}

processBatch <- function(
  batch, batchsize, batchcount,
  markers,
  data, diplo,
  select, diploselect,
  parentalPriors,
  samplePriors,
  startmeans,
  maxiter, maxn.bin, nbin,
  sd.threshold,
  p.threshold,
  call.threshold,
  peak.threshold,
  try.HW,
  dip.filter,
  sd.target,
  plot, plot.type, plot.dir,
  ploidy,
  savelog, savediplo, saveallmodels,
  filePrefix,
  sMMnow, sMMinfo) {
  #the processing of one batch in the foreach loop in fitMarkers
  #print(paste("batch =",batch,"; start_mem =",pryr::mem_used()))
  minmrkix <- (batch - 1) * batchsize + 1
  maxmrkix <- min(batch * batchsize, length(markers))
  batchdatarows <- data$MarkerName %in%
    sMMinfo$markernames[markers[minmrkix:maxmrkix]]
  batchdata <- data[batchdatarows,]
  batchselect <- select[batchdatarows]
  if (is.data.frame(diplo)) {
    batchdiplorows <- diplo$MarkerName %in%
      sMMinfo$markernames[markers[minmrkix:maxmrkix]]
    batchdiplo <- diplo[batchdiplorows, ]
    batchdiploselect <- diploselect[batchdiplorows]
  } else batchdiplo <- NULL
  batchmrknames <- sMMinfo$markernames[markers[minmrkix:maxmrkix]]
  if (is.null(parentalPriors) || is1NA(parentalPriors)) batchParentalPriors <- NULL else {
    batchParentalPriors <- parentalPriors[parentalPriors[,1] %in% batchmrknames,]
    if (nrow(batchParentalPriors) == 0) batchParentalPriors <- NULL
  }
  if (is.null(samplePriors) || is1NA(samplePriors)) batchSamplePriors <- NULL else {
    batchSamplePriors <- samplePriors[samplePriors[,1] %in% batchmrknames,]
    if (nrow(batchSamplePriors) == 0) batchSamplePriors <- NULL
  }
  if (is.null(startmeans) || is1NA(startmeans)) batchStartmeans <- NULL else {
    batchStartmeans <- startmeans[startmeans[,1] %in% batchmrknames,]
    if (nrow(batchStartmeans) == 0) batchStartmeans <- NULL
  }
  batchscores <- data.frame()
  batchdiploscores <- data.frame()
  batchmodeldata <- data.frame()
  batchallmodeldata <- data.frame()
  batchlog <- character(0)
  for (mrkix in minmrkix:maxmrkix) {
    mrk <- markers[mrkix]
    mrkresult <- fitOneMarker(
      ploidy=ploidy,
      marker=mrk,
      data=batchdata, diplo=batchdiplo,
      select=batchselect, diploselect=batchdiploselect,
      parentalPriors=batchParentalPriors,
      samplePriors=batchSamplePriors,
      startmeans=batchStartmeans,
      maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
      sd.threshold=sd.threshold,
      p.threshold=p.threshold,
      call.threshold=call.threshold,
      peak.threshold=peak.threshold,
      try.HW=try.HW,
      dip.filter=dip.filter,
      sd.target=sd.target,
      plot=plot, plot.type=plot.type, plot.dir=plot.dir,
      sMMinfo=sMMinfo) #list as expected by fitOneMarker

    #add the results of this marker to the batch results:
    if (is.data.frame(mrkresult$scores)) { # test for NA
      batchscores <- rbind(batchscores, mrkresult$scores)
    }
    #print(paste("scores added for mrk=", mrk))
    if (savediplo && is.data.frame(mrkresult$diploscores)) {
      batchdiploscores <- rbind(batchdiploscores, mrkresult$diploscores)
    }
    #print(paste("diploscores added for mrk=", mrk))
    if (saveallmodels && is.data.frame(mrkresult$allmodeldata)) {
      batchallmodeldata <- rbind(batchallmodeldata, mrkresult$allmodeldata)
    }
    batchmodeldata <- rbind(batchmodeldata, mrkresult$modeldata)
    if (savelog) batchlog <- c(batchlog, mrkresult$log)
    #print(paste("mrk=", mrk, "finished"))
  } # for mrk
  #write the batch results files:
  #we keep the (all)modeldata and the scores separate as they will become
  #very big after combining
  fname <- paste(sMMnow, "_batch",
                 padded(batch, batchcount),
                 ".RData", sep="")
  batdat <- list(batchscores=batchscores,
                 batchdiploscores=batchdiploscores)
  save(batdat, file=paste(filePrefix, "_", "sMMscores", fname, sep=""))
  batdat <- list(batchmodeldata=batchmodeldata,
                 batchallmodeldata=batchallmodeldata,
                 batchlog=batchlog)
  suppressWarnings(rm(batchmodeldata, batchallmodeldata, batchlog))
  for (i in 1:2) {
    batdat[[i]]$MarkerName <- factor(batdat[[i]]$MarkerName)
    batdat[[i]]$model <- factor(batdat[[i]]$model)
    batdat[[i]]$message <- factor(batdat[[i]]$message)
  }
  save(batdat, file=paste(filePrefix, "_", "sMMmodels", fname, sep=""))
  #let foreach save the number of rows of data:
  return(c(batch, nrow(batchscores), nrow(batchdiploscores),
    nrow(batdat[[1]]), nrow(batdat[[2]]), length(batdat[[3]])))
} #processBatch

combineScoreBatches <- function(ng, sMMnow, rowcounts, outfnames, ncores,
                                mrklevels, samplevels, filePrefix) {
  batchcount <- nrow(rowcounts)
  batchfnames <- paste(
    paste(filePrefix, "_","sMMscores", sMMnow, "_batch", sep=""),
    padded(1:batchcount, batchcount),
    ".RData", sep="")
  totrow <- colSums(rowcounts)
  firstline <- matrix(NA, nrow=batchcount+1, ncol=2)
  lst <- list() #will contain scores and diploscores
  for (i in 1:2) {
    if (outfnames[i] != "") {
      lst[[i]] <- makeScoresDF(totrow[i], ng, mrklevels, samplevels[[i]])
      firstline[,i] <- c(1, 1+cumsum(rowcounts[,i]))
    }
  }
  for (batch in 1:batchcount) {
    load(file=batchfnames[batch])
    #  loaded batdat: list with batchscores and possibly batchdiploscores
    for (i in 1:2) if (outfnames[i] != "" && nrow(batdat[[i]]) > 0) {
      lst[[i]][firstline[batch, i]:(firstline[batch+1, i] - 1),] <-
        batdat[[i]]
    }
    suppressWarnings(rm(batdat))
  } #for batch
  if (ncores > 1) for (i in 1:2) {
    if (outfnames[i] != "" && totrow[i] > 0)
      lst[[i]] <- lst[[i]][order(lst[[i]]$marker, lst[[i]]$SampleName),]
  }
  #for savedata we need actual variable names:
  scores <- lst[[1]]
  savedata(scores, outfnames[1])
  if (outfnames[2] != "") {
    diploscores <- lst[[2]]
    savedata(diploscores, outfnames[2])
  }
  file.remove(batchfnames)
} #combineScoreBatches

combineModelBatches <- function(ng, sMMnow, rowcounts, outfnames, ncores,
                                mrklevels, pop.parents, filePrefix) {
  batchcount <- nrow(rowcounts)
  batchfnames <- paste(
    paste(filePrefix, "_", "sMMmodels", sMMnow, "_batch", sep=""),
    padded(1:batchcount, batchcount),
    ".RData", sep="")
  totrow <- colSums(rowcounts)
  firstline <- matrix(NA_integer_, nrow=batchcount+1, ncol=3)
  lst <- list() #will contain modeldata, allmodeldata and log
  for (i in 1:3) {
    if (outfnames[i] != "") {
      if (i<3) {
        lst[[i]] <- makeModelsDF(totrow[i], ng, mrklevels, pop.parents)
        lst[[i]]$message <- factor(lst[[i]]$message) #to reduce memory space,
        #  but this means we must adjust the levels in the loop below
      } else lst[[i]] <- character(totrow[i])
      firstline[,i] <- c(1, 1+cumsum(rowcounts[,i]))
    }
  }
  for (batch in 1:batchcount) {
    load(file=batchfnames[batch])
    #  loaded batdat: list with batchmodeldata and possibly batchallmodeldata and/or batchlog
    for (i in 1:2) {
      if (outfnames[i] != "") {
        levels(lst[[i]]$message) <-
          c(levels(lst[[i]]$message),
          setdiff(levels(batdat[[i]]$message), levels(lst[[i]]$message)))
        lst[[i]][firstline[batch, i]:(firstline[batch+1, i] - 1),] <-
          batdat[[i]]
      }
    }
    if (outfnames[3] != "") {
      lst[[3]][firstline[batch, 3]:(firstline[batch+1, 3] - 1)] <- batdat[[3]]
    }
    suppressWarnings(rm(batdat))
  } #for batch

  for (i in 1:2) {
    if (outfnames[i] != "" && totrow[i] > 0) {
      if (ncores > 1)
        lst[[i]] <- lst[[i]][order(lst[[i]]$marker, lst[[i]]$m),]
    }
  }
  #for savedata we need actual variable names:
  modeldata <- lst[[1]]
  savedata(modeldata, outfnames[1])
  if (outfnames[2] != "") {
    allmodeldata <- lst[[2]]
    savedata(allmodeldata, outfnames[2])
  }
  if (outfnames[3] != "" ) {
    log <- lst[[3]]
    savedata(log, outfnames[3], append=TRUE)
  }
  file.remove(batchfnames)
} #combineModelBatches

makeScoresDF <- function(nrow, ng, mrklevels, samplevels) {
  rows0 <- nrow == 0
  if (rows0) nrow <- 1 #else construction of the data.frame fails
  df <- data.frame(
    marker=integer(nrow),
    MarkerName=factor(NA, levels=mrklevels),
    SampleName=factor(NA, levels=samplevels),
    #model=factor(NA, levels=getAllModelNames(ng)),
    #nsamp=NA_integer_,
    #select=NA_integer_,
    matrix(NA_real_, nrow=nrow, ncol=ng+1),
    maxgeno=NA_integer_,
    maxP=NA_real_,
    geno=NA_integer_
  )
  names(df)[4:(4+ng)] <- c("ratio", paste("P", 0:(ng-1), sep=""))
  if (rows0) df <- df[0,]
  df
} #makeScoresDF

makeModelsDF <- function(nrow, ng, mrklevels, pop.parents) {
  rows0 <- nrow == 0
  if (rows0) nrow <- 1 #else construction of the data.frame fails
  nm <- getResultnames(ng, pop.parents)
  df  <- data.frame(
    marker=rep(NA_integer_, nrow),
    MarkerName=factor(NA, levels=mrklevels),
    m=NA_integer_,
    model=factor(NA, levels=getAllModelNames(ng)),
    matrix(NA_integer_, nrow=nrow, ncol=5), #nsamp, nsel, npar, iter, dip
    matrix(NA_real_,  nrow=nrow, ncol=length(nm)-10),
    message=NA_character_, #we don't know all possible messages so we cannot
    #                       make it a factor (which would save on memory)
    stringsAsFactors=FALSE
  )
  names(df) <- nm
  if (rows0) df <- df[0,]
  df
} #makeModelsDF

checkFilePrefix <- function(fpf) {
  #fpf (file prefix) is a string that may have a (relative or absolute) path
  #plus the start of a valid filename
  #if the start of the filename is missing it is replaced by "sMM"
  #later, file names will be appended with _ as separator
  if (length(fpf) > 1 || is.na(fpf)) fpf <- ""
  fpf <- gsub("(^ +)|( +$)", "", fpf) #remove leading and trailing blanks
  if (substring(fpf, nchar(fpf), nchar(fpf)) %in% c("/", "\\")) {
      fpfdir <- fpf
      fpfname <- "sMM"
  } else {
    fpfdir <- dirname(fpf)
    if (fpfdir %in% c("", ".")) fpfdir <- "" else fpfdir <- paste(fpfdir, "/", sep="")
    fpfname <- gsub("(^ +)|( +$)", "", basename(fpf))
    if (fpfname == "") fpfname <- "sMM"
  }
  fpf <- paste(fpfdir, fpfname, sep="")
  if (!check.filename(
    paste(fpf,
          paste(sprintf("%x", sample(0:255, 8)), collapse=""), sep="_")))
    stop("filePrefix must be a valid (path +) filename")
  fpf
} #checkFilePrefix

check.filename <- function(filename, overwrite=TRUE) {
  #Check if a filename is valid for writing by trying to create the file
  #If the filename contains a path, result will only be TRUE if the whole
  #path already exists.
  #filename: a single text string to be interpreted as file (path +) name
  #overwrite: if TRUE (default) an existing file of that name will be deleted
  #           (if it is not protected by being locked, read-only, owned by
  #           another user etc)
  if (length(filename) != 1) return(FALSE)
  if (filename != gsub("(^ +)|( +$)", "", filename)) {
    #filename contains leading and trailing blanks)
    return(FALSE)
  }
  if (!overwrite && file.exists(filename)) return(FALSE)
  tryCatch(
{ suppressWarnings(
  if (file.create(filename)) {
    file.remove(filename)
    TRUE
  } else FALSE)
}, error = function(e) { FALSE }
  )
} #check.filename

savedata <- function(data, filename, append=FALSE) {
  #The purpose of this function is to save an  object under its own name;
  #if the function is called with data=<a value> then it is saved as an object
  #with name 'data'
  if (tolower(substr(tools::file_ext(filename), 1, 3)) == "rda") {
    #save data in an RData file:
    if (append) {
      stop(paste("savedata", filename,
                ": append not supported for RData files", sep=""))
    } else {
      pars <- as.list(match.call())
      if (is.name(pars[[2]])) {
        #savedata was called with data=<object name>;
        #first we create a local object with that name:
        command <- paste(as.character(pars[[2]]), "<- data")
        #print(paste("savedata: command=", command))
        commandexp <- parse(text=command)
        eval(commandexp)
        #we save the object under its original name:
        command <- paste("save(", as.character(pars[[2]]), ", file=filename)", sep="")
        #print(paste("savedata: command=", command))
        commandexp <- parse(text=command)
        eval(commandexp)
        #this might go wrong if this function is called with data=data
        #but that won't happen in fitMarkers
      } else {
        #savedata was called with data=<values>;
        #the values are saved as an object with the name 'data':
        save(data, file=filename)
      }
    } #RData, not append
  } else {
    #save data as a text file:
    if (is.data.frame(data)) {
      write.table(data, file=filename, col.names=!append, row.names=FALSE,
                     na="", quote=FALSE, sep="\t", append=append)
    } else {
      #not a data frame, we assume it is a vector:
      write(data, ncolumns=1, file=filename, append=append)
    }
  }
} #savedata

checkInputData <- function(data) {
  pars <- as.list(match.call()) #parameters of call to this function
  if (! is.name(pars[[2]]))
    #should never occur!
    stop("checkInputData may only be called with a named argument")
  if (!is.data.frame(data) ||
      nrow(data) == 0 ||
      length(which(names(data) == "MarkerName")) != 1 ||
      length(which(names(data) == "SampleName")) != 1 ||
      length(which(names(data) == "ratio")) != 1 )
    stop(paste(as.character(pars[[2]]), "is not a valid data frame"))
  if (any(data$ratio < 0, na.rm=TRUE) ||
      any(data$ratio > 1, na.rm=TRUE))
    stop(paste(as.character(pars[[2]]), "$ratio contains invalid values",
               sep=""))
} #checkInputData

checkSelect <- function(select, data) {
  pars <- as.list(match.call()) #parameters of call to this function
  if (!is.name(pars[[2]]))
    #should never occur!
    stop("checkSelect may only be called with a named argument")
  if (any(is.na(select)))
    stop(paste(as.character(pars[[2]]), "may not contain NA values"))
  if (!is.logical(select)) {
    if (!all(select %in% 0:1)) {
      stop(paste(as.character(pars[[2]]), "must be a logical vector"))
    } else select <- as.logical(select)
  }
  if (length(select) != nrow(data)) {
    select <- rep(select, times=ceiling(nrow(data)/length(select)))
    if (length(select) > nrow(data)) select <- select[1:nrow(data)]
  }
  select
} # checkSelect

#allNA <- function(x) { all(is.na(x)) }

getIntendedSamples <- function(select, data) {
  #intended samples are those that occur at least once in data and that are
  #not de-selected or have ratio set to NA for ALL markers
  #(the idea being that samples that are not relevant are omitted by
  #deleting their rows, de-selecting them for all markers, or setting all their
  #ratios to NA; e.g. samples with bad DNA, samples not relevant for the
  #current analysis etc.)
  missing <- !select | is.na(data$ratio)
  tap <- tapply(missing, data$SampleName, all) #all missing: TRUE or FALSE
  sort(names(tap)[!tap])
}

getSortedLevels <- function(x) {
  #get the sorted unique values of vector x;
  #is x is a character of numeric, sorted on the values,
  #if x is a factor, sorted on the levels instead of the values themselves
  result <- unique(x)
  if (is.factor(result)) result <- as.character(result)
  sort(result)
} #getSortedLevels

checkSamplePriors <- function(ploidy, samplePriors, origpopstruct) {
  #samplePriors: a data.frame with first column MarkerName followed by one
  #column for each sample that has priors; the sample names are the column names
  #Note: when reading the data frame with read.table, set check.names=FALSE
  #      so sample names are not changed
  if (!is.data.frame(samplePriors) ||
      length(samplePriors) < 2)
    stop("invalid samplePriors")
  if (names(samplePriors)[1] != "MarkerName")
    stop("name of 1st column of samplePriors must be MarkerName")
  if (anyNA(samplePriors$MarkerName) ||
      anyDuplicated(samplePriors$MarkerName) > 0)
    stop("no missing or duplicated MarkerNames allowed in samplePriors")
  for (i in 2:length(samplePriors)) {
    if (!is.numeric(samplePriors[,i])) stop("invalid samplePriors")
  }
  x <- as.matrix(samplePriors[, -1])
  xm <- suppressWarnings( max(x - trunc(x), na.rm=TRUE) )
  if (xm > 0) stop("all priors in samplePriors must be integers")
  suppressWarnings({
    if (!all(is.na(x)) &&
        (min(x, na.rm=TRUE) < 0 ||
         max(x, na.rm=TRUE) > ploidy))
      stop("invalid values in samplePriors")
  })
  if (!is.null(origpopstruct$orig_parPriorCols)) {
    #check if none of the samples in samplePriors is a parent
    samplepops <- origpopstruct$orig_population
    samplepops <-
      unique(samplepops[samplepops[,1] %in% names(samplePriors)[-1], 2])
    if (length(intersect(samplepops, origpopstruct$orig_parents)) > 0)
      stop("sample sets of samplePriors and parentalPriors overlap")
  }
} #checkSamplePriors

checkStartmeans <- function(ploidy, startmeans) {
  #startmeans: a data.frame with first column MarkerName followed by (ploidy+1)
  #columns for the means
  #On each row all means must be NA or none, and means must be in strictly
  #ascending order
  if (!is.data.frame(startmeans) ||
      length(startmeans) != ploidy+2)
    stop("invalid startmeans")
  if (names(startmeans)[1] != "MarkerName")
    stop("name of 1st column of startmeans must be MarkerName")
  if (anyNA(startmeans$MarkerName) || anyDuplicated(startmeans$MarkerName) > 0)
    stop("no missing or duplicated MarkerNames allowed in startmeans")
  for (i in 2:length(startmeans)) {
    if (!is.numeric(startmeans[,i])) stop("invalid startmeans")
  }
  x <- as.matrix(startmeans[, -1])
  suppressWarnings({
    if (min(x, na.rm=TRUE) < 0 ||
        max(x, na.rm=TRUE) > 1)
      stop("invalid values in startmeans")
  })
  NAsums <- rowSums(is.na(x))
  if (any(NAsums > 0 & NAsums < ploidy+1))
    stop("invalid startmeans: all or none of means must be NA")
  d <- diff(t(x)) #diff works within columns, so transpose x
  if (any(!is.na(d) & d <= 0))
    stop("invalid startmeans: means must be in increasing order")
} #checkStartmeans

getPolysomicSegr <- function(ploidy, ploidy2=ploidy) {
  #produce an array with 3 dimensions: dosages of parent 1, parent 2
  #and offspring, with for each combination of parental dosages the
  #polysomic segregation ratios (fractions adding to 1) in the progeny.
  #ploidy and ploidy2 each must be a positive even number, may be different
  #(not checked)

  polygamfrq <- function(ploidy, dosage) {
    #calculate the probabilities of the gamete dosages given the parental
    #ploidy and dosage
    gp <- ploidy / 2
    dhyper(0:gp, dosage, ploidy-dosage, gp)
  }

  popploidy <- (ploidy+ploidy2)/2
  result <- array(NA_real_, dim = c(ploidy+1, ploidy2+1, popploidy+1),
                  dimnames=list(P1=0:ploidy, P2=0:ploidy2, F1=0:popploidy))
  gp1 <- ploidy / 2 #ploidy of parent1 gametes
  gp2 <- ploidy2 / 2 #ploidy of parent2 gametes
  for (p1 in 0:ploidy) {
    gam1 <- polygamfrq(ploidy, p1)
    for (p2 in 0:ploidy2) {
      gam2 <- polygamfrq(ploidy2, p2)
      outmat <- outer(gam1, gam2)
      result[p1+1, p2+1,] <- tapply(outmat, col(outmat) + row(outmat) - 2, sum)
      #array of 1 row and 5 columns, with expected probabilities,
      #      with dosages (0..popploidy) as colnames
    }
  }
  result
} #getPolysomicSegr

addErrorProb <- function(segrArray, errorprob=0.04) {
  #segrArray: array as produced by getPolysomicSegr,
  #           or a vector of ng probabilities summing to 1
  #errorprob: the probability that a score is obtained randomly
  #result: an array or vector as segrArray with the probabilities combined
  #        with an error distribution

  d <- length(dim(segrArray))
  if (d == 0) ng <- length(segrArray) else ng <- dim(segrArray)[d]

  #simple version: uniform error distribution:
  #errorprob/ng + (1-errorprob) * segrArray

  #more realistic: the probability of a misscore is halved with each
  #dosage step from the true value:
  #error coeff:
  tmp <- 2^ -c((ng-1):1, 100, 1:(ng-1))
  errcoeff <- matrix(NA_real_, nrow=ng, ncol=ng)
  for (i in 1:ng) errcoeff[i,] <- tmp[(ng-i+1):(2*ng-i)]
  errcoeff <- (1/rowSums(errcoeff)) * errcoeff # each row sums to 1

  if (d == 0) {
    #vector
    errprbs <- errorprob * segrArray * errcoeff
    segrArray <- (1-errorprob) * segrArray + colSums(errprbs)
  } else {
    #assume 3-dim array with dim 1 and 2 the parental dosages:
    for (p1 in 1:dim(segrArray)[1]) for (p2 in 1:dim(segrArray)[2]) {
      errprbs <- errorprob * segrArray[p1, p2, ] * errcoeff
      segrArray[p1, p2, ] <- (1-errorprob) * segrArray[p1, p2, ] +
        colSums(errprbs)
    }
  }
  segrArray
} #addErrorProb

getOrigPopstruct <- function(samplenames,
                             pop.parents, population, parentalPriors,
                             try.HW, ploidy) {
  #This function checks if the necessary data are available to define
  #subpopulations, and if this information needs to be reformatted
  #samplenames: a sorted character vector with all samplenames occurring in the
  #             polyploid data for any of the markers
  #pop.parents: either a data.frame (user supplied) with 3 columns: the first
  #             with a population ID (name or number), the 2nd and 3rd with the
  #             population IDs of the parents of these populations (if F1's)
  #             or NA (if not),
  #             or a matrix with 2 columns and rownames (already processed),
  #             or NA (no popstruct, all samples belong to the same population)
  #population: a data.frame with two columns, the first containing the
  #            SampleName (containing at least all SampleNames occurring in
  #            polyploid data), the second column containing population IDs
  #            that all match pop.parents;
  #            if pop.parents is a pre-processed matrix the second column of
  #            population should contain integers indexing the pop.parents rows.
  #            In both columns no NA may occur.
  #parentalPriors: a data frame with one column MarkerName followed by one
  #            column for each F1 parent (or optionally two columns if there is
  #            only one F1 population), one row per marker.
  #            Column names (except first) are population IDs matching
  #            the parental populations in pop.parents.
  #      Note: when reading the data frame with read.table, set check.names=FALSE
  #            so column names are not changed
  #try.HW: logical; if TRUE, the ptype for panels will be "p.HW", else "p.free"
  #ploidy: any even ploidy level > 0
  #Result value:
  # either an error message
  # or a list with components:
  #    $orig_population: data frame with in column 1 (SampleNames) the
  #                      samplenames (sorted) and in column 2 (Population)
  #                      integers indexing the rows of $orig_pop.parents;
  #    $orig_pop.parents: matrix with 2 columns, and population IDs as rownames,
  #                       one row for each population, each parent as the index
  #                       to its rwo in pop.parents or twice NA,
  #                       sorted such that each F1 population appears before its
  #                       parents
  #    $segrArray: 3D array with for all dosages of parent1 and parent2 (the
  #                first 2 dims) the expected frequencies of all dosages in an
  #                F1 (note that in all 3 dims the dosage is 0:ploidy but the
  #                array indices are 1:(ploidy+1))
  #    $orig_ptype: character vector with the ptype for all populations
  #                 (including "p.nodata" for populations with no data)
  #    $orig_parents: an integer vector with all rows in pop.parents that are
  #                   parents of F1 populations
  #    $orig_nonparents: an integer vector with all rows in pop.parents that are
  #                      not parents of F1 populations (i.e. F1's or panels)
  #    $orig_F1s: an integer vector with all rows in pop.parents that are
  #               F1 populations
  #    $orig_panels: an integer vector with all rows in pop.parents that are
  #                  panels
  #    $orig_parPriorCols: NULL if no parentalPriors, else a matrix with one row
  #                  per F1 population and one column per parent, containing the
  #                  column in parentalPriors with the priors of that parent.
  #                  Optionally, with only one F1 population (i.e. the matrix
  #                  has only 1 row), there may be 2 columns for each parent.
  #                  The row names are the population IDs of the F1s. All
  #                  elements of the matrix are specified (no NA's)

  # In case population and pop.parents are not specified these components
  # get the values:
  #    $orig_population: data frame with the SampleNames and a column with only
  #                      1's, containing only the samples in dat
  #    $orig_pop.parents: matrix(c(NA, NA)), ncol=2) #1 row, no rownames
  #    $orig_parentalPriors: NULL
  #    other elements are not present in that case.
  #Note that it is an error if population contains values that do not match
  #any row of pop.parents, but that not all rows of pop.parents need to be
  #indexed in populations (i.e. pop.parents may specify additional populations
  #that are not used in population)

  nopopstruct <- function(samplenames) {
    list(orig_population = data.frame(SampleName = samplenames,
                                      Population = rep(1, length(samplenames))),
         orig_pop.parents = matrix(c(NA, NA), ncol=2,
                                   dimnames= list("1", c("P1", "P2"))),
         orig_ptype = ifelse(try.HW, "p.HW", "p.free"),
         orig_parents = integer(0),
         orig_nonparents = 1,
         orig_F1s = integer(0),
         orig_panels = 1,
         parPriorCols = NULL,
         segrArray = NA)
  } # nopopstruct within get.popstruct

  if (is.null(pop.parents) || is1NA(pop.parents)) {
    if ((!is.null(population) && !is1NA(population)) ||
        (!is.null(parentalPriors) && !is1NA(parentalPriors))) {
      return("population and parentalPriors may not be specified without pop.parents")
    } else return(nopopstruct(samplenames))
  }

  #pop.parents not NULL or NA
  if (!(class(pop.parents) %in% c("data.frame", "matrix")))
    return("pop.parents should be a data.frame or matrix")
  if (is.null(population) || is1NA(population))
    return("population not specified while pop.parents is specified")
  if (!is.data.frame(population) ||
      length(population) != 2)
    return("population must be a data frame with 2 columns")
  if (anyNA(population[,1]) || anyNA(population[,2]))
    return("population contains NA values")
  if (length(unique(population[,1])) < nrow(population))
    return("population contains duplicate SampleNames")
  if (is.factor(population[,2])) {
    population[,2] <- as.character(population[,2]) #population IDs
    suppressWarnings( pop2num <- as.numeric(population[,2]) )
    if (!anyNA(pop2num)) population[,2] <- pop2num
  }

  if (is.matrix(pop.parents)) {
    if (ncol(pop.parents) != 2)
      return("matrix pop.parents should have 2 columns")
    if (!all(population[,2] %in% 1:nrow(pop.parents)))
      return("population indexes non-existing rows in matrix pop.parents")
    #population and pop.parents are already matching integer / matrix;
    #therefore pop.parents should be ok without sorting,
    #but we check anyway:
    for (p in nrow(pop.parents)) {
      if(xor(is.na(pop.parents[p, 1]), is.na(pop.parents[p, 2])) ||
         (!is.na(pop.parents[p, 1]) && min(pop.parents[p,]) <= p) ) {
        #last line checks that all parents appear below their F1
        return("invalid matrix pop.parents")
      }
    }
    if (is.null(rownames(pop.parents)))
      rownames(pop.parents) <- 1:nrow(pop.parents)
    popstruct <- list(
      orig_population=data.frame(
        SampleName = samplenames,
        Population = population[match(samplenames, population[,1]), 2]),
      orig_pop.parents=pop.parents)
    if (anyNA(popstruct$orig_population[, 2]))
      return("not all samples specified in population")
  } else {
    #pop.parents is a data.frame
    if (length(pop.parents) != 3 )
      return("data.frame pop.parents should have 3 columns")
    if (anyNA(pop.parents[,1]))
      return("NAs among population ID's in data.frame pop.parents")
    for (i in 1:3) {
      if (is.factor(pop.parents[,i]))
        pop.parents[,i] <- as.character(pop.parents[,i])
    }
    if (length(unique(pop.parents[,1])) != nrow(pop.parents))
      return("duplicated population IDs in data.frame pop.parents")
    for (i in 2:3) {
      emptylines <- !is.na(pop.parents[,i]) & pop.parents[,i] == ""
      pop.parents[emptylines, i] <- NA
      nonna <- pop.parents[!is.na(pop.parents[, i]), i]
      if (anyNA(match(nonna, pop.parents[,1])))
      return("not all parental population IDs occur in data.frame pop.parents")
    }
    if (any(xor(is.na(pop.parents[,2]), is.na(pop.parents[,3]))))
      return("all populations in data.frame pop.parents must have 0 or 2 parents")
    pop.parents <- sortPedigree(pop.parents, 1, 2, 3, parentsFirst=FALSE)
    if (is.character(pop.parents)) {
      #an error message was returned
      return(paste("invalid pop.parents:", pop.parents))
    }
    #convert population and pop.parents to integer / matrix:
    poppar <- matrix(integer(2 * nrow(pop.parents)), ncol=2)
    rownames(poppar) <- pop.parents[,1]
    for (r in 1: nrow(pop.parents)) {
      if (is.na(pop.parents[r, 2])) {
        poppar[r,] <- c(NA, NA)
      } else {
        poppar[r, 1] <- which(pop.parents[,1] == as.character(pop.parents[r, 2]))
        poppar[r, 2] <- which(pop.parents[,1] == as.character(pop.parents[r, 3]))
      }
    }
    population[,2] <- match(population[,2], rownames(poppar))
    if (anyNA(population[,2]))
      return("not all population IDs in population occur in pop.parents")
    popstruct <- list(
      orig_population=data.frame(
        SampleName = samplenames,
        Population = population[match(samplenames, population[, 1]), 2]),
      orig_pop.parents=poppar)
    if (anyNA(popstruct$orig_population[, 2]))
      return("not all samples specified in population")
  }
  pop.parents <- popstruct$orig_pop.parents
  population <- popstruct$orig_population

  # add ptype:
  ptype <- rep(ifelse(try.HW, "p.HW", "p.free"),
               nrow(pop.parents)) #all pop are panels unless changed later
  ptype[unique(pop.parents)] <- "p.free" #parents of F1 populations
  for (r in 1: nrow(pop.parents)) {
    if (!is.na(pop.parents[r,1])) {
      #F1 population
      if (!any(population == pop.parents[r, 1]) ||
          !any(population == pop.parents[r, 2])) {
        return(paste("missing parental samples for F1 population",
                     rownames(pop.parents)[r]))
      }
      ptype[r] <- "p.F1"
    }
    if (!any(population == r)) ptype[r] <- "p.nodata"
  }
  popstruct$orig_ptype <- ptype

  #add some info about the original population structure:
  origpopcount <- nrow(popstruct$orig_pop.parents)
  popstruct$orig_parents <- sort(unique(as.integer(popstruct$orig_pop.parents)))
  popstruct$orig_parents <- popstruct$orig_parents[!is.na(popstruct$orig_parents)]
  popstruct$orig_nonparents <- setdiff(1:origpopcount, popstruct$orig_parents)
  popstruct$orig_F1s <- which(!is.na(popstruct$orig_pop.parents[,1]))
  popstruct$orig_panels <- setdiff(popstruct$orig_nonparents, popstruct$orig_F1s)
  #Note that some of the panels and F1s, but not the parents, may have "p.nodata";
  #ALL of the parents in pop.parents are sure to have some samples at this point

  if ("p.F1" %in% popstruct$orig_ptype) {
    if (ploidy %% 2 != 0) {
      return("F1 populations not allowed with odd ploidy")
    } else popstruct$segrArray <- addErrorProb(getPolysomicSegr(ploidy))
  } else popstruct$segrArray=NA

  #next we check the parentalPriors:
  if (is.null(parentalPriors) || is1NA(parentalPriors)) {
    #add a NULL orig_parPriorCols element to list:
    #popstruct$orig_parPriorCols <- NULL (does not work, deletes the element)
    popstruct <- c(popstruct, list(orig_parPriorCols=NULL))
  } else {
    if (length(popstruct$orig_F1s) == 0)
      return("parentalPriors may not be specified without F1 population")
    if (!is.data.frame(parentalPriors))
      return("parentalPriors should be a data.frame")
    if (names(parentalPriors)[1] != "MarkerName")
      return("name of 1st column of parentalPriors must be MarkerName")
    if (anyNA(parentalPriors$MarkerName) ||
        anyDuplicated(parentalPriors$MarkerName) > 0)
      return("no missing or duplicated MarkerNames allowed in parentalPriors")
    priorvals <- unique(as.vector(as.matrix(parentalPriors[,-1])))
    priorvals <- priorvals[!is.na(priorvals)]
    if (anyNA(match(priorvals, 0:ploidy)))
      return("parentalPriors contains values not in 0:ploidy")

    parnames <- rownames(pop.parents)[popstruct$orig_parents]
    parpriornames <- unique(names(parentalPriors)[-1])
    if (!setequal(parnames, parpriornames))
      return("names of parentalPriors don't match the F1 parents")
    #Note that with only 1 F1 population there may be 2 columns of priors per parent
    if (anyDuplicated(names(parentalPriors)[-1]) > 0) ncol <- 4 else ncol <- 2
    parPriorCols <- matrix(NA_integer_, ncol=ncol,
                           nrow=length(popstruct$orig_F1s))
    if (ncol == 4) {
      if (length(popstruct$orig_F1s) > 1) {
        return("with multiple F1 populations, parentalPriors may have only one prior per parent")
      } else {
        pop <- popstruct$orig_F1s #one F1
        parent1cols <- which(names(parentalPriors) ==
                               row.names(pop.parents)[pop.parents[pop, 1]])
        parent2cols <- which(names(parentalPriors) ==
                               row.names(pop.parents)[pop.parents[pop, 2]])
        if (length(parent1cols) == 2 && length(parent2cols) == 2) {
          parPriorCols[1,] <- c(parent1cols, parent2cols)
        } else {
          return("with one F1 population, both parents should have 1 or both 2 columns in parentalPriors")
        }
      }
    } else {
      #ncol==2, each parent has 1 column in parentalPriors
      for (i in seq_along(popstruct$orig_F1s)) {
        pop <- popstruct$orig_F1s[i]
        for (par in 1:2) {
          parentcol <- which(names(parentalPriors) ==
                               row.names(pop.parents)[pop.parents[pop, par]])
          if (length(parentcol) == 1) parPriorCols[pop, par] <- parentcol
        }
        if (anyNA(parPriorCols[pop,]))
          return("all F1 parents must have a column in parentalPriors")
        # else we may get lookup errors or need to do additional checking later
      }
    }
    popstruct$orig_parPriorCols <- parPriorCols
  }
  popstruct
} #getOrigPopstruct

getMarkerPopstruct <- function(seldat, origpopstruct, try.HW, ng) {
  #Here we assign an alternative p.type where needed:
  #populations that have no data will get a ptype="p.nodata",
  #and F1 populations where one or both parents have no data will get
  #a ptype "p.HW" or "p.free" instead of "p.F1", changing them to panels for the
  #current marker;
  #also a p.start for each population is calculated
  #seldat: data frame as "data" argument of fitOneMarker with only the data
  #        for the current marker, with the samples with NA ratio already
  #        removed
  #origpopstruct: a list as returned by getOrigPopstruct
  #try.HW: logical; if TRUE, the ptype for panels will be "p.HW", else "p.free"
  #ng: number of genotypes (= ploidy + 1)
  #return value: a list with 3 components:
  # $curr_pop: the orig_population including only the samples in seldat
  # $curr_ptype: the orig_ptype modified if a population has no samples
  #              or no parent(s)
  # $curr_pstart: uniform 1/ng for all pop and all dosages

  pop <- origpopstruct$orig_population[match(seldat$SampleName,
                                         origpopstruct$orig_population[,1]), 2]
  #check - should never happen:
  if (anyNA(pop))
    stop(paste("pop contains NA at marker", seldat$MarkerName[1]))
  #pop now one column, in same order as seldat
  ptype <- origpopstruct$orig_ptype
  for (r in 1:nrow(origpopstruct$orig_pop.parents)) {
    if (ptype[r] != "p.nodata" && !any(pop == r)) ptype[r] <- "p.nodata"
  }
  for (r in origpopstruct$orig_F1s) {
    if (ptype[r] == "p.F1" &&
        (ptype[origpopstruct$orig_pop.parents[r, 1]] == "p.nodata") ||
        (ptype[origpopstruct$orig_pop.parents[r, 2]] == "p.nodata"))
      # one or both parents no data, treat F1 as panel for this marker:
      ptype[r] <- ifelse(try.HW, "p.HW", "p.free")
  }
  mrkpopstruct <- list()
  mrkpopstruct$curr_pop <- pop
  mrkpopstruct$curr_ptype <- ptype
  #default p.start, will be replaced if priors available:
  mrkpopstruct$curr_p.start <-
    matrix(1/ng, ncol=ng, nrow=nrow(origpopstruct$orig_pop.parents))
  mrkpopstruct
} #getMarkerPopstruct

getModelPopstruct <- function(model, origPopstruct, mrkPopstruct, ng) {
  #get population structure for current model to be used for CodomMarker:
  #if (nrow(origPopstruct$orig_pop.parents) == 1 || ((model - 1) %% 8) < 4) {
  if (((model - 1) %% 8) < 4) {
    #we use all populations (could be one as in old times), use popstruct$curr:
    list(population =  mrkPopstruct$curr_pop,
         pop.parents = origPopstruct$orig_pop.parents,
         nonparents =  origPopstruct$orig_nonparents,
         ptype =       mrkPopstruct$curr_ptype,
         p.start =     mrkPopstruct$curr_p.start,
         segrArray =   origPopstruct$segrArray,
         allcombined = FALSE)
  } else {
    #we have multiple populations (possibly one left) but here we combine all
    #samples into one:
    #ng <- dim(origPopstruct$segrArray)[3] #doesn't work as segrArray may not be calculated
    list(population =  rep(1, length(mrkPopstruct$curr_pop)),
         pop.parents = matrix(c(NA, NA), nrow=1),
         nonparents =  1,
         ptype =       "p.free",
         p.start =     matrix(1/ng, ncol=ng, nrow=1),
         segrArray =   origPopstruct$segrArray,
         allcombined = TRUE)
  }
} #getModelPopstruct

find.dips.populations <- function(mat, pops, minflank=0.01, minfrac=0.1) {
  #mat: a matrix of ng fractions per row, one row for each population
  #pops: vector indicating which populations to analyze (the non-parent ones)
  #returns a vector with the dip positions in any of the non-parent populations

  find.dips <- function(x, minflank=0.01, minfrac=0.1) {
    #x is a vector of ng fractions summing to 1, of length >=3, without NA's
    #result: a vector with the positions of the dips in x, or integer(0) if none,
    #where a dip is calculated:
    #take the maximum fraction left and the maximum fraction right of the position;
    #take the smallest of these two flanking peaks, -> mintop
    #we recognize a dip at position i if
    #  mintop > minflank (i.e. at both sides of i there is a sizeable peak), and
    #  x[i] < (1-minfrac) * mintop (i.e. x[i] is substantially smaller
    #                             than the lowest flanking peak)
    lx <- length(x)
    top1 <- rep(NA, lx-2); top2 <- top1
    for (i in 2:(lx-1)) {
      top1[i-1] <- max(x[1:(i-1)])  #max x below i
      top2[i-1] <- max(x[(i+1):lx]) #max x above i
    }
    mintop <- pmin(top1, top2)  #pmin is "parallel minimum", length=lx-2
    depth <- (mintop - x[2:(lx-1)]) / mintop # depth as fraction of mintop; <=0: no dip
    dips <- (mintop > minflank) & (depth > minfrac)
    which(dips) + 1 #add 1 since vector dips corresponds to x[2 : (lx-1)]
  } #find.dips within find.dips.populations

  dips <-integer(0)
  #print(paste("pops=",paste(pops, collapse=" ")," nrow(mat)=",nrow(mat)))
  #print(mat)
  for (p in pops) dips <- union(dips, find.dips(mat[p,], minflank, minfrac))
  sort(dips)
} #find.dips.populations

is1NA <- function(x) {
  class(x)[1] %in% c("logical", "integer", "numeric", "character") &&
    length(x) == 1 && is.na(x)
}

#'@title A function to convert a set of mixture means from one ploidy to another
#'
#'@description convertStartmeans takes a set of means at one ploidy level (e.g.
#'the fitted means for a tetraploid data set)
#'and uses them to generate a set of means for another ploidy level (e.g. as
#'startmeans for fitting triploid data for the same markers).
#'
#'@usage convertStartmeans(ploidy, origmeans)
#'@param ploidy The ploidy to which the means must be converted.
#'@param origmeans A data.frame with a first column MarkerName, followed
#'by <oldploidy+1> columns (names are ignored) that contain the ratio means
#'for dosages 0 to <oldploidy>. Column MarkerName may not contain missing values.
#'On each row the other columns must either all contain NA, or only non-NA
#'values between 0 and 1 in strictly ascending order.
#'@return A data.frame like origmeans with the same column MarkerName, now
#'followed by <ploidy+1> columns with the new means.
#'@details The new means are calculated by linear interpolation between the old
#'means on the asin(sqrt(x)) transformed scale and back-transformed to the
#'original scale; the new means for dosage 0 are equal to the old, and the
#'new means for dosage <ploidy> are equal to the old means for dosage
#'<oldploidy>.
#'@examples
#'# means from tetraploid data set:
#'tetrameans <- data.frame(MarkerName=c("mrk1", "mrk2"), mu0=c(0.02, 0.0),
#'mu1=c(0.2, 0.25), mu2=c(0.3, 0.5), mu3=c(0.4, 0.75), mu4=c(0.6, 1.0))
#'# convert to means for triploid data set:
#'trimeans <- convertStartmeans(ploidy=3, origmeans=tetrameans)
#'tetrameans
#'trimeans
#'
#'@export
convertStartmeans <- function(ploidy, origmeans) {
  #This is an exported user function that allows to obtain e.g. startmeans for
  #a triploid sample set based on fitted models from a tetraploid sample set.
  #obtain the new (ploidy+1) start means by linear interpolation among the
  #(oldploidy+1) original means (interpolation on the asin(sqrt(x) transformed
  #scale,input and output on non-transformed scale)
  #ploidy: target ploidy
  #origmeans: data frame with columns MarkerName and (oldploidy+1) means
  #           (on non-transformed scale)
  #return value: as origmeans but now with (ploidy+1) instead of (oldploidy+1)
  #              means after the MarkerName column
  if (!is.data.frame(origmeans) || length(origmeans) < 3)
    stop("convertStartmeans: invalid origmeans")
  if (names(origmeans)[1] != "MarkerName")
    stop("convertStartmeans: name of first column must be MarkerName")
  if (length(unique(origmeans[, 1])) != nrow(origmeans))
    stop("convertStartmeans: duplicate marker names in origmeans")
  ploidy <- as.integer(ploidy)
  if (length(ploidy) != 1 || ploidy < 2)
    stop("convertStartmeans: invalid ploidy")
  origmeans[, 2:length(origmeans)] <-
    asin(sqrt(origmeans[, 2:length(origmeans)]))
  origploidy <- length(origmeans)-2
  origdose <- (0:origploidy)/origploidy
  result <- matrix(NA_real_, ncol=ploidy+1, nrow=nrow(origmeans))
  result[, 1] <- origmeans[, 2]
  result[, ploidy+1] <- origmeans[, length(origmeans)]
  for (i in 2:(ploidy)) {
    dose <- (i-1)/ploidy - 1e-8
    startinterval <- which.max(origdose > dose) - 1
    result[, i] <- origmeans[, startinterval+1] +
      (dose-origdose[startinterval]) * origploidy *
      (origmeans[, startinterval+2] - origmeans[, startinterval+1])
  }
  sinmeans <- sin(result)
  data.frame(MarkerName=origmeans[, 1],
             sinmeans*sinmeans)
} #convertStartmeans

get.pHW <- function(ploidy, allelefrqB) {
  #allelefrqB: allele frequency of the B (or Y) allele
  #return value: numeric vector of length (ploidy+1) with the HW probabilities
  #of dosages 0:ploidy
  dbinom(0:ploidy, ploidy, allelefrqB)
} #get.pHW

get.ps <- function(origPopstruct, mrkPopstruct, parentalPriors){
  #parentalPriors is (here) a matrix with one row per F1 and two columns with
  #the priors (one for each parent, one or both may be NA)
  #return value: a modified version of mrkPopstruct with:
  # $curr_p.start: a matrix with for each population in origPopstruct$orig_pop.parents
  #                one row with the ploidy+1 p's to be used as start.p
  #                For parents with a prior available: p=0.98 for the prior,
  #                rest for other dosages
  #                For F1's with both parents known: use the polysomic prediction,
  #                but if that is monomorphic, do as for the parents
  #                For F1's with one parental prior known: use HW ratios based
  #                on the allele freq of the parents, assuming the unknown parent
  #                has dosage ploidy/2
  #                For all other cases (parents without prior, F1's with both
  #                parents unknown, panels): use uniform p's (1/(ploidy+1)
  #                for all dosages)
  # $curr_ptype:   For all parents with a known prior (even if the other parent
  #                of an F1 has no known prior), change from p.free to p.fixed,
  #                Also for F1s with 2 parents with known prior change to p.fixed
  ng <- ncol(mrkPopstruct$curr_p.start)
  ploidy <- ng - 1
  for (f1 in seq_along(origPopstruct$orig_F1s)) {
    parpri <- parentalPriors[f1,] #becomes a vector
    F1pop <- origPopstruct$orig_F1s[f1]
    parents <- origPopstruct$orig_pop.parents[F1pop,] #becomes a vector
    #set the p.start and ptype for the parents:
    for (par in 1:2) {
      if (is.na(parpri[par])) {
        if (mrkPopstruct$curr_ptype[parents[par]] != "p.nodata")
          mrkPopstruct$curr_ptype[parents[par]] <- "p.free"
      } else {
        mrkPopstruct$curr_p.start[parents[par],] <- #as.integer(0:ploidy == parpri[par])
          #getPriorP(ploidy, parpri[par])
          addErrorProb(0 + ((0:ploidy) == parpri[par]))
        if (mrkPopstruct$curr_ptype[parents[par]] != "p.nodata")
          mrkPopstruct$curr_ptype[parents[par]] <- "p.fixed"
      }
    }
    #set the p.start for the F1 populations:
    #whether parental data (ratios) are available or not, we can still use
    #the parental priors to get a better p.start:
    if (sum(is.na(parpri)) == 1) {
      #one parent prior known, set p.start for F1 to HW ratios
      #where the allele freq is calculated with unknown prior set to ploidy/2:
      mrkPopstruct$curr_p.start[F1pop,] <-
        get.pHW(ploidy=ploidy,
                allelefrqB=(sum(parpri, na.rm=TRUE) + (ploidy/2)) / (2*ploidy))
    } else if (!anyNA(parpri)) {
      #priors for both parents known
      mrkPopstruct$curr_p.start[F1pop,] <-
        origPopstruct$segrArray[parpri[1]+1, parpri[2]+1,]
      mrkPopstruct$curr_ptype[F1pop] <- "p.fixed"
    }
  } #for f1
  #we don't assign HW probabilities to panels as we have no idea of the
  #allele freq (so leave them as uniform probabilities):
  mrkPopstruct
} # get.ps

calcMus <- function(seldat, priorcomb, priorcomb_index,
                   origpopstruct, mrkpopstruct,
                   samplepriors, startmus, ng) {
  #seldat: data frame with the selected data (no NA ratios)
  #priorcomb: NA, or the list of parental priors matrices returned by
  #           getPriorCombinations: each with one row per F1, 1 column
  #           per parent
  #priorcomb_index: which of these matrices to use, (ignored if priorcomb is NA)
  #origpopstruct: a list as returned by getOrigPopstruct
  #mrkpopstruct: a list as returned by getMarkerPopstruct; the $population
  #              member matches the ratios in seldat (same order)
  #samplepriors: NA, or the row of the original samplePriors data.frame for the
  #              current marker, with the first column (MarkerName) omitted
  #startmus: NA, or a vector of ng mus for the current marker
  #ng: ploidy+1
  #If non-missing startmus are given these are used instead of deriving
  #startmus from priors and ratios.
  #Else first the priors and ratios from priorcomb and samplepriors are
  #combined into one data frame (where the samplepriors overwrite the priorcomb)
  #and this is used to calculate means for the dosages.
  #If there is only one prior dosage it is still used to fix one of the mus
  #If the priors are too strange or there are no priors, NA is returned

  if (!is.null(startmus) && !is1NA(startmus)) return(startmus)
  #assemble the parentalPriors and samplePriors into two matching vectors:
  dosages <- ratios <- numeric(0)
  #first parentalPriors:
  if (length(priorcomb) > 0) {
    prico <- priorcomb[[priorcomb_index]]
    for (f1 in seq_along(origpopstruct$orig_F1s)) {
      F1pop <- origpopstruct$orig_F1s[f1]
      for (par in 1:2) {
        dos <- prico[f1, par] #prior dosage of parent par of F1-population f1
        if (!is.na(dos)) {
          parpop <- origpopstruct$orig_pop.parents[F1pop, par]
          parrat <- seldat$ratio[mrkpopstruct$curr_pop == parpop]
          dosages <- c(dosages, rep(dos, length(parrat)))
          ratios <- c(ratios, parrat)
        }
      }
    }
  }
  #next samplepriors:
  if (!is.null(samplepriors) && !is1NA(samplepriors) &&
      (!is.list(priorcomb) || length(priorcomb) == 1)) {
    #we use (also) the samplepriors (we can be sure that
    #even if there are 2 priors per parent in the data set, only one (max) of
    #each is not NA)
    for (i in 1:length(samplepriors)) {
      samprat <- seldat$ratio[seldat$SampleName == names(samplepriors)[i]]
      if (length(samprat) == 1) {
        ratios <- c(ratios, samprat)
        dosages <- c(dosages, samplepriors[1, i])
      }
    }
  }
  noNA <- !is.na(ratios) & !is.na(dosages)
  if (sum(noNA) == 0) return(NA)

  #calculate startmus from priors and dosages:
  dosages <- dosages[noNA]
  ratios <- asin(sqrt(ratios[noNA]))
  pi2002 <- pi/2 - 0.02
  skewlimit <- 0.3 #if all mus below skewlimit or above pi/2-skewlimit return NA
  uniqdos <- sort(unique(dosages))
  if (length(uniqdos) > 1) {
    #based on regression, take out-of-range estimates from initmus
    #Alternative: calculate mean for each dosage, interpolate dosages
    #that have no priors. Advantage: better fit for non-linear response. But
    #we don't do that as some means may be represented by more samples,
    #regression takes weights into account.
    coefs <- tryCatch(lm(ratios~dosages)$coefficients,
                      error=function(e) c(-1, -1))
    if (coefs[2] <= 0) return(NA) #negative slope
    mus <- coefs[1] + (0:(ng-1) * coefs[2])
    if (mus[1] > pi/2-skewlimit || mus[ng] < skewlimit) return(NA)
    lo <- mus < 0.02 #first mu should be >= 0.02
    if (any(lo))
      mus[lo] <- seq(0.02, mus[!lo][1], length.out=sum(lo)+1)[1:sum(lo)]
    hi <- mus > pi2002 #last mu should be <= p2002
    if (any(hi))
      mus[hi] <- seq(mus[ng-sum(hi)], pi2002, length.out=sum(hi)+1)[-1]
    return(sin(mus)*sin(mus))
  } else {
    #length(uniqdos) == 1)
    mus <- numeric(ng)
    mus[uniqdos+1] <- mean(ratios)
    if ((mus[uniqdos+1] <= 0.02) ||
        (mus[uniqdos+1] >= pi2002) ||
        (uniqdos == 0 && mus[1] > pi/2 - skewlimit) ||
        (uniqdos == ng-1 && mus[ng] < skewlimit))
      return(NA)
    #we distribute the mus to each side equally.
    if (uniqdos > 0)
      mus[1:uniqdos] <-
        seq(0.02, mus[uniqdos+1], length.out=uniqdos+1)[1:uniqdos]
    if (uniqdos < (ng - 1))
      mus[(uniqdos+2):ng] <-
        seq(mus[uniqdos+1], pi2002, length.out = ng - uniqdos)[-1]
    return(sin(mus)*sin(mus))
  }
} #calcMus

getPriorCombinations <- function(priorCols, parPriors) {
  #priorCols: matrix parPriorCols from popstruct: for each F1 one row,
  #           2 or 4 columns (1 or 2 per parent)
  #parPriors: the row of parentalPriors for the current marker
  #return value: a list with matrices, each matrix has one row per (original)
  #F1 and one column for each parents. Contents are the prior
  #dosages of each parent.
  #If there are no original F1's, there are no priorCols and this function
  #is not called.
  #If there are no parPriors (NULL, or data frame with 0 rows) an empty list
  #is returned.
  #If priorCols has only one column per parent, the results is either an empty
  #list (if all priors for all F1 parents are NA) or a list with one matrix.
  #Note that on each row of the matrix 0, 1 or both priors can be NA.
  #If priorCols has two columns per parent, there can be only one F1 so
  #priorCols has only one row. The results is an empty list (if both priors for
  #both parents are NA) or a list with 1 to 4 one-row matrices. Only the
  #combinations where not both parents have an NA prior are generated.
  #TODO: we now allow rows with one missing prior, check if it is needed to have
  #      only both priors present or both absent
  result <- list()
  if (is.null(parPriors) || nrow(parPriors) != 1) return(result)
  if (ncol(priorCols) == 2) {
    #one prior per parent, multiple F1s(rows of priorCols) possible:
    result[[1]] <- matrix(NA_integer_, nrow=nrow(priorCols), ncol=2)
    for (f1 in 1:nrow(priorCols)) for (p in 1:2) {
      result[[1]][f1, p] <- parPriors[, priorCols[f1, p]]
    }
    #result[[1]][rowSums(is.na(result[[1]])) > 0,] <- c(NA, NA) #only if we don't allow one known prior
    if (all(is.na(result[[1]]))) result <- list()
    return(result)
  }
  # one F1, each parent 2 columns of priors (which might be NA for this marker)
  # if one parent has one or two priors and the other none, give these;
  # else give only the combinations of known priors
  parpri <- list()
  for (parent in 1:2) {
    if (parent == 1) parcols <- priorCols[,1:2] else parcols <- priorCols[,3:4]
    parpri[[parent]] <- unlist(parPriors[,parcols])
    parpri[[parent]] <- parpri[[parent]][!is.na(parpri[[parent]])]
  }
  if (length(parpri[[1]]) == 0) {
    #if we only allow both parents non-missing, omit this
    for (i in seq_along(parpri[[2]])) {
      result[[i]] <- matrix(c(NA, parpri[[2]][i]), nrow=1)
    }
  } else {
    if (length(parpri[[2]]) == 0) {
      #if we only allow both parents non-missing, omit this
      for (i in seq_along(parpri[[1]])) {
        result[[i]] <- matrix(c(parpri[[1]][i], NA), nrow=1)
      }
    } else {
      #both parents have non-missing priors
      i <- 0
      for (p1 in seq_along(parpri[[1]])) for (p2 in seq_along(parpri[[2]])) {
        i <- i+1
        result[[i]] <- matrix(c(parpri[[1]][p1], parpri[[2]][p2]), nrow=1)
      }
    }
  }
  result #can be empty list but not NA
} #getPriorCombinations


#'@title Function to fit multiple mixture models to signal ratios of a single
#'bi-allelic marker
#'
#'@description This function takes a data frame with allele signal ratios for
#'multiple bi-allelic markers and samples, and fits multiple mixture models to
#'a selected marker. It returns a list, reporting on the performance of these
#'models, selecting the best one based on the BIC criterion, optionally
#'plotting results.


#'
#'@usage fitOneMarker(ploidy, marker, data, diplo=NULL, select=TRUE,
#'diploselect=TRUE, pop.parents=NULL, population=NULL, parentalPriors=NULL,
#'samplePriors=NULL, startmeans=NULL, maxiter=40, maxn.bin=200, nbin=200,
#'sd.threshold=0.1, p.threshold=0.9, call.threshold=0.6, peak.threshold=0.85,
#'try.HW=TRUE, dip.filter=1, sd.target=NA,
#'plot="none", plot.type="png", plot.dir, sMMinfo=NULL)
#'
#'@param ploidy The ploidy level, 2 or higher: 2 for diploids, 3 for triploids
#'etc.
#'@param marker A marker name of number. Used to select the data for one marker,
#'referring to the MarkerName column of parameter data. If a number, the number
#'of the marker based on alphabetic order of the MarkerNames in data.
#'@param data A data frame with the polyploid samples, with (at least) columns
#'MarkerName, SampleName and ratio, where ratio is the Y-allele signal
#'divided by the sum of the X- and Y-allele signals: ratio == Y/(X+Y)
#'@param diplo NULL or a data frame like data, with the diploid samples and (a subset
#'of) the same markers as in data. Genotypic scores for diploid samples are
#'calculated according to the best-fitting model calculated for the polyploid
#'samples and therefore may range from 0 (nulliplex) to <ploidy>, with the
#'expected dosages 0 and <ploidy> for the homozygotes and <ploidy/2> for the
#'heterozygotes.\cr
#'Note that diplo can also be used for any other samples that need to be
#'scored, but that should not affect the fitted models.
#'@param select A logical vector, recycled if shorter than nrow(data):
#'indicates which rows of data are to be used (default TRUE, i.e. keep all rows)
#'@param diploselect A logical vector like select, matching diplo instead of data
#'@param pop.parents NULL or a data.frame specifying the population structure. The
#'data frame has 3 columns: the first containing population IDs, the 2nd and 3rd
#'with the population IDs of the parents of these populations (if F1's) or NA
#'(if not). The poopulation IDs should match those in parameter population. If
#'pop.parents is NULL all samples are considered to be in one population, and
#'parameter population should also be NULL (default).
# Alternative, in case fitOneMarker is called from fitMarkers:
# a matrix with 2 columns and 1 row per population; the cells contain the row
# numbers of the parental populations in case of an F1 and NA otherwise. The
# rownames are the population ISd, and the rows must be sorted such that all
# F1s occur above their parental populations.
#'@param population NULL or a data.frame specifying to which population each
#'sample belongs. The data frame has two columns, the first containing
#'the SampleName (containing all SampleNames occurring in data),
#'the second column containing population IDs that match pop.parents. In both
#'columns NA values are not allowed. Parameters pop.parents and population
#'should both be NULL (default) or both be specified.
# Alternative, in case fitOneMarker is called from fitMarkers and pop.parents
# is a matrix: then the second column of population should contain integers
# indexing the pop.parents rows.
#'@param parentalPriors NULL or a data frame specifying the prior dosages for
#'the parental populations. The data frame has one column MarkerName
#'followed by one column for each F1 parental population. Column names (except
#'first) are population IDs matching the parental populations in pop.parents.
#'In case there is just one F1 population in pop.parents, it is possible to
#'have two columns for both parental populations instead of one (allowing two
#'specify two different prior dosages); in that case both columns for each
#'parent have the same caption. Each row specifies the priors for
#'one marker. The contents of the data frame are dosages, as integers from 0
#'to <ploidy>; NA values are allowed.\cr
#'Note: when reading the data frame with read.table or read.csv, set
#'check.names=FALSE so column names (population IDs) are not changed.
#'@param samplePriors NULL or a data.frame specifying prior dosages for individual
#'samples. The first column called MarkerName is followed by one column per
#'sample; not all samples in data need to have a column here, only
#'those samples for which prior dosages for one or more markers are available.
#'Each row specifies the priors for one marker. The contents of the data frame
#'are dosages, as integers from 0 to <ploidy>; NA values are allowed.\cr
#'Note: when reading the data frame with read.table or read.csv, set
#'check.names=FALSE so column names (population IDs) are not changed.
#'@param startmeans NULL or a data.frame specifying the prior means of
#'the mixture distributions. The data frame has one column MarkerName,
#'followed by <ploidy+1> columns with the prior means on the original
#'(untransformed) scale. Each row specifies the
#'means for one marker in strictly ascending order (all means NA is allowed, but
#'markers without start means can also be omitted).
#'@param maxiter A single integer, passed to CodomMarker, see there for explanation
#'@param maxn.bin A single integer, passed to CodomMarker, see there for explanation
#'@param nbin A single integer, passed to CodomMarker, see there for explanation
#'@param sd.threshold The maximum value allowed for the (constant) standard
#'deviation of each peak  on the arcsine - square root transformed scale,
#'default 0.1. If the optimal model has a larger standard deviation the marker
#'is rejected. Set to a large value (e.g. 1) to disable this filter.
#'@param p.threshold The minimum P-value required to assign a genotype (dosage)
#'to a sample; default 0.99. If the P-value for all possible genotypes is less
#'than p.threshold the sample is assigned genotype NA. Set to 1 to disable
#'this filter.
#'@param call.threshold The minimum fraction of samples to have genotypes
#'assigned ("called"); default 0.6. If under the optimal model the fraction of
#'"called" samples is less than call.threshold the marker is rejected. Set to 0
#'to disable this filter.
#'@param peak.threshold The maximum allowed fraction of the scored samples that
#'are in one peak; default 0.85. If any of the possible genotypes (peaks in the
#'ratio histogram) contains more than peak.threshold of the samples the marker
#'is rejected (because the remaining samples offers too little information for
#'reliable model fitting).
#'@param try.HW Logical: if TRUE (default), try models with and without a
#'constraint on the mixing proportions according to Hardy-Weinberg equilibrium
#'ratios. If FALSE, only try models without this constraint. Even when the HW
#'assumption is not applicable, setting try.HW to TRUE often still leads to
#'a better model. For more details on how try.HW is used see the Details
#'section.
#'@param dip.filter if 1 (default), select best model only from models
#'that do not have a dip (a lower peak surrounded by higher peaks: these are not
#'expected under Hardy-Weinberg equilibrium or in cross progenies). If all
#'fitted models have a dip still the best of these is selected. If 2, similar,
#'but if all fitted models have a dip the marker is rejected. If 0, select best
#'model among all fitted models, including those with a dip.
#'@param sd.target If the fitted standard deviation of the peaks on the
#'transformed scale is larger than sd.target a penalty is given (see Details);
#'default NA i.e. no penalty is given.
#'@param plot String, "none" (default), "fitted" or "all". If "fitted" a plot
#'of the best fitting model and the assigned genotypes is saved with filename
#'<marker number><marker name>.<plot.type>, preceded by "rejected_" if the
#'marker was rejected. If "all", small plots of all models are saved to files
#'(8 per file) with filename
#'<"plots"><marker number><A..F><marker name>.<plot.type> in addition to the
#'plot of the best fitting model.
#'@param plot.type String, "png" (default), "emf", "svg" or "pdf". Indicates
#'format for saving the plots.
#'@param plot.dir String, the directory where to save the plot files. Must be
#'specified if plot is not "none". Set this to "" to save plot files
#'in the current working directory.
#'@param sMMinfo NULL (default), for internal use only. Prevents unneeded checking
#'and recalculation of input parameters when called from fitMarkers.
#'@details fitOneMarker fits a series of mixture models for the given marker by
#'repeatedly calling CodomMarker and selects the optimal one. The initial
#'models vary according to the values of try.HW, pop.parents,
#'parentalPriors, samplePriors and startmeans:
#'\itemize{
#'  \item no pop.parents, try.HW FALSE: 4 models with different constraints
#'   on the means (different or equal X and Y background signal, ratio a linear or
#'   quadratic function of dosage), no restrictions on the mixing proportions
#'   (the fractions of samples in each dosage peak)
#'  \item no pop.parents, try.HW TRUE: The previous 4 models are fitted and
#'  also 4 models with the same restrictions on the means and the mixing
#'  proportions restricted to Hardy-Weinberg ratios (assuming polysomic
#'  inheritance)
#'  \item pop.parents specified, no parentalPriors / samplePriors / startmeans,
#'  try.HW FALSE: 4 models
#'  are fitted with the same restrictions on the means as above, but with
#'  different restrictions on the mixing proportions for each population:
#'  no restriction on parental populations, none on accession panels, polysomic
#'  F1 segregation ratios on F1 populations. Additionally 4 models are fitted
#'  with all samples considered as one population, with the same 4 models for
#'  the means and no restrictions on mixing proportions.
#'  \item pop.parents specified, no parentalPriors / samplePriors / startmeans,
#'  try.HW TRUE: 4 models
#'  are fitted with the same restrictions on the means as above, but with
#'  different restrictions on the mixing proportions for each population:
#'  no restriction on parental populations, HW-ratios for accession panels,
#'  polysomic F1 segregation ratios on F1 populations. Additionally 4 models
#'  are fitted with all samples considered as one population, with the same 4
#'  models for the means and mixing proportions according to HW ratios.
#'  \item pop.parents and parentalPriors specified, try.HW FALSE: 4 models
#'  are fitted with the same restrictions on the means as above, but with
#'  different restrictions on the mixing proportions for each population:
#'  no restriction on parental populations and the accession panels,
#'  polysomic F1 segregation ratios on F1 populations ignoring the parental
#'  priors. Additionally 4 models are fitted with the same restrictions on the
#'  means and mixing proportions of the accession panels, but where the
#'  mixing proportions of the parental populations are set to (almost) 1 for
#'  the prior dosage and (almost) 0 for all other dosages, and those for the
#'  F1 populations to the polysomic segregation ratios expected for the
#'  parental priors.
#'  \item pop.parents and parentalPriors specified, try.HW TRUE: same as with
#'  try.HW FALSE, except that the mixing proportions of accession panels are now
#'  restricted to HW ratios.
#'  \item if parentalPriors and/or samplePriors are specified, these and the
#'  signal ratios of the corresponding samples are (also) used to estimate starting
#'  values of the mixture component means in the EM algorithm. Alternatively
#'  startmeans can be specified directly.
#'}
#'Because convergence to the optimal solution often fails, the models are fitted
#'with several start values for the <ploidy+1> means of the mixture
#'distributions: (1) based on initial clustering of the ratios, (2) based on
#'a uniform distribition from 0.02 to pi/2-0.02 on the asin(sqrt(x)) scale,
#'and (3) if startmeans are specified or can be calculated from samplePriors
#'and/or parentalPriors these are used for a third set of model fits.\cr
#'The main difference between parentalPriors and samplePriors is that
#'parentalPriors are treated as fixed (and if both parents of an F1 population
#'have priors, the F1 segregation is also fixed) while samplePriors are only
#'used to calculate starting ratio means for each dosage. Depending on the
#'confidence the user has in the prior dosages of the parents they can be
#'supplied as parentalPriors or samplePriors.
#'In some cases an additional fit is performed with a modified set of initial
#'means.\cr
#'An optimal model is selected based on the Bayesian Information Criterium
#'(BIC), which takes into account the Log-Likelihood and the number of fitted
#'parameters of the models. If sd.target is specified and the standard
#'deviation of the mixture model components is larger than this target a
#'penalty is applied, making is less likely that that model is selected.\cr
#'The plots consist of one histogram per (non-parent) population showing the
#'frequency distribution of the signal ratios of the samples in that population.
#'The fitted model is shown in green (density and means), and for F1 populations
#'the samples of parent 1 and 2 are shown as red and blue triangles.\cr
#'If diplodata are present, a histogram for the diploid samples is plotted
#'in the top histogram (diploid bars are narrower and gray). The diploid bars
#'are scaled so the maximum bar is half the maximum polyploid bar. At the
#'bottom of the plot for the fitted model a rug plot shows the scores of each
#'sample, while the bottom (red) samples are unscored.
#'
#'@return A list with components:
#'\describe{
#'  \item{log}{A character vector with the lines of the log text.}
#'  \item{modeldata}{A data frame as allmodeldata (see below) with only the
#'  one row with data on the selected model.}
#'  \item{allmodeldata}{A data frame with for each tried model one row with
#'  the marker number, marker name, number of samples and (if the marker is
#'  not rejected) data of the fitted model (see below).}
#'  \item{scores}{A data frame with the name and data for all samples
#'  (including NA's for the samples that were not selected, see parameter
#'  select), with columns:\cr
#'  marker (the sequential number of the marker (based on alphabetic
#'  order of the marker names in data)\cr
#'  MarkerName\cr
#'  SampleName\cr
#'  ratio (the given ratio from parameter data)\cr
#'  P0 .. P<ploidy> (the probabilities that this sample belongs to each of the
#'  <ploidy+1> mixture components)\cr
#'  maxgeno (0..ploidy, the genotype = mixture component with the highest P
#'  value)\cr
#'  maxP (the P value for this genotype)\cr
#'  geno (the assigned genotype number: same as maxgeno, or NA if
#'  maxP < p.threshold).}
#'  \item{diploscores}{A data frame like scores for the samples in the data
#'  frame supplied with argument diplo. If diplo is NA also diploscores will be
#'  NA.}
#'}
#'The modeldata and allmodeldata data frames present data on a fitted model.
#'modeldata presents data on the selected model; allmodeldata lists all
#'attempted models. Both data frames contain the following columns:
#'\describe{
#'  \item{marker}{the sequential number of the marker (based on alphabetic
#'  order of the marker names in data)}
#'  \item{MarkerName}{the name of the marker}
#'  \item{m}{the number of the fitted model}
#'  \item{model}{the type of the fitted model. Possible values are "b1", "b2", "b1,q",
#'  "b2,q", each by itself or followed by "HW" or "pop". The first 4 refer to
#'  the models for the mixture means: b1 and b2 indicate 1 or 2
#'  parameters for signal background, q indicates that a quadratic term in the
#'  signal response was fitted as well. HW and pop refer to the restrictions on
#'  the mixing proportions: HW indicates that the mixing proportions were
#'  constrained according to Hardy-Weinberg equilibrium ratios in case of only
#'  one population, pop indicates that multiple populations were fitted (see
#'  Details section). For more details see Voorrips et al (2011),
#'  doi:10.1186/1471-2105-12-172.}
#'  \item{nsamp}{the number of samples for this marker for which select==TRUE,
#'  i.e. the number on which the call rate is based.}
#'  \item{nsel}{the number of these samples that have a non-NA ratio value}
#'  \item{npar}{the number of free parameters fitted}
#'  \item{iter}{the number of iterations to reach convergence}
#'  \item{dip}{whether the model had a dip (a smaller peak surrounded by
#'  larger peaks): 0=no, 1=yes}
#'  \item{LL}{the log-likelihood of the model}
#'  \item{AIC}{Akaike's Information Criterion}
#'  \item{BIC}{Bayesian Information Criterion}
#'  \item{selcrit}{the selection criterion; the model with the lowest selcrit
#'  is selected. If argument sd.target is NA selcrit is equal to BIC, else
#'  selcrit is larger than BIC if the standard deviation of the mixture
#'  components is larger than sd.target; see Details for details.}
#'  \item{minsepar}{a measure of the minimum peak separation. Each difference
#'  of the means of two successive mixture components is divided by the average
#'  of the standard deviations of the two components. The minimum of the
#'  values is reported. All calculations are on the arcsine-square root
#'  transformed scale.}
#'  \item{meanP}{For each sample the maximum probability of belonging to any
#'  mixture component is calculated. The average of these P values is reported
#'  in meanP}
#'  \item{P80 .. P99}{the fraction of samples that have a probability of at
#'  least 0.80 .. 0.99 to belong to one of the mixture components (by default a
#'  level of 0.99 is required to assign a genotype score to a sample)}
#'  \item{muact0 ..}{the actual means of the samples in each of
#'  the  mixture components for dosages 0 .. <ploidy> on the transformed scale}
#'  \item{sdact0 ..}{the actual standard deviations of the samples in each of
#'  the  mixture components for dosages 0 .. <ploidy> on the transformed scale}
#'  \item{mutrans0 ..}{the means of the mixture components for dosages 0 ..
#'  <ploidy> on the transformed scale}
#'  \item{sdtrans0 ..}{the standard deviations of the mixture components for
#'  dosages 0 .. <ploidy> on the transformed scale}
#'  \item{P0 ..}{the mixing proportions of the mixture components for dosages
#'  0 to <ploidy>. If multiple populations are specified there are two
#'  possibilities: (1) the specified population structure is used in the
#'  current model; then for each population the mixing proportions are given
#'  as <npop> sequences of <ploidy+1> fractions, or (2) the population
#'  structure is ignored for the current model, the mixing proportions are given
#'  in the first sequence of <ploidy+1> fractions and all following sequences
#'  are filled with NA. The the item names are adapted to have the
#'  population names between the P and the dosage}
#'  \item{mu0 ..}{the model means of the <ploidy+1> mixture components
#'  back-transformed to the original scale}
#'  \item{sd0 ..}{the model standard deviations of the <ploidy+1> mixture
#'  components back-transformed to the original scale}
#'  \item{message}{if no model was fitted or the model was rejected, the reason
#'  is reported here}
#'}
#'@examples
#' \donttest{
#'  # These examples run for a total of about 9 sec.
#'
#'  data(fitPoly_data)
#'
#'  # triploid, no specified populations
#'  fp <- fitOneMarker(ploidy=3, marker="mrk039",
#'                     data=fitPoly_data$ploidy3$dat3x)
#'
#'  # tetraploid, specified populations
#'  # plot of the fitted model saved in tempdir()
#'  fp <- fitOneMarker(ploidy=4, marker=2,
#'                     data=fitPoly_data$ploidy4$dat4x,
#'                     population=fitPoly_data$ploidy4$pop4x,
#'                     pop.parents=fitPoly_data$ploidy4$pop.par4x,
#'                     plot="fitted",
#'                     plot.dir=paste0(tempdir(),"/fpPlots4x"))
#'
#'  # hexaploid, specified populations, start values for means,
#'  # plot of the fitted model saved in tempdir()
#'  fp <- fitOneMarker(ploidy=6, marker=1,
#'                     data=fitPoly_data$ploidy6$dat6x,
#'                     population=fitPoly_data$ploidy6$pop6x,
#'                     pop.parents=fitPoly_data$ploidy6$pop.par6x,
#'                     startmeans=fitPoly_data$ploidy6$startmeans6x,
#'                     plot="fitted", plot.dir=paste0(tempdir(),"/fpPlots6x"))
#' }
#'
#'@export
fitOneMarker <- function (
  ploidy,
  marker, #one integer (index to markernames) or one character string
  #        or, if called from fitMarkers, a list with number and name
  #        of one marker in data
  data, diplo=NULL,
  select=TRUE, diploselect=TRUE, # by default all samples of polyploids and diploids selected
  pop.parents=NULL,
  population=NULL,
  parentalPriors=NULL,
  samplePriors=NULL,
  startmeans=NULL,
  maxiter=40, maxn.bin=200, nbin=200,
  sd.threshold=0.1, p.threshold=0.9,
  call.threshold=0.6, peak.threshold=0.85,
  try.HW=TRUE, dip.filter=1,
  sd.target=NA, #not NULL
  plot="none", plot.type="png", plot.dir,
  sMMinfo=NULL) {

  # version with populations:
  # try.HW has a different effect depending on whether there is a
  # population / pop.parents combination with more than one population.
  # If there is not, the effect of try.HW remains unchanged:
  # the first 4 of each set of 8 models (the ones with p.HW) are skipped if
  # try.HW is FALSE
  # If the population / pop.parents combination has more than one
  # population: in each set of 8 models the first four take the populations
  # into account, and here the "association panels" are analyzed with p.HW if
  # try.HW is TRUE, else with p.free. The second four models do not use the
  # populations, and use p.free for the total set of samples irrespective of
  # the try.HW value.

  mrkresult <- NULL
  rejected <- FALSE
  nselsamp <- NA
  seldat <- NULL
  diplonames <- NULL

  if (is.null(sMMinfo)) {
    #not called from fitMarkers, check input

    if (!(class(ploidy) %in% c("numeric", "integer")) ||
        length(ploidy) != 1 ||
        ploidy < 2) {
      stop("ploidy must be  >= 2")
    }

    # check (diplo)data and select
    # and get all markernames, samplenames and diplonames:
    # marker numbers refer to markers actually present in data,
    # in the alphabetic order according to the current OS and R version
    checkInputData(data)
    select <- checkSelect(select, data)
    markernames <- getSortedLevels(data$MarkerName)
    samplenames <- getSortedLevels(data$SampleName)
    IntendedSamples <- getIntendedSamples(select, data)
    if (!is.null(diplo) && !is1NA(diplo)) {
      checkInputData(diplo)
      diploselect <- checkSelect(diploselect, diplo)
      diplonames <- getSortedLevels(diplo$SampleName)
      #IntendedDiplosamples <- getIntendedSamples(diploselect, diplo)
    }

    # check population structures:
    origpopstruct <- getOrigPopstruct(samplenames=samplenames,
                                      pop.parents=pop.parents,
                                      population=population,
                                      parentalPriors=parentalPriors,
                                      try.HW=try.HW,
                                      ploidy=ploidy)
    if (is.character(origpopstruct))
      stop(paste("invalid population structure:", origpopstruct))

    if (!is.null(startmeans) && !is1NA(startmeans)) checkStartmeans(ploidy, startmeans)
    if (!is.null(samplePriors) && !is1NA(samplePriors))
      checkSamplePriors(ploidy, samplePriors, origpopstruct)

    if (!(dip.filter %in% 0:2))
      stop("dip.filter must be in 0:2")

    # arrange plotting output:
    if (missing(plot.dir)) plot.dir <- NULL
    plot <- checkPlot(plot, plot.type, plot.dir, "")
    if (plot[3] != "") message(plot[3])
    plot.type <- plot[[2]]
    plot.dir <- plot[[4]]
    plot <- plot[[1]]

    # check current marker:
    if (length(marker) != 1 || is.na(marker)) stop("invalid marker")
    if (is.character(marker)) {
      mrknr <- match(marker, markernames)
      if (length(mrknr) == 1) {
        markername <- marker
        marker <- mrknr
      } else stop(paste(marker, "does not occur in data"))
    } else if (marker %in% 1:length(markernames)) {
      markername <- markernames[marker]
    } else stop("invalid marker")

    # end if (!is.list(sMMinfo))
  } else {
    # called from fitMarkers
    markernames <- sMMinfo$markernames
    markername <- markernames[marker]
    samplenames <- sMMinfo$samplenames
    IntendedSamples <- sMMinfo$IntendedSamples
    if ("diplonames" %in% names(sMMinfo)) {
      diplonames <- sMMinfo$diplonames
      #IntendedDiplosamples <- sMMinfo$IntendedDiplosamples
    }
    origpopstruct <- sMMinfo$origpopstruct
  }

  ng <- ploidy + 1 # ng = number of genotypic classes
  plotall <- plot == "all"
  plotfitted <- plotall || plot == "fitted"
  optmodel1 <- NA
  optmodel2 <- NA
  mrknr <- padded(marker, length(markernames)) #adds leading zeroes if needed
  # start log
  log <- ""
  log <- append(log, paste(marker, "\tname=", markername, sep=""))

  # first: select the diploid data (if any) for this marker,
  # to plot into the polyploid histogram:
  if (is.data.frame(diplo)) {
    sampthismarker <- diplo$MarkerName==markername & diploselect &
      !is.na(diplo$ratio)
    drr <- diplo$ratio[sampthismarker]
    names(drr) <- diplo$SampleName[sampthismarker]
    #diploid sample names, should be unique within marker
    #TODO: select the entire rows from diplo, as we do below for the polyploids?
  } else drr <- numeric(0) #for later test on length > 0

  #now select the polyploids for this marker:
  sampthismarker <- (data$MarkerName == markername) & select &
    !is.na(data$ratio) #last term could be due to previous selection for minimum
  #                     intensity
  seldat <- data[sampthismarker,]
  seldat <- seldat[order(seldat$SampleName),]
  #seldat contains only this marker, and only samples with non-NA ratios,
  #and excludes any rows for which select is FALSE,
  #and is ordered by SampleName

  #Now calculate nselsamp: the total number of INTENDED samples
  nselsamp <- length(IntendedSamples)
  if (nrow(seldat) < 10 * ng ||   # in CodomMarker at least 10*ng observations required
      nrow(seldat) < nselsamp * call.threshold) {
    rejected <- TRUE
    log <- append(log, paste(marker, "", "", "rejected: too few data",
                             sep="\t"))
  } else { #test if no samples occur more than once for this marker
    tabsamp <- table(seldat$SampleName)
    if (max(tabsamp) > 1) {
      rejected <- TRUE
      log <- append(log, paste(marker, "", "",
                               "rejected: some samples occur more than once",
                               sep="\t"))
    } else if (is.data.frame(diplo)) {
      allsamp <- diplo$SampleName[diplo$MarkerName==markername]
      tabsamp <- tabulate(allsamp)
      if (max(tabsamp)>1) {
        rejected <- TRUE
        log <- append(log, paste(marker, "", "",
                                 "rejected: some diploid samples occur more than once",sep="\t"))
      }
    }
  }

  mrkpopstruct <- NULL
  if (!rejected) {
    #data and origpopstruct ok, fit initial series of models -> optmodel1
    mrkpopstruct <- getMarkerPopstruct(seldat, origpopstruct, try.HW, ng)
    startmus <- NA
    if (!is.null(startmeans) && !is1NA(startmeans)) {
      i <- which(startmeans[,1] == markername)
      if (length(i) == 1 && !is.na(startmeans[i, 2]))
        startmus <- as.numeric(startmeans[i, -1])
    }
    samppri <- NA
    if (!is.null(samplePriors) && !is1NA(samplePriors)) {
      i <- which(samplePriors[,1] == markername)
      if (length(i) == 1) {
        samppri <- samplePriors[i, -1]
        if (all(is.na(samppri))) samppri <- NA
      }
    }

    priorcomb <- list()
    if (!is.null(origpopstruct$orig_parPriorCols)) {
      mrkparpri <- which(parentalPriors[,1] == markername)
      if (length(mrkparpri) == 1)
        priorcomb <-
          getPriorCombinations(
            origpopstruct$orig_parPriorCols,
            parentalPriors[mrkparpri,, drop=FALSE])
    }

    #Which model to try: new version / phrasing 20160810:
    #
    # Without priors, and in all cases for calculating optmodel2:
    #
    #   for each startmeans (clustered, equidistant or specified) 8 models:
    #   - 1st 4 models: mumodel 1,2,5,6; if 1 pop, p.HW, else use ptype vector
    #   - 2nd 4 models: mumodel 1,2,5,6; all samples together, p.free
    #
    #   effect of try.HW:
    #     if 1 pop: omit 1st 4 models if try.HW FALSE
    #     else: use p.HW or p.free for panels when try.HW TRUE resp FALSE
    #
    # with parentalPriors, when calculating for optmodel1:
    #   modification 20161116: only when more than 1 combinations of
    #   parental priors; with one combination of parental priors
    #   we use the normal approach: clusterd, equidistant and based on
    #   startmeans or parentalPriors and/or samplePriors
    #   (approach of Konrad Zych)
    #   4 mumodels with populations without priors
    #     (i.e. p.start uniform for F1's and parents,
    #     mustart from clustering)
    #   n * 4 mumodels with priors (n combinations of priors)
    #     (i.e. p.start based on F1 seg and parental priors,
    #     mustart calculated from priors and ratios of parents)
    #
    #   effect of try.HW TRUE: use p.HW instead of p.free for panels

    # try several models with ng peaks
    nmodel <- 64 #without parentalPriors: models 1:24 and 41:64 (max);
    #             with priors: models 1:4, 9:12, 17:20, 25:28, 33:36
    #             and 41:64 (max)
    resultnames <-
      getResultnames(ng=ng, pop.parents=origpopstruct$orig_pop.parents)
    mrkresult <- list()
    # mrkresult$stats <- data.frame(rep(NA,nmodel))
    # for (i in 2:length(resultnames)) mrkresult$stats[[i]] <- rep(NA,nmodel)
    # names(mrkresult$stats) <- resultnames
    mrkresult$stats <- makeModelsDF(nrow=nmodel, ng=ng, mrklevels=markernames,
                                    pop.parents=origpopstruct$orig_pop.parents)
    mrkresult$probs <- list()
    mrkresult$Pact <- list()

    # initialise the clustered means once
    # use kmeans with fixed random seed for reproducible results
    # but don't affect the random or non-random number sequence for
    # the calling program:
    if (exists(".Random.seed", .GlobalEnv)) {
      savedseed <- .GlobalEnv$.Random.seed } else savedseed <- NULL
    set.seed(3)
    clusinit <- ClusterInit(asin(sqrt(seldat$ratio)), ng=ng,
                            closedips=TRUE) #returns clus.mu
    #                                        and clus.sd on transformed scale
    if (!is.null(savedseed)) {
      .GlobalEnv$.Random.seed <- savedseed
    } else rm(".Random.seed", envir = .GlobalEnv)
    clus.mu <- sin(clusinit$clus.mu)^2  # transform clus.mu back to 0-1 scale,
    #                                    but keep clus.sd on transformed scale

    #the number of models to calculate with a startmeans
    #(not used for the initial models with parental priors):
    modelcount <- ifelse(try.HW || origpopstruct$orig_nonparents > 1, 8, 4)

    if (length(priorcomb) > 1 && (is.null(startmus) || is1NA(startmus))) {
      #multiple combinations of parentalPriors are specified.
      #In this case we first fit 4 models with the population structure but
      #without the priors information, and next for each (2-4) combination
      #of parentalPriors we fit again 4 models.
      #Since the models should be always the first 4 of a group of 8 (due to
      #the way getModelPopstruct determines which population structure to apply)
      #we fit each time 4 models but then increase modelsfitted by 8.
      preferFrom <- 9 #we prefer models based on priors
      for (model in 1:8) {
        #use popstruct but no parentalPriors, use clustered means instead:
        mrkresult <- MarkerResult(marker=marker, markername=markername,
                                  ratio=seldat$ratio,
                                  model=model,
                                  origPopstruct=origpopstruct,
                                  mrkPopstruct=mrkpopstruct,
                                  mutype=getMutype(model), ng=ng,
                                  maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
                                  mustart=clus.mu, sd=clusinit$clus.sd,
                                  mrkresult=mrkresult, nsamp=nselsamp,
                                  plothist=plotall,
                                  plot.type=plot.type, plot.dir=plot.dir,
                                  modelcount=8)
      }
      modelsfitted <- 8
      tmp_mrkpopstruct <- mrkpopstruct
      for(j in seq_along(priorcomb)) {
        curPriors <- priorcomb[[j]]
        mus <- calcMus(seldat=seldat, priorcomb=priorcomb, priorcomb_index=j,
                       origpopstruct=origpopstruct, mrkpopstruct=mrkpopstruct,
                       samplepriors=samppri, startmus=startmus, ng=ng)
        mrkpopstruct <- get.ps(origpopstruct, mrkpopstruct, curPriors)
        #               modify curr_p.start and curr_ptype
        for (model in 1:4) {
          mrkresult <- MarkerResult(marker=marker, markername=markername,
                                    ratio=seldat$ratio,
                                    model=modelsfitted+model,
                                    origPopstruct=origpopstruct,
                                    mrkPopstruct=mrkpopstruct,
                                    mutype=getMutype(model), ng=ng,
                                    maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
                                    mustart=mus, #sd=clusinit$clus.sd,
                                    mrkresult=mrkresult, nsamp=nselsamp,
                                    plothist=plotall,
                                    plot.type=plot.type, plot.dir=plot.dir,
                                    modelcount=4)
        }
        mrkpopstruct <- tmp_mrkpopstruct
        modelsfitted <- modelsfitted + 8
      }
    } else {
      #0 or 1 combination of parentalPriors are specified
      preferFrom <- 17 #we prefer models based on priors or startmu
      # modelcount models with clustering:
      for (model in 1:8) if (model > 4 || modelcount == 8) {
        mrkresult <- MarkerResult(marker=marker, markername=markername,
                                  ratio=seldat$ratio,
                                  model=model, origPopstruct=origpopstruct,
                                  mrkPopstruct=mrkpopstruct,
                                  mutype=getMutype(model), ng=ng,
                                  maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
                                  mustart=clus.mu, sd=clusinit$clus.sd,
                                  mrkresult=mrkresult, nsamp=nselsamp,
                                  plothist=plotall,
                                  plot.type=plot.type, plot.dir=plot.dir,
                                  modelcount=modelcount)
      }
      #modelcount models without clustering:
      equimu <- seq(asin(sqrt(0.02)), asin(sqrt(0.98)), length.out=ng)
      #         equidistant mu's on transformed scale
      equimu <- sin(equimu)^2 #on non-transformed scale
      for (model in 9:16) if (model > 12 || modelcount == 8) {
        mrkresult <- MarkerResult(marker=marker, markername=markername,
                                  ratio=seldat$ratio,
                                  model=model, origPopstruct=origpopstruct,
                                  mrkPopstruct=mrkpopstruct,
                                  mutype=getMutype(model), ng=ng,
                                  maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
                                  mustart=equimu,
                                  mrkresult=mrkresult, nsamp=nselsamp,
                                  plothist=plotall,
                                  plot.type=plot.type, plot.dir=plot.dir,
                                  modelcount=modelcount)
      }
      modelsfitted <- 16
      #modelcount models with startmeans and/or samplePriors:
      #(priorcomb_index is ignored if length(priorcomb)==0)
      sm <- calcMus(seldat=seldat, priorcomb=priorcomb, priorcomb_index=1,
                    origpopstruct=origpopstruct, mrkpopstruct=mrkpopstruct,
                    samplepriors=samppri, startmus=startmus, ng=ng)
      if (!is1NA(sm)) {
        tmp_mrkpopstruct <- mrkpopstruct
        if (length(priorcomb) == 1)
          mrkpopstruct <- get.ps(origpopstruct, mrkpopstruct, priorcomb[[1]])
        for (model in 17:24) if (model > 20 || modelcount == 8) {
          mrkresult <- MarkerResult(marker=marker, markername=markername,
                                    ratio=seldat$ratio,
                                    model=model, origPopstruct=origpopstruct,
                                    mrkPopstruct=mrkpopstruct,
                                    mutype=getMutype(model), ng=ng,
                                    maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
                                    mustart=sm,
                                    mrkresult=mrkresult, nsamp=nselsamp,
                                    plothist=plotall,
                                    plot.type=plot.type, plot.dir=plot.dir,
                                    modelcount=modelcount)
        }
        mrkpopstruct <- tmp_mrkpopstruct
        modelsfitted <- 24
      }
    } # length(priorcomb) <= 1
    #Now from the fitted models we select the best: optmodel1
    selcritTolerance <- 5 #TODO: make a parameter of fitOneMarker?
    mrkresult$stats <-
      calcSelcrit(mrkresult$stats,
                  rep(TRUE, nrow(mrkresult$stats)), sd.target)
    optmodel1 <- which.min(mrkresult$stats$selcrit)
    if (length(optmodel1) == 0 || is.na(optmodel1)) {
      rejected <- TRUE
      optmodel2 <- optmodel1 <- NA
      log <- append(log,paste(marker,"optmodel1=NA", "",
                              "rejected: optmodel1=NA", sep="\t"))
    } else {
      #optmodel1 not rejected
      if (modelsfitted > preferFrom) {
        optmodel1a <-
          which.min(mrkresult$stats$selcrit[preferFrom:modelsfitted]) +
          preferFrom - 1
        #optmodel1a is the best among the preferred models
        if (length(optmodel1a) == 1 && !is.na(optmodel1a) &&
                 mrkresult$stats$selcrit[optmodel1a] -
                 mrkresult$stats$selcrit[optmodel1] < selcritTolerance)
        optmodel1 <- optmodel1a
      }
    }
  } # optmodel 1 fitted or rejected TRUE

  if (!rejected) {
    #optmodel1 != NA

    #optmodel.ps <- getModelPopstruct(optmodel1, origpopstruct, mrkpopstruct)
    #opt.mutype <- getMutype(optmodel1)
    #opt_ptype <- getPtype(optmodel1)
    logline <- paste(marker, "\toptmodel1=", optmodel1, "\t",
                     mrkresult$stats$model[optmodel1], sep="")

    #First we check for two situations indicating a possibly wrong fit:
    # - the nulli or <ploidy>plex peak is very small while the simplex or
    #   <ploidy-1>plex is big ("shifts"), or
    # - there are small peaks surrounded by higher peaks ("dips")
    #In these cases we run the 8 models again with one or more new mustart
    # check for dips and shifts and calculate corrected start.mus:
    nwmu <- shift_mus(ng=ng, mean_ratio=mean(asin(sqrt(seldat$ratio))),
                      stats=mrkresult$stats[optmodel1,],
                      probs=mrkresult$probs[[optmodel1]])
    #fit modelcount models with each of the new mustart:
    if (length(nwmu) > 0) {
      modelsfitted <- 40 # above any of the earlier models
      logline <- paste(logline, nwmu[[length(nwmu)]], sep="\t")
      for (munr in 1:(length(nwmu)-1)) {
        for (model in (modelsfitted+1):(modelsfitted+8))
          if (model > modelsfitted+4 || modelcount == 8) {
          mrkresult <- MarkerResult(marker=marker, markername=markername,
                                    ratio=seldat$ratio,
                                    model=model, origPopstruct=origpopstruct,
                                    mrkPopstruct=mrkpopstruct,
                                    mutype=getMutype(model), ng=ng,
                                    maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
                                    mustart=nwmu[[munr]],
                                    mrkresult=mrkresult, nsamp=nselsamp,
                                    plothist=plotall,
                                    plot.type=plot.type, plot.dir=plot.dir,
                                    modelcount=modelcount)
        } # for model
        modelsfitted <- modelsfitted + 8
      } # for munr
    } #length(nwmu) > 0
    log <- append(log, logline)

    #now select the best model (optmodel2) taking dip into account:
    if (dip.filter > 0) {
      dipok <- !(is.na(mrkresult$stats$dip) | mrkresult$stats$dip)
      #        only rows with model without dip
      mrkresult$stats <- calcSelcrit(mrkresult$stats, dipok, sd.target)
      optmodel2 <- which.min(mrkresult$stats$selcrit[1:modelsfitted])
      optmodel2a <- which.min(mrkresult$stats$selcrit[41:max(41,modelsfitted)])
      if (length(optmodel2) == 0 || is.na(optmodel2)) {
        mrkresult$stats <- calcSelcrit(mrkresult$stats,
                                       rep(TRUE, nrow(mrkresult$stats)),
                                       sd.target) #all rows
        optmodel2 <- which.min(mrkresult$stats$selcrit[1:modelsfitted])
        optmodel2a <- which.min(mrkresult$stats$selcrit[41:max(41,modelsfitted)])
      }
      #now we have optmodel2 NA if no model converged; else optmodel2 is best
      #model without dip if these occur, else best model overall (with dip)
    } else {
      #dip.filter==0, select best model without considering dip:
      mrkresult$stats <- calcSelcrit(mrkresult$stats,
                                     rep(TRUE, nrow(mrkresult$stats)),
                                     sd.target) #all rows
      optmodel2 <- which.min(mrkresult$stats$selcrit[1:modelsfitted])
      optmodel2a <- which.min(mrkresult$stats$selcrit[41:max(41,modelsfitted)])
    }
    if (length(optmodel2)==0 || is.na(optmodel2)) {
      optmodel2 <- NA
      # reject this marker
      rejected <- TRUE
      log <- append(log, paste(marker, "optmodel2=NA", "",
                               "marker rejected: optmodel2=NA", sep="\t"))
    } else {
      # not rejected, optmodel2 found
      if (length(optmodel2a) == 1 &&
          mrkresult$stats$selcrit[40 + optmodel2a] - mrkresult$stats$selcrit[optmodel2]
          < selcritTolerance) {
        # a model obtained by shifting the start means is (almost) as good as an
        #original model; in this case we prefer the new model:
        optmodel2 <- 40 + optmodel2a
      }
      #opt_ptype <- getPtype(optmodel2)
      #opt_mutype <- getMutype(optmodel2)
      logline <- paste(marker, "\toptmodel2=", optmodel2, "\t",
                       mrkresult$stats$model[optmodel2], sep="")
      # geno <- maxgeno for samples above p.threshold, else NA:
      mrkresult$probs[[optmodel2]]$geno <- mrkresult$probs[[optmodel2]]$maxgeno
      sel.maxP <- mrkresult$probs[[optmodel2]]$maxP < p.threshold
      mrkresult$probs[[optmodel2]]$geno[sel.maxP] <- NA
      probs <- makeprobs(marker=marker, markername=markername,
                         samplenames=samplenames,
                         modelname=mrkresult$stats$model[optmodel2],
                         select=select,
                         markerlines=(data$MarkerName == markername),
                         seldat=seldat,
                         resultprobs=mrkresult$probs[[optmodel2]],
                         dat=data)
      # check several reasons to still reject optmodel2:
      # (in these cases we have an optmodel2 to output and plot,
      #  but it is rejected)
      if (sum(!is.na(mrkresult$probs[[optmodel2]]$geno)) <
          call.threshold*nselsamp) {
        # note that we compare to the samples with select=TRUE,
        #not to the total number of samples
        rejected <- TRUE
        mess <- "rejected: less than minimum number of samples called"
        logline <- paste(logline, mess, sep="\t")
        mrkresult$stats$message[optmodel2] <- mess
      } else if (mrkresult$stats$sdtrans0[optmodel2] > sd.threshold) {
        rejected <- TRUE
        mess <- "rejected: sd>sd.threshold"
        logline <- paste(logline, mess, sep="\t")
        mrkresult$stats$message[optmodel2] <- mess
      } else  if (mrkresult$stats$dip[optmodel2] && dip.filter==2) {
        rejected <- TRUE
        mess <- "rejected: no models without dip"
        logline <- paste(logline, mess, sep="\t")
        mrkresult$stats$message[optmodel2] <- mess
      } else { #optmodel2 not yet rejected
        gsamp <- hist(mrkresult$probs[[optmodel2]]$geno,
                      breaks=(0:ng)-0.5, plot=FALSE)$counts #based on called samples
        if (max(gsamp) / sum(gsamp) > peak.threshold) {
          rejected <- TRUE
          mess <- "rejected: more than maximum fraction of samples in one peak"
          logline <- paste(logline, mess, sep="\t")
          mrkresult$stats$message[optmodel2] <- mess
        }
      }
      if (!rejected && mrkresult$stats$dip[optmodel2] && dip.filter==1) {
        #if dip.filter==2 marker was already rejected; if 0 there might be
        #a model without dip but worse selcrit, now we just give a message:
        logline <- paste(logline, "no models without dip", sep="\t")
      }
      log <- append(log, logline)
    }
  } # optmodel1 not NA

  # Now 3 situations possible:
  # rejected FALSE, optmodel2 found: produce full output and
  #                                  plot fitted model if plotfitted
  # rejected TRUE, optmodel2 found: produce full output but set all probs$geno
  #                                 to NA; plot fitted model if plotfitted
  # rejected TRUE, optmodel2 NA: produce limited output and plot histogram(s)
  #                              if plotfitted

  if (plotfitted) {
    if (is.na(optmodel2)) {
      plotfitted.fname <- paste(plot.dir, "rejected_", mrknr, " ",
                                markernames[marker], sep="")
      startPlot(plot.type=plot.type, plot.filename=plotfitted.fname)
      par(cex=1, mex=1)
      if (is.character(origpopstruct) || is.null(mrkpopstruct)) {
        #invalid population structure
        errorplot(main=paste(marker, markername),
                  message="too few observations or invalid input")
      } else {
        drawRejectedPlot(rr=seldat$ratio,
                         #rrpop=origpopstruct$orig_population,
                         rrpop=mrkpopstruct$curr_pop,
                         drr=drr,
                         ploidy=ploidy,
                         pop.parents=origpopstruct$orig_pop.parents,
                         nonparents=origpopstruct$orig_nonparents,
                         maintitle=paste(marker, markername, "(no fit)"))
      }
      savePlot()
    } else {
      #plotfitted, optmodel2 found
      plotfitted.fname <- paste(plot.dir,
                                ifelse(rejected, "rejected_", ""),
                                mrknr, " ",
                                markernames[marker], sep="")
      startPlot(plot.type=plot.type, plot.filename=plotfitted.fname)
      nfirst <- which(names(mrkresult$stats) == "mutrans0")
      psi <- list()
      psi$mu <- as.numeric(mrkresult$stats[optmodel2, nfirst:(nfirst+ng-1)])
      psi$sigma <- as.numeric(mrkresult$stats[optmodel2,
                                              (nfirst+ng):(nfirst+2*ng-1)])
      optmodel.ps <- getModelPopstruct(optmodel2, origpopstruct,
                                       mrkpopstruct, ng)
      npop <- nrow(optmodel.ps$pop.parents)
      psi$p <- matrix(
        as.numeric(
          mrkresult$stats[optmodel2,
                          (nfirst+2*ng):(nfirst+(2+npop)*ng-1)]),
        nrow=npop, byrow=TRUE)
      drawFittedPlot(rr=seldat$ratio, rrpop=optmodel.ps$population,
                     drr=drr, pop.parents=optmodel.ps$pop.parents,
                     nonparents=optmodel.ps$nonparents,
                     geno=mrkresult$probs[[optmodel2]]$geno,
                     psi, maintitle=paste(marker, markername))
      savePlot()
    }
  } #plotfitted

  # Finally produce the output list:

  result <- list(log=log)
  result$modeldata <- NA
  result$allmodeldata <- NA
  result$scores <- NA
  result$diploscores <- NA
  if (is.null(mrkresult)) {
    #only if no fitting attempted (invalid data, invalid popstruct,
    #too few or duplicate samples)
    # result$modeldata <- data.frame(
    #   marker=marker,
    #   markername=markername,
    #   m=NA,
    #   model="none",
    #   nsamp=nselsamp,
    #   nsel=ifelse(is.null(seldat), NA, nrow(seldat)))
    # for (i in 7:length(resultnames)) result$modeldata[[i]] <- NA
    # names(result$modeldata) <- resultnames
    result$modeldata <- makeModelsDF(nrow=1, ng=ng, mrklevels=markernames,
                                     pop.parents=origpopstruct$orig_pop.parents)
    result$modeldata$marker <- marker
    result$modeldata$MarkerName <- markername
    result$modeldata$nsamp <- nselsamp
    result$modeldata$nsel <- ifelse(is.null(seldat), NA, nrow(seldat))
    result$modeldata$message <- "rejected: too few observations or invalid input"
    result$allmodeldata <- result$modeldata
  } else { #valid mrkresult
    result$allmodeldata <- mrkresult$stats[!is.na(mrkresult$stats$m),]
    if (is.na(optmodel2)) {
      result$modeldata <- mrkresult$stats[1,]
      result$modeldata$m <- NA  #no model selected
      result$modeldata$model <- NA  #no model selected
      for (i in which(names(mrkresult$stats) == "npar"):
           (length(mrkresult$stats)-1))
        result$modeldata[1, i] <- NA
      result$modeldata$message <- paste("tried", nrow(result$allmodeldata),
                                        "models, no convergence")
      result$scores <- NA
      result$diploscores <- NA
    } else {
      #optmodel2 found (possibly rejected)
      #from 2015-12-21 we output all info also for rejected markers
      #(with geno NA in that case):
      result$modeldata <- mrkresult$stats[optmodel2,]
      rownames(probs) <- NULL
      result$scores <- probs
      if (rejected) result$scores$geno <- NA
      if (is.data.frame(diplo)) { #not NA
        # if diplo exists, also export scores for all samples in diplo:
        mu0 <- which(names(result$modeldata) == "mutrans0")
        sd0 <- which(names(result$modeldata) == "sdtrans0")
        #p0 <- which(names(result$modeldata) == "P0")
        psi <- list(
          mu=as.numeric(result$modeldata[, mu0:(mu0+ng-1)]),
          sigma=as.numeric(result$modeldata[, sd0:(sd0+ng-1)]),
          #p=result$modeldata[, p0:(p0+ng-1)])
          p=rep(1/ng, ng)) #no prior for "diploids" (can be any ploidy actually)
        diploprobs <- EMGaussExp.vectorized(asin(sqrt(drr)), psi) #matrix,
        #             per diplo sample ng numbers:
        #             probs that sample belongs to each peak
        #calculate maxgeno, maxP, geno:
        maxPna <- apply(diploprobs, 1, max) #find the maxP of each row of p-values
        maxP <- maxPna[!is.na(maxPna)] # omit missing values
        probs <- data.frame(diploprobs)
        rownames(probs) <- names(drr) #note: does not include samples rejected
        #                  because select=FALSE
        names(probs) <- paste("P", 0:(ng-1), sep="")
        #list for each sample the maximum p value and the maxgeno:
        probs$maxgeno <- max.col(diploprobs) - 1  # for every row (sample) the
        #                max peak (0..ploidy instead of 1..ng)
        probs$maxP <- maxPna # for every row the P at the max
        #add extra columns and add rows for samples not in drr:
        result$diploscores <-
          makeprobs(marker=marker, markername=markername,
                    samplenames=diplonames,
                    modelname=mrkresult$stats$model[optmodel2],
                    select=diploselect,
                    markerlines=(diplo$MarkerName == markername),
                    seldat=data.frame(SampleName=names(drr), ratio=drr),
                    resultprobs=probs,
                    dat=diplo)
        rownames(result$diploscores) <- NULL
        result$diploscores$geno <- NA
        if (!rejected) {
          sel.maxP <- which(result$diploscores$maxP >= p.threshold)
          result$diploscores$geno[sel.maxP] <-
            result$diploscores$maxgeno[sel.maxP]
        }
      } #is.data.frame(diplo)
    }
  }
  result
} # fitOneMarker


# makeprobs takes a probs data frame that for each selected, non-NA ratio in data or diplo
# contains the ng P-values, maxgeno, maxP and geno
# and makes a similar data frame but now with rows for all samples (not only the selected, non-NA ones)
# and with additional columns for marker number and name, sample name, model name,
# select and ratio
makeprobs <- function(marker, markername, samplenames, modelname, select,
                      markerlines, seldat, resultprobs, dat) {
  nsamp <- length(samplenames)
  probs <- data.frame(marker=rep(marker, nsamp),
                      MarkerName=markername,
                      SampleName=samplenames,
                      #model=modelname,
                      #nsamp=nsamp,
                      #select=0,
                      ratio=NA_real_)
  #probs has a row for each sample while rr has only the selected non-NA
  #sample ratios
  #sel <- select[dat$MarkerName==markername] #note that selection is less
  #      stringent than that of seldat, so sel may be longer than seldat
  #names(sel) <- dat$SampleName[dat$MarkerName==markername]
  #m <- match(probs$SampleName, names(sel))
  #probs$select <- as.numeric(sel[m])
  m <- match(probs$SampleName, seldat$SampleName) #has NA for samples not in seldat
  probs$ratio <- seldat$ratio[m] #has NA for samples where m=NA
  probs <- cbind(probs, resultprobs[m, 1:length(resultprobs)]) #idem
  probs
} #makeprobs

shift_mus <- function(ng, mean_ratio, stats, probs) {
  #This function generates up to 3 alternative configurations of mu's
  #to be used as start.mu, attempting to correct for dips and/or shifts
  #in optmodel1
  # ng: number of dosage classes = ploidy + 1
  # mean_ratio: mean of the ratios of all samples
  # stats: mrkresult$stats[optmodel1,]
  # probs: mrkresult$probs[[optmodel1]]
  # result: a list of length 0 if no shift, else a list of length 2 to 4
  #         where the last element is a comment for the log file
  #         and all earlier elements are vectors of ng mu's on the
  #         non-transformed scale
  result <- list()
  remark_dips <- ""
  shifts <- FALSE
  #get a vector of the model mu's on the transformed scale:
  mu0 <- which(names(stats)=="mutrans0")
  mu <- unlist(stats[, mu0:(mu0+ng-1)])
  gsamp <- hist(probs$maxgeno, breaks=(0:ng) - 0.5, plot=FALSE)$counts
  #        based on maxgeno: all samples
  origpeaks <- 1:ng
  #check for dips:
  if (stats$dip) {
    dip_threshold <- 1/ng/4 # 25% of peak height if all peaks same height;
    #                         same as in MarkerResult
    tops <- which(gsamp >= dip_threshold * sum(gsamp)) #at least 1:
    #       tops are all peak nrs (1:ng) with fraction of samples above
    #       dip_threshold
    #       a valley is the series of peaks between two tops
    if (max(tops) - min(tops) >= length(tops)) {
      #there are small peaks within the range containing all larger peaks
      #identify all the valleys:
      topnrs_before_valley <-
        which(tops[2:length(tops)] > tops[1:(length(tops)-1)] + 1)
      startvalleys <- tops[topnrs_before_valley] + 1
      endvalleys <- tops[topnrs_before_valley + 1] -1
      #first and last gsamp of each valley
      #from each valley identify the lowest peak
      dels <- numeric(length(startvalleys))
      for (i in seq_along(startvalleys)) {
        dels[i] <-
          startvalleys[i] - 1 + which.min(gsamp[startvalleys[i]:endvalleys[i]])
      }
      #next: before we delete these peaks, we add their samples to the neighbours:
      #(and we know that the neighbours will not be deleted)
      for (d in seq_along(dels)) {
        #we roughly divide the samples and model freq 50/50 over the two neighbours
        #(this can be done more sophisticated if needed)
        gsamp[dels[d]-1] <- gsamp[dels[d]-1] + gsamp[dels[d]]/2
        gsamp[dels[d]+1] <- gsamp[dels[d]+1] + gsamp[dels[d]]/2
        #p[dels[d]-1] <- p[dels[d]-1] + p[dels[d]]/2
        #p[dels[d]+1] <- p[dels[d]+1] + p[dels[d]]/2
      }
      #delete the dips:
      gsamp <- gsamp[-dels]
      mu <- mu[-dels]
      #p <- p[-dels]
      origpeaks <- origpeaks[-dels]
      remark_dips <- paste("remove",length(dels),"dip(s)")
    }
  }
  #now we have deleted the deepest dip from each valley (if any)
  #and we have ng or fewer (but at least 2) remaining peaks.
  #Any small outer peaks are still present, and the original peak numbers (1..ng)
  #of the remaining peaks are in origpeaks.

  #next, we check if there is evidence that low peaks at one end
  #should be deleted:
  nonEmpty <- which(gsamp > 0.02 * sum(gsamp)) #the current non-almost empty peaks
  lne <- length(nonEmpty)
  if (lne == 1) {
    #monomorphic, and there was no internal valley from which dips could have
    #been removed.
    if (!(origpeaks[nonEmpty] %in% c(1, ng))) {
      #(else the peak was already at dosage 0 or ploidy, no new model needed)
      #make dosage 0 or ploidy:
      result[[1]] <- numeric(ng)
      if (mean_ratio < pi/4) {
        result[[1]][1] <- mu[nonEmpty]
        result[[1]][2:ng] <- mu[nonEmpty] +
          (1:(ng-1)) * (pi/2 - mu[nonEmpty]) / (ng-1)
      } else {
        result[[1]][ng] <- mu[nonEmpty]
        result[[1]][1:(ng-1)] <- (0:(ng-2)) * mu[nonEmpty] / (ng-1)
      }
      shifts <- TRUE
    }
  } else if (lne == 2) {
    #two non-almost empty peaks
    if (sum(origpeaks[nonEmpty]) != 3 &&
        sum(origpeaks[nonEmpty]) != 2*ng-1) {
      #the two peaks were not the first two or the last two
      if (mean_ratio <= pi/4 && gsamp[nonEmpty[1]] >= 0.8*gsamp[nonEmpty[2]]) {
        #put the two peaks only at the low side
        result[[1]] <- numeric(ng)
        result[[1]][1:2] <- mu[nonEmpty]
        result[[1]][3:ng] <- mu[nonEmpty][2] + (1:(ng-2)) * (pi/2-mu[nonEmpty[2]])/(ng-2)
      } else if (mean_ratio >= pi/4 && gsamp[nonEmpty[2]] >= 0.8*gsamp[nonEmpty[1]]) {
        #put the two peaks only at the high side
        result[[1]] <- numeric(ng)
        result[[1]][(ng-1):ng] <- mu[nonEmpty]
        result[[1]][1:(ng-2)] <-  (0:(ng-3)) * mu[nonEmpty[1]]/(ng-2)
      } else {
        #put the two peaks at the low and high side
        result[[2]] <- numeric(ng)
        result[[2]][(ng-1):ng] <- mu[nonEmpty]
        result[[2]][1:(ng-2)] <-  (0:(ng-3)) * mu[nonEmpty[1]]/(ng-2)
        result[[1]] <- numeric(ng)
        result[[1]][1:2] <- mu[nonEmpty]
        result[[1]][3:ng] <- mu[nonEmpty][2] + (1:(ng-2)) * (pi/2-mu[nonEmpty[2]])/(ng-2)
      }
      shifts <- TRUE
    }
  } else if (lne > 0 && lne < ng) {
    #three or more non-almost empty peaks
    #while for 1 and 2 peaks we accepted the original version if the peaks were
    #already at one of the ends, here we first determine the desired
    #combinations and then discard the one that was already done (if any).
    if (mean_ratio <= pi/4 && gsamp[nonEmpty[1]] >= gsamp[nonEmpty[lne]]) {
      #put the (lne) peaks only at the low side
      if (sum(origpeaks[nonEmpty]) != sum(1:lne)) {
        result[[1]] <- numeric(ng)
        result[[1]][1:lne] <- mu[nonEmpty]
        result[[1]][(lne+1):ng] <-
          mu[nonEmpty][lne] + (1:(ng-lne)) * (pi/2 - mu[nonEmpty[lne]]) / (ng-lne)
        shifts <- TRUE
      }
    } else if (mean_ratio >= pi/4 && gsamp[nonEmpty[lne]] >= gsamp[nonEmpty[1]]) {
      #put the (lne) peaks only at the high side
      if (sum(origpeaks[nonEmpty]) != sum((ng-lne+1):ng)) {
        result[[1]] <- numeric(ng)
        result[[1]][(ng-lne+1):ng] <- mu[nonEmpty]
        result[[1]][1:(ng-lne)] <-  (0:(ng-lne-1)) * mu[nonEmpty[1]] / (ng-lne)
        shifts <- TRUE
      }
    } else {
      #put the (lne) peaks at the low and high side
      if (sum(origpeaks[nonEmpty]) != sum(1:lne)) {
        result[[1]] <- numeric(ng)
        result[[1]][1:lne] <- mu[nonEmpty]
        result[[1]][(lne+1):ng] <-
          mu[nonEmpty][lne] + (1:(ng-lne)) * (pi/2 - mu[nonEmpty[lne]]) / (ng-lne)
        shifts <- TRUE
      }
      if (sum(origpeaks[nonEmpty]) != sum((ng-lne+1):ng)) {
        i <- length(result) + 1
        result[[i]] <- numeric(ng)
        result[[i]][(ng-lne+1):ng] <- mu[nonEmpty]
        result[[i]][1:(ng-lne)] <-  (0:(ng-lne-1)) * mu[nonEmpty[1]]  /(ng-lne)
        shifts <- TRUE
      }
    }
    #If there are ng-2 non-(almost) empty peaks and the distribution is
    #more or less symmetric and unimodal, a central position is also likely;
    #don't do this for tetraploids, where 1:8:18:8:1 is likely to have at least
    #one of the extreme peaks (each 2.8%) detected and an extra model for
    #1:2:1 segregation of not wanted
    if (ng > 5 && lne == ng-2 && sum(origpeaks[nonEmpty]) != sum(2:(ng-1))) {
      mx <- which.max(gsamp[nonEmpty])
      if (mx == (lne+1) / 2 && #central peak
          gsamp[mx] < (gsamp[mx-1] + gsamp[mx+1])) {#excludes most 1:2:1 and 1:4:1
        i <- length(result) + 1
        result[[i]] <- numeric(ng)
        result[[i]][1] <- 0
        result[[i]][ng] <- pi/2
        result[[i]][2:(ng-1)] <- mu[nonEmpty]
        shifts <- TRUE
      }
    }
  }
  if (length(result) > 0) {
    #back-transform mu's to original scale:
    for (i in seq_along(result)) {
      result[[i]] <- sin(result[[i]]) * sin(result[[i]])
    }
    #add remark for log as last element of result
    i <- length(result) + 1
    if (remark_dips != "") {
      if (shifts) {
        result[[i]] <- paste(remark_dips, "and apply shift(s)")
      } else result[[i]] <- remark_dips
    } else result[[i]] <- "apply shift(s)"
  }
  result
} #shift_mus

checkCodomMarkerData <- function(y, ng, pop, pop.parents, ptype) {
  #returns an error message or ""
  if (ng < 1) return("ng < 1")
  if (sum(!is.na(y)) < 10 * ng)
    return("y should have at least 10 * ng non-NA elements")
  if (any(y < 0 | y > 1, na.rm=TRUE))
    return("all y values should be between 0 and 1 (inclusive)")
  if (length(pop) != length(y) || anyNA(pop) ||
      !(class(pop)[1] %in% c("factor", "character", "integer")))
    return("pop must be a vector matching y without missing values")
  if (!is.matrix(pop.parents) ||
      nrow(pop.parents) == 0 || ncol(pop.parents) != 2 ||
      any(!(as.integer(pop.parents) %in% c(NA, 1:nrow(pop.parents))))) {
    return("invalid pop.parents")
  }
  if (is.null(row.names(pop.parents)))
    row.names(pop.parents) <- 1:nrow(pop.parents)
  for (p in 1:nrow(pop.parents)) {
    if (xor(is.na(pop.parents[p, 1]), is.na(pop.parents[p, 2])) ||
        (!is.na(pop.parents[p, 1]) && min(pop.parents[p,]) <= p) ) {
      #last line checks that parents occur below their progeny
      return("invalid pop.parents")
    }
  }
  #change factor or character pop to integer
  nwpop <- FALSE
  if (is.factor(pop)) pop <- as.character(pop)
  if (is.character(pop)) {
    if (!all(unique(pop) %in% row.names(pop.parents)))
      return("population ID's in pop don't match rownames of pop.parents")
    pop <- match(pop, rownames(pop.parents))
    nwpop <- TRUE
  } else {
    if (!all(pop %in% 1:nrow(pop.parents)))
      return("pop doesn't match rows of pop.parents")
  }
  nwPtype <- FALSE
  if (is1NA(ptype)) {
    #set the ptype according to pop.parents: p.F1 for F1's, p.free
    #for all other populations; check later for missing parental values
    ptype <- rep("p.free", nrow(pop.parents))
    ptype[!is.na(pop.parents[,1])] <- "p.F1"
    nwPtype <- TRUE
  }
  if (length(ptype) != nrow(pop.parents))
    return("length(ptype) different from nrow(pop.parents)")
  #for populations with no data, ptype must be set to "p.nodata"
  #and F1 parents may not be all missing
  if (ng %% 2 == 0 && "p.F1" %in% ptype)
    return("ptype p.F1 not allowed with odd ploidy")
  for (popi in nrow(pop.parents):1) {
    if (ptype[popi] == "p.nodata" && any(!is.na(y[pop == popi])))
      return(paste("population", popi,
                   "has non-missing values but its ptype is p.nodata"))
    if (ptype[popi] != "p.nodata" && all(is.na(y[pop == popi])))
      { ptype[popi] <- "p.nodata"; nwPtype <- TRUE }
    if (ptype[popi] == "p.F1") {
      if (anyNA(pop.parents[popi,]))
        return(paste("population", popi,
                     "has ptype p.F1 but one or both parents are unspecified"))
      #check that both parents have samples - they have been checked earlier:
      if ("p.nodata" %in% ptype[pop.parents[popi,]])
        return(paste("population", popi,
                     "has ptype p.F1 but one or both parents have no data"))
    } else {
      #ptype[popi] not p.F1
      if (sum(!is.na(pop.parents[popi,])) > 0)
        return(paste("population", popi,
                     "doesn't have ptype p.F1 but one or two parents are specified"))
    }
  }
  result <- list()
  if (nwpop) result$pop <- pop
  if (nwPtype) result$ptype <- ptype
  result
} #checkCodomMarkerData


#'@title Function to fit a multiple mixture model to a vector of signal ratios
#'of a single bi-allelic marker
#'
#'@description This function fits a specified mixture model to a vector of
#'signal ratios of multiple samples for a single bi-allelic marker.
#'Returns a list with results from the fitted mixture model.
#'
#'@usage CodomMarker(y, ng, pop.parents=matrix(c(NA,NA), nrow=1),
#'pop=rep(1, length(y)), mutype=0, sdtype="sd.const", ptype=NA,
#'clus=TRUE, mu.start=NA, sd=rep(0.075, ng), p=NA,
#'maxiter=500, maxn.bin=200, nbin=200, plothist=TRUE, nbreaks=40,
#'maintitle=NULL, closeScreen=TRUE, fPinfo=NA)
#'
#'@param y the vector of signal ratios (each value is from one sample,
#'vector y contains the values for one marker). All values must be between
#'0 and 1 (inclusive), NAs are not allowed. The minimum length of y is 10*ng.
#'@param ng the number of possible genotypes (mixture components) to be fitted:
#'one more than the ploidy of the samples.
#'@param pop.parents a matrix with 2 columns and 1 row per population;
#'the cells contain the row numbers of the parental populations in case of an
#'F1 and NA otherwise. The rows must be sorted such that all F1s occur above
#'their parental populations. By default 1 row with elements NA, i.e. all
#'samples belong to a single non-F1 population. If parameter pop is a factor or
#'character vector, its levels or elements must correspond to the rownames of
#'pop.parents.
#'@param pop an integer vector specifying the population to which each sample
#'in y belongs. All values must index rows of pop.parents. By default a vector
#'of 1's, i.e. all samples belong to a single non-F1 population. Alternatively
#'pop can be a factor or character vector of which the levels or elements
#'match the rownames of pop.parents
#'@param mutype	an integer in 0:6; default 0. Describes how to fit the means of the
#'components of the mixture model: with mutype=0 the means are not constrained,
#'requiring ng degrees of freedom. With mutype in 1:6 the means are constrained
#'based on the ng possible allele ratios according to one of 6 models;
#'see Details.
#'@param sdtype	one of "sd.const", "sd.free", "sd.fixed"; default "sd.const".
#'Describes how to fit the standard deviations of the components of the mixture
#'model: with "sd.const" all standard deviations (on the transformed scale)
#'are equal (requiring 1 degree of freedom); with "sd.free" all standard
#'deviations are  fitted separately (ng d.f.); with "sd.fixed" all sd's ON
#'THE TRANSFORMED SCALE are equal to parameter sd (0 d.f.).
#'@param ptype a character vector of length nrow(pop.parents) containing for
#'each population one of "p.free", "p.fixed", "p.HW" or "p.F1". The
#'default NA is interpreted as "p.F1" for F1 populations and "p.free" for all
#'other populations; this is not necessarily the best choice for GWAS panels
#'where "p.HW" may be more appropriate. Describes per population how to fit
#'the mixing proportions of the components of the mixture model:
#'with "p.free", the proportions are not constrained (and require ng-1 degrees
#'of freedom per population); with "p.fixed" the proportions given in
#'parameter p are fixed; with "p.HW" the proportions are calculated per
#'population from an estimated allele frequency, requiring only 1 degree of
#'freedom per population; with "p.F1" polysomic (auto-polyploid) F1
#'segregation ratios are calculated based on the fitted dosages of the F1
#'parents and require no extra d.f.
#'@param clus boolean. If TRUE, the initial means and standard deviations are
#'based on a kmeans clustering of all samples into ng or fewer groups. If FALSE,
#'the initial means are equally spaced on the transformed scale between the
#'values corresponding to 0.02 and 0.98 on the original scale and the initial
#'standard deviations are 0.075 on the transformed scale.
#'@param mu.start	vector of ng values. If present, gives the start values of mu
#'(the means of the mixture components) on the original (untransformed) scale.
#'Must be strictly ascending (mu[i] > mu[i-1]) between 0 and 1 (inclusive).
#'Overrides the start values determined by clus TRUE or FALSE.
#'@param sd vector of ng values. If present, gives the initial (or fixed,
#'if sd.fixed is TRUE) values of sd (the standard deviations of the mixture
#'components) ON THE TRANSFORMED SCALE. Overrides the start values determined
#'by clus TRUE or FALSE.
#'@param p a matrix of nrow(pop.parents) rows and ng columns, each row summing
#'to 1. If present, specifies the initial (or fixed, for populations where
#'ptype is "p.fixed") mixing proportions of the mixture model components.
#'@param maxiter a single integer: the maximum number of times the nls function
#'is called (0 = no limit, default=500).
#'@param maxn.bin a single integer, default=200: if the length of y is larger
#'than max.nbin the values of y (after arcsine square root transformation) are
#'binned (i.e. the range of y (0 to pi/2) is divided into nbin bins of equal
#'width and the number of y values in each bin is used as the weight of the
#'midpoints of each bin). This results in significant speed improvement with
#'large numbers of samples without noticeable effects on model fitting.
#'@param nbin a single integer, default=200: the number of bins (see maxn.bin).
#'@param plothist if TRUE (default) a histogram of y is plotted with the fitted
#'distributions superimposed
#'@param nbreaks number of breaks (default 40) for plotting the histogram;
#'does not have an effect on fitting the mixture model.
#'@param maintitle string, used as title in the plotted histogram.
#'@param closeScreen logical, only has an effect if plothist is TRUE.
#'closeScreen should be TRUE (default) unless CodomMarker will plot on a
#'device that is managed outside CodomMarker.
#'@param fPinfo NA (default), for internal use only. Prevents unneeded checking
#'and recalculation of input parameters when called from fitOneMarker.
#'@details This function takes as input a vector of ratios of the signals of
#'two alleles (a and b) at one genetic marker locus (ratios as b/(a+b)), one for
#'each sample, and fits a mixture model with ng components (for a tetraploid
#'species: ng=5 components representing the nulliplex, simplex, duplex, triplex
#'and quadruplex genotypes). Ideally these signal ratios should reflect the
#'possible allele ratios (for a tetraploid: 0, 0.25, 0.5, 0.75, 1) but in real
#'life they show a continuous distribution with a number of more or less clearly
#'defined peaks. The samples can represent multiple populations, each with
#'their own segregation type (polysomic F1 ratios, Hardy-Weinberg ratios or
#'free ratios). Multiple arguments specify what model to fit and with what
#'values the iterative fitting process should start.\cr
#'Parameter mutype determines how the means of the mixture model components are
#'constrained based on the possible allele ratios, as follows
#'\describe{
#'  \item{0}{all means are fitted without restrictions (ng parameters)}
#'  \item{1}{a basic model assuming that both allele signals have a linear
#'           response to the allele dosage; one parameter for the ratio of the
#'           slopes of the two signal responses, and two parameters for the
#'           background levels (intercepts) of both signals (total 3
#'           parameters)}
#'  \item{2}{as 1, but with the same background level for both signals
#'           (2 parameters)}
#'  \item{3}{as 1, with two parameters for a quadratic effect in the signal
#'           responses (5 parameters)}
#'  \item{4}{as 3, but with the same background level for both signals
#'           (4 parameters)}
#'  \item{5}{as 3, but with the same quadratic parameter for both signal
#'           responses (4 parameters)}
#'  \item{6}{as 5, but with the same background level for both signals
#'           (3 parameters)}
#'}
#'
#'@return A list; if an error occurs the only list component is
#'\describe{
#'  \item{message}{the error message}
#'}
#'If no error occurs the list has the following components:
#'\describe{
#'  \item{loglik}{the optimized log-likelihood}
#'  \item{npar}{the number of fitted parameters}
#'  \item{AIC}{Akaike's Information Criterion}
#'  \item{BIC}{Bayesian Information Criterion}
#'  \item{psi}{a list with components mu, sigma and p: mu and sigma each
#'             a vector of length ng with the means and standard deviations
#'             of the components of the fitted mixture model ON THE TRANSFORMED
#'             SCALE. p a matrix with one row per population and ng columns:
#'             the mixing proportions of the mixture components for each
#'             population}
#'  \item{post}{a matrix of ng columns and length(y) rows; each row r gives the
#'              ng probabilities that y[r] belongs to the ng components}
#'  \item{nobs}{the number of observations in y (excluding NA's)}
#'  \item{iter}{the number of iterations}
#'  \item{message}{an error message, "" if no error}
#'  \item{back}{a list with components mu.back and sigma.back: each a vector
#'              of length ng with the means and standard deviations of the
#'              mixture model components back-transformed to the original
#'              scale}
#'}
#'@examples
#' data(fitPoly_data)
#' mrkdat <- fitPoly_data$ploidy6$dat6x[fitPoly_data$ploidy6$dat6x$MarkerName == "mrk001",]
#'
#' # hexaploid, without specified populations
#' cdm <- CodomMarker(mrkdat$ratio, ng=7)
#' names(cdm)
#'
#' # hexaploid, with specified populations (4 F1 populations and a cultivar panel)
#' # first set the ptype for each population: p.F1 for F1 populations,
#' # p.HW for the panel, p.free for the F1 parents
#' ptype <- rep("p.HW", nrow(fitPoly_data$ploidy6$pop.parents))
#' ptype[!is.na(fitPoly_data$ploidy6$pop.parents[,1])] <- "p.F1"
#' ptype[unique(fitPoly_data$ploidy6$pop.parents)] <- "p.free" #all F1 parents
#' cdm <- CodomMarker(y=mrkdat$ratio, ng=7,
#'                    pop=fitPoly_data$ploidy6$pop,
#'                    pop.parents=fitPoly_data$ploidy6$pop.parents,
#'                    mutype=5, ptype=ptype)
#'
#'@export
CodomMarker <- function(y, ng,
                        pop.parents=matrix(c(NA,NA), nrow=1), #no population structure
                        pop=rep(1, length(y)), ##no population structure
                        mutype=0, sdtype="sd.const", ptype=NA,
                        clus=TRUE, mu.start=NA,
                        sd=rep(0.075, ng), p=NA,
                        maxiter=500, maxn.bin=200, nbin=200,
                        plothist=TRUE, nbreaks=40,
                        maintitle=NULL,
                        closeScreen=TRUE, fPinfo=NA) {

  res <- tryCatch( {
    if (is1NA(fPinfo)) {
      #if CodomMarker is called from fitOneMarker etc these checks are not needed
      #here we call checkCodomMarkerData to find problems that need an
      #error message,
      #and locally we adjust some variables if needed
      if (is.null(row.names(pop.parents)))
        row.names(pop.parents) <- 1:nrow(pop.parents)
      if (class(pop)[1] == "numeric") pop <- as.integer(round(pop))
      chkCMD <- checkCodomMarkerData(y, ng, pop, pop.parents, ptype)
      if (is.character(chkCMD)) stop(chkCMD) else {
        if (!is.null(chkCMD$pop)) pop <- chkCMD$pop
        if (!is.null(chkCMD$ptype)) ptype <- chkCMD$ptype
      }
      if (!is1NA(mu.start) &&
          (length(mu.start) != ng ||
           anyNA(mu.start) ||
           any(diff(mu.start) <= 0) ||
           mu.start[1] < 0 || mu.start[ng] > 1))
        stop("invalid mu.start")
      if (is1NA(sd)) sd <- rep(0.075, ng) else
        if (anyNA(sd) || any(sd <= 0) || !(length(sd) %in% c(1,ng)))
          stop("invalid sd") else
            if (length(sd) == 1) sd <- rep(sd, ng)
      #set p for each population in case it is NA or otherwise not matching to
      #number of populations, and make sure p sums to 1 for each population:
      if (is1NA(p)) {
        p <- matrix(1/ng, ncol=ng, nrow=nrow(pop.parents))
      } else {
        if (is.vector(p)) p <- matrix(p, nrow=1) #vector to matrix
        if (!is.matrix(p) || nrow(p) != nrow(pop.parents) || ncol(p) != ng ||
            anyNA(p) || any(p < 0) || any(p > 1) ||
            any(abs(1 - rowSums(p)) > 1e-6))
          stop("invalid p")
      }
      #p is now a matrix with nrow=nrow(pop.parents) (possibly 1) rows
      #and ng columns, each row summing to 1
      if ("p.F1" %in% ptype)
        segrArray <- getPolysomicSegr(ploidy=ng-1)
    } else {
      #called from fitOneMarker
      segrArray <- fPinfo$segrArray
    }

    # transform data to arcsine(square root):
    yw <- asin(sqrt(y))
    if (!is.na(mu.start[1])) mu.start <- asin(sqrt(mu.start))

    minyw <- min(yw, na.rm=TRUE)
    maxyw <- max(yw, na.rm=TRUE)

    #determine the number of free parameters of the model:
    if (mutype %in% 0:6) {
      mupar <- getMuModelNpar(mutype, ng)
    } else  stop("invalid mutype")

    sigmapar = switch (sdtype,
      sd.free = ng,
      sd.const = 1,
      sd.fixed = 0,
      stop("invalid sdtype")
    )

    ppar <- 0
    for (popi in 1:nrow(pop.parents)) {
      if (popi %in% pop.parents && ptype[popi] != "p.nodata") {
        # if population popi is F1 parent, then a single p parameter is assumed
        ppar <- ppar + 1
      } else {
        ppar <- ppar + switch (ptype[popi],
                               p.free = ng - 1,   # this may need modification
                               # if some classes do not contain data
                               p.HW = 1,
                               p.fixed = 0,
                               p.part = sum(p[popi,] > 0) - 1,
                               p.F1 = 0,
                               p.nodata = 0,
                               stop("invalid ptype") ) #end switch
      }
    }
    npar <- mupar + sigmapar + ppar
    result <- list(message="Error in CodomMarker") #default, to be replaced

    #check for ng==1 (1 peak):
    if (ng==1) {
      result <- fitOneDist(y=yw, sdtype=sdtype, sd.fixed=sd, npar=npar)
      #         population structure not relevant
    } else {
      if (is1NA(mu.start)) {
        # no starting values for mu were give, so calculate these first

        # initialize parameters through clustering if asked for
        if (clus) {
          yh <- yw[!is.na(yw)]
          # alternative (3-5-2012): use kmeans clustering, much faster
          # .Random.seed (affects and is affected by kmeans) is saved, set and
          # restored in fitOneMarker, but CodomMarker just uses the random
          # numbers as they come
          clus.init <- ClusterInit(yh, ng, closedips=TRUE) # returns mu's and sd on transformed scale
          mu.start <- clus.init$clus.mu
          sd.start <- clus.init$clus.sd #no need for rep, already ng values
          sd.start[sd.start < 0.01] <- 0.01  # don't allow sds smaller than 0.01
          if (sdtype == "sd.const") {
            sd.start <- rep(mean(sd.start), ng)
          } else if (sdtype == "sd.fixed") sd.start <- sd
        } else { # mu.start==NA and clus==FALSE; spread evenly over range
          ylo <- asin(sqrt(0.02))
          yhi <- asin(sqrt(0.98))
          mu.start <- seq(ylo, yhi, length.out=ng)
          sd.start <- sd
        }
      } else {
        # mu.start given, use also given or builtin sd.start:
        sd.start <- sd
      }

      result <- EMGaussMix(y=yw, ng=ng, pop=pop, ptype=ptype, mutype=mutype,
                           sdtype=sdtype, sd.fixed=sd,
                           p.start=p, mu.start=mu.start, sd.start=sd.start,
                           npar=npar, maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
                           pop.parents=pop.parents, segrArray=segrArray)
    } # if ng==1 else
    if ("message" %in% names(result)) {
      if (result$message == "") {
        # transform mu and sigma back to original scale -
        # note that the naive y <- sin(x)^2 is not correct here!
        #mu.back <- sin(result$psi$mu)^2 + 0.5*(2*(cos(result$psi$mu)^2 - sin(result$psi$mu)^2))*(result$psi$sigma^2)   # Use E(f(y)) ~ f(Ey) + 0.5*f''(Ey) var(y)
        # Use E(f(y)) ~ f(Ey) + 0.5*f''(Ey) var(y) :
        mu.back <- sin(result$psi$mu)^2 +
          (cos(result$psi$mu)^2 - sin(result$psi$mu)^2) * result$psi$sigma^2
        # Use sd(f(y)) ~ f'(Ey) sd(y) :
        sigma.back <- 2*sin(result$psi$mu)*cos(result$psi$mu)*result$psi$sigma
        back <- list(mu.back=mu.back, sigma.back=sigma.back)
        result$back <- back

        # prepare for histogram visualization
        if (plothist) {
          tryCatch({
            if (closeScreen) close.screen(all.screens=TRUE)
            drawCompositePlot(rr=y, rrpop=pop, drr=numeric(0),
                              pop.parents=pop.parents,
                              nonparents=setdiff(1:nrow(pop.parents),
                                                 pop.parents),
                              #geno=rep(NA, length(y)),
                              psi=result$psi,
                              maintitle=maintitle, nbreaks=40,
                              drawGeno=FALSE)
            if (closeScreen) close.screen(all.screens=TRUE)
          }, error = function(ex) {
            warning(paste("Plotting error in CodomMarker:", ex$message, "\n"),
                    call.=FALSE)
          })
        }
      } else { #result$message!=""
        result$message <- paste("Error1 in CodomMarker: ",
                                result$message, sep="")
      }
      #end of "message" in names(result)
    } else {
      result$message <- "Error3 in CodomMarker"
    }
    result
  }, error = function(ex) {
    list(message=paste("Error2 in CodomMarker: ", ex$message, sep=""))
  } )
  res
} # CodomMarker

# ClusterInit ******************************
# Initialisation of parameters by clustering
# ******************************************
# parameters:
# yh: ratios on transformed arcsine(sqrt) scale; may not contain NA
# ng: total number of possible clusters (one more than ploidy level)
# closedips: if TRUE any clusters with less than (1/ng/4) of observations,
#            between larger clusters, are discarded
# return value: list with 2 elements:
# $clus.mu: vector of length ng with cluster means on transformed scale
# $clus.sd: vector of length ng, each number is the same within-cluster standard deviation
#           on transformed scale

ClusterInit <- function(yh, ng, closedips=TRUE) {
  tryCatch ({
    seq2 <- function(low,up,n) seq(low,up, length.out=n+2)[-c(1,n+2)]  # help function for filling gaps

    # Try 1 to ng clusters, decide on number of clusters, fill in gaps at ends
    km <- list(); wss <- numeric(ng); BIC <- numeric(ng)
    for (i in 1:ng){
      km[[i]] <- kmeans(yh, centers=i, nstart=10*i)  # try different starts, increasing numbers for more clusters
      wss[i] <- km[[i]]$tot.withinss

      clus.mu <- km[[i]]$centers             #i*1 matrix of cluster centers (transformed scale)
      clus.sd <- sqrt(wss[i]/(length(yh)-i)) #one numeric value: overall within-cluster sd (transformed scale)
      clus.p <- km[[i]]$size/length(yh)      #vector of fraction of observations in each cluster

      if (i == 1) loglik <- sum(log(sapply(yh, dnorm, clus.mu, clus.sd)))
      else loglik <- sum(log(t(sapply(yh, dnorm, clus.mu, clus.sd)) %*% clus.p))
      BIC[i] <- -2*loglik + log(length(yh))*(7*i)   # use more heavy penalty than ordinary BIC to avoid spurious clusters
    }

    nrclust <- which.min(BIC)

    clus.mu <- km[[nrclust]]$centers #nrclust*1 matrix
    clus.sd <- rep(sqrt(wss[nrclust]/(length(yh)-nrclust)), ng)
    o <- order(clus.mu)
    clus.mu <- clus.mu[o] #vector of length nrclust
    if (closedips) {
      dip_threshold <- 1/ng/4 #clusters with less than this fraction of samples are discarded
      clus.p <- km[[nrclust]]$size/length(yh)
      clus.p <- clus.p[o]
      keep <- which(clus.p >= dip_threshold) #at least one
      #we delete any small clusters within this range (not at the outside,
      #and they are larger than dip_threshold):
      if (max(keep) - min(keep) >= length(keep)) {
        firstTopAt <- 2
        while (firstTopAt <= length(keep) &&
                 clus.p[keep[firstTopAt]] >= clus.p[keep[firstTopAt-1]])
          firstTopAt <- firstTopAt + 1
        firstTopAt <- firstTopAt - 1
        lastTopAt <- length(keep) - 1
        while (lastTopAt >= 1 &&
                 clus.p[keep[lastTopAt]] >= clus.p[keep[lastTopAt+1]])
          lastTopAt <- lastTopAt - 1
        lastTopAt <- lastTopAt + 1
        if (lastTopAt > firstTopAt +1) {
          gap_threshold <- min(clus.p[keep[c(firstTopAt, lastTopAt)]]) / 4
          del <- which(clus.p[min(keep):max(keep)] < gap_threshold) + min(keep) - 1
          if (length(del) > 0) clus.mu <- clus.mu[keep[-del]]
        }
      }
    }
    extraclust <- ng - length(clus.mu)

    if (extraclust > 0) {
      #The remaining peaks are considered one consecutive block;
      #we add additional peaks (up to ng) at the outsides, based on the
      #distances between the first mu and 0 vs last mu and pi/2:
      maxy <- pi/2
      lo.mu <- min(clus.mu)
      hi.mu <- max(clus.mu)
      lo.empty <- lo.mu
      hi.empty <- pi/2-hi.mu
      frac.lo <- lo.empty/(lo.empty+hi.empty)
      frac.hi <- 1-frac.lo

      n.lo <- round(frac.lo * extraclust)
      n.hi <- round(frac.hi * extraclust)
      if ((n.lo+n.hi) > extraclust) {
        # happens only if n.lo and n.hi both rounded up, i.e. exact n in both cases x.5
        # decreasing the smallest one tends to shift distribution to one side,
        # decreasing the largest tends to center the distribution
        # we go for the shift:
        if (n.lo < n.hi) n.lo <- n.lo - 1 else n.hi <- n.hi - 1
      } else {
        if ((n.lo+n.hi) < extraclust) {
          # happens only if n.lo and n.hi both rounded down, i.e. exact n in both cases x.5
          # increasing the largest one tends to shift distribution to one side,
          # increasing the smallest tends to center the distribution
          # we go for the shift:
          if (n.hi > n.lo) n.hi <- n.hi + 1 else n.lo <- n.lo + 1
        }
      }
      clus.mu <- c(seq2(0, lo.mu, n.lo), clus.mu, seq2(hi.mu, maxy, n.hi))
    }

    list(clus.mu=clus.mu, clus.sd=clus.sd) #both on transformed scale
  }, error = function(x) {
    #if clustering error (because less than ng distinct ratio values,
    #as when (almost) all ratios are exactly 0 or 1)
    #do as when mu.start=NA and clus=FALSE, but now with extreme values:
    list (clus.mu=seq(0.02, (pi/2) - 0.02,length.out=ng),
          clus.sd <- rep(0.075,ng))
  })
} #ClusterInit

# ***********************
# Supporting functions
# ***********************


# PlotHistDensity **************************
# The main function for plotting the results

dnormasr <- function(x, mu, sigma) { dnorm(asin(sqrt(x)), mean=mu, sd=sigma) }
#dnormlogit <- function(x, mu, sigma) { dnorm(log(x/(1-x)), mean=mu, sd=sigma) }

PlotHistDensity <- function(rr, drr=numeric(0),
                            rrP1=numeric(0), rrP2=numeric(0),
                            mu, sigma, p,
                            breaks,
                            #trafo="asr",
                            maintitle="", xaxis="n",  xlabel=NULL,
                            nbreaks=40, fitcol="darkgreen") {
  #This function draws a plot for one population and one fitted marker.
  #The plot consists of a histogram of the (polyploid) F1 population or
  #association panel, with superimposed a histogram of diploid samples. If
  #the population is an F1 population the positions of the parental samples
  #are indicated. Over the histograms the distribution of the fitted mixture
  #model is drawn.
  #rr: vector of signal ratios (asin(sqrt(x)) transformed) of all polyploid
  #    samples for one marker for the F1 population or association panel.
  #drr: vector of signal ratios (asin(sqrt(x)) transformed) of all diploid
  #     samples for one marker
  #rrP1, rrP2: signal ratios for the samples of parent 1 and parent 2
  #mu: vector of length <ploidy+1> with means of the mixture components
  #    (on transformed scale)
  #    In case mu has not length 5 or the first element of mu is NA
  #    the model is not drawn
  #sigma: vector of length 5 with stdev of the mixture components
  #      (on transformed scale)
  #p: vector of length 5 with the mixing proportions of the 5 mixture components
  #   for the population
  #breaks: vector of breaks for hist on the transformed ratio scale
  #maintitle: the title to print in the plot
  #xaxis: allows to draw ("s") or suppress ("n") the x-axis
  #xlabel: NULL or x-axis title

  htet <- hist(rr, breaks=breaks, plot=FALSE)
  if (length(drr) > 0) {
    hdip <- hist(drr, breaks=breaks, plot=FALSE)
  } else hdip <- list(counts=0)
  ylabel <- "frequency"
  maxy <- max(htet$counts)
  maxdip <- max(hdip$counts)
  if (maxy <= 0) {
    #for populations with no non-missing data
    if (maxdip == 0) maxy <- 1 else maxy <- maxdip
  } else {
    #scale diploid counts so the maximum bar is half the max polyploid bar:
    hdip$counts <- hdip$counts * maxy / maxdip / 2
  }
  xlim <- c(0, 1)
  htet$ylim <- c(0, 1.05 * maxy) #so ylim will be part of the return value
  barplot (htet$counts,
           width=(htet$breaks[length(htet$breaks)] - htet$breaks[1]) /
             (length(htet$breaks) - 1),
           space=0, xlim=xlim, ylim=htet$ylim, col="white",
           main=paste("\n\n",maintitle,sep=""),
           xaxt=xaxis, xlab=xlabel, ylab=ylabel, yaxt="n")
  if (xaxis=="s") {
    #plot an x-axis and an y-axis with "0"
    axis(side=1, labels=TRUE, at=seq(0, 1, by=0.2))
    axis(side=2, labels=TRUE)
  } else {
    #plot only an y-axis, with the "0" at the lowest tick omitted:
    yticks <- axTicks(side=2)
    ylabs <- as.character(yticks) ; ylabs[1] <- ""
    axis(side=2, at=yticks, labels=ylabs)
  }
  #area <- 1
  classwidth <- breaks[2] - breaks[1]
  area <- length(rr[!is.na(rr)]) * classwidth

  if (length(drr) > 0) {
    #plot diploids histogram into tetraploids plot:
    nbreaks <- length(breaks) - 1
    binwidth <- 1 / nbreaks
    barwidth <- 1/3 #relative to binwidth
    par(new=TRUE) # start new plot in existing barplot
    barplot (hdip$counts,
             width=binwidth * barwidth,
             space=1/barwidth - 1,
             xlim=c(binwidth * (1-barwidth)/2, 1 + binwidth * (1 - barwidth)/2),
             ylim=htet$ylim, col="gray", main="", xlab="", xaxt="n", yaxt="n")
  }
  oldcex <- par("cex")
  if (is.null(mu)) ng <- 0 else ng <- length(mu)
  if (ng > 0 && length(sigma) == ng && length(p) == ng &&
      !anyNA(c(mu, sigma, p))) {
    #model fitted, draw model:
    #draw the model density lines:
    r <- seq(from=xlim[1], to=xlim[2], length.out=301) #by=((xlim[2]-xlim[1])/300))
    for (i in 1:ng) {
      area.i <-
        integrate(dnormasr, lower=0.001, upper=0.999 , mu=mu[i],
                  sigma=sigma[i], subdivisions=300)$value
      #the area is changed due to the y=asin(sqrt(x)) transformation,
      #correct for this
      lines(cbind(r, (area / area.i) * p[i] *
                     dnorm(asin(sqrt(r)), mean=mu[i], sd=sigma[i])),
            lwd=1, col=fitcol)
    }
    #plot the means (at back-transformed positions):
    muback <- sin(mu)
    muback <- muback * muback
    #points(rep(0, ng) ~ muback, pch=17, cex=3*oldcex, col=fitcol) #the means
    text(x=muback, y=rep(0, ng), labels=as.character(0:(ng-1)),
         cex=2*oldcex, adj=c(0.5, 0), col=fitcol) #the means, as numbers with
    #                                              bottom on X axis
  }
  #plot the parental samples:
  if (length(rrP1) > 0) points(rrP1, rep(0,length(rrP1)), pch=17,
                               cex=3*oldcex, col="red")
  if (length(rrP2) > 0) points(rrP2, rep(0,length(rrP2)), pch=17,
                               cex=3*oldcex, col="blue")
  par(cex=oldcex)
  box("plot")
} # PlotHistDensity

getResultnames <- function(ng, pop.parents) {
  #these are the names of the modeldata component of the return value of fitTetra
  result <-
    c("marker", "MarkerName", "m", "model",
      "nsamp", "nsel", "npar", "iter", "dip",
      "LL", "AIC", "BIC",
      "selcrit", "minsepar", "meanP", "P80", "P90", "P95", "P975", "P99",
      paste(rep("muact", ng), 0:(ng-1), sep=""), #actual means of the samples in
      #                              each peak on arcsine-sqrt transformed scale
      paste(rep("sdact", ng), 0:(ng-1), sep=""), #actual sd's of the samples in
      #                              each peak on arcsine-sqrt transformed scale
      paste(rep("Pact", ng), 0:(ng-1), sep="")) #actual freqs of the samples in each peak
  #   Note that Pact is over all populations together, NOT per population
  #   like the fitted P's (see below). In fitOneMarker we make use of the combined Pact!
  result <- c(result,
              paste(rep("mutrans", ng), 0:(ng-1), sep=""), #mu's on arcsine-sqrt transformed scale
              paste(rep("sdtrans", ng), 0:(ng-1), sep="")) #sd's on arcsine-sqrt transformed scale
  if (nrow(pop.parents) == 1) {
    result <-c(result,
               paste(rep("P", ng), 0:(ng-1), sep="")) #mixture proportions
  } else for (r in 1:nrow(pop.parents)) {
    result <-
      c(result,
        paste(rep(paste("P", rownames(pop.parents)[r], sep="_"), ng),
              0:(ng-1), sep="_")) #mixture proportions for each pop
  }
  c(result,
    paste(rep("mu", ng), 0:(ng-1), sep=""),  #backtransformed mu's
    paste(rep("sd", ng), 0:(ng-1), sep=""), #backtransformed sd's
    "message")
} #getResultnames


# calcSelcrit calculates the selcrit (selection criterion based on the values
# of the BIC and sdtrans0 components of allmodeldata and on sd.target
# include is a boolean vector of length nrow(allmodeldata); where include=FALSE
# the allmodeldata$selcrit becomes NA
calcSelcrit <- function(allmodeldata, include, sd.target) {
  #different 20161116: we now scale according to the observed BIC range:
  #for sd <= sd.target: selcrit == BIC
  #for sd = 2* sd.target: selcrit <- BIC + BICrange + 5
  allmodeldata$selcrit <- rep(NA, nrow(allmodeldata))
  include[is.na(include)] <- FALSE
  allmodeldata$selcrit[include] <- allmodeldata$BIC[include]
  if (!is1NA(sd.target) && any(!is.na(allmodeldata$BIC[include]))) {
    BICrange <- max(allmodeldata$BIC[include], na.rm=TRUE) -
      min(allmodeldata$BIC[include], na.rm=TRUE)
    sdvalues <- allmodeldata$sdtrans0 #all sdtrans are equal, take the first
    mulf <- rep(1, nrow(allmodeldata))
    changerows <- include & !is.na(sdvalues) & sdvalues > sd.target
    mulf[changerows] <- sdvalues[changerows] / sd.target #larger sdvalue -> larger mulf
    corr <- (mulf - 1) * (BICrange + 5)
    allmodeldata$selcrit[changerows] <- (allmodeldata$BIC + corr)[changerows]
  }
  allmodeldata
} #calcSelcrit

#MarkerResult collects the output of CodomMarker
#for comparisons of models, and performs some initial calculations
#for genotype allocation
#the return value is a list with
#stats: a data frame, and
#probs: itself a list of data frames
#MarkerResult adds the results of the current marker to row [model] of stats
#and adds element[[model]] of prob
#the whole list of model results to which these data are added is passed
#as parameter mrkresult
# the model determines how the means and mixing proportions of the mixture
# component distributions are modelled

MarkerResult <- function(marker, markername, ratio,
                         model, origPopstruct, mrkPopstruct,
                         mutype,
                         ng, mrkresult, nsamp,
                         mustart, sd=rep(0.075, ng),
                         maxiter, maxn.bin, nbin,
                         plothist,
                         plot.type, plot.dir,
                         modelcount) {
  #model: model number, determines how the means and mixing proportions of the
  #       mixture component distributions are modelled
  #modelcount: 4 or 8, the number of models fitted with a given mustart,
  #            used only for plotting (if plothist TRUE)
  model.ps    <- getModelPopstruct(model, origPopstruct, mrkPopstruct, ng)
  modelresult <- mrkresult$stats
  probslist <- mrkresult$probs
  Pactlist    <- mrkresult$Pact #matrices with the actual fractions of samples
  #                              in each peak, for each population
  modelname <- getModelName(mutype, model.ps, ng)
  result <- tryCatch({
    CodomMarker(y=ratio, ng=ng,
                pop=model.ps$population,
                pop.parents=model.ps$pop.parents,
                mutype=mutype, sdtype="sd.const",
                ptype=model.ps$ptype,
                clus=FALSE, #as we always supply a mu.start
                mu.start=mustart, sd=sd, p=model.ps$p.start,
                maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
                plothist=FALSE, # not plothist!
                #plot.type=plot.type, plot.dir=plot.dir,
                maintitle=paste(marker, markername, modelname),
                closeScreen=FALSE, fPinfo=list(segrArray=model.ps$segrArray))
    #DONE plothist=FALSE in call to CodomMarker,
    #TODO: check separate plotting with popstruct if plothist TRUE
  }, error = function(x) {x} )
  modelresult$marker[model] <- marker
  modelresult$MarkerName[model] <- markername
  modelresult$m[model] <- model
  modelresult$model[model] <- modelname
  modelresult$nsamp[model] <- nsamp #including the ratio=NA but
  #                                  excluding the select=FALSE samples
  modelresult$nsel[model] <- length(ratio)
  if (class(result)[1] == "character") {
    modelresult$message[model] <- result[1]
  } else {
    modelresult$message[model] <- gsub("[\r\n]", " ", result$message) #removes
    #  any newlines in message: message "$ operator is invalid for atomic vectors"
    #  seems to contain a newline
  }

  if (plothist) {
    startPlotall(model, origPopstruct$orig_nonparents,
                  plot.type, plot.dir, marker, markername)
  }

  if (class(result)[1] == "character" || result$message != "") {
    #CodomMarker error
    modelresult$npar[model] <- NA
    modelresult$iter[model] <- NA
    modelresult$dip[model] <- NA
    modelresult$LL[model] <- NA
    modelresult$AIC[model] <- NA
    modelresult$BIC[model] <- NA
    if (plothist) {
      errorplot(paste(marker, markername, modelname))
    }
  } else {
    #CodomMarker successful
    modelresult$npar[model] <- result$npar
    modelresult$iter[model] <- result$iter
    modelresult$LL[model] <- result$loglik
    modelresult$AIC[model] <- result$AIC # Akaike's Information Criterion,
    #                                      better than BIC for small models
    modelresult$BIC[model] <- result$BIC # Bayesian Information Criterion,
    #                                 optimal for complex models
    #                                 more severe penalty than AIC if nsel>=8
    #                                 we use BIC for model selection
    #the following statistics may also reflect the quality of the fit
    #but are not used any more:
    mudiff <- diff(result$psi$mu)
    sigavs <- (result$psi$sigma[2:ng] + result$psi$sigma[1:(ng-1)]) / 2
    #         the average of successive sigma's
    separation <- mudiff/sigavs
    #             a measure of the (ng-1) separations between the peaks
    modelresult$minsepar[model] <- min(separation)
    #the following calculates the fraction of genotyped samples at
    #various P thresholds. These fractions relate to the total number of
    #ratio values passed to CodomMarker; this includes only the non-NA samples
    #for which select is TRUE (in fitOneMarker)
    probs <- data.frame(result$post)
    names(probs) <- paste("P", 0:(ng-1), sep="") # "P0" .. "P<ploidy>"
    #list for each sample the maximum p value and the maxgeno corresponding to it:
    probs$maxgeno <- max.col(result$post)  # for every row (sample) the max peak
    #                                            (0..(ng-1) instead of 1..ng)
    probs$maxP <- result$post[cbind(1:nrow(result$post), probs$maxgeno)]
    #             for every row the P at the max
    if (anyNA(probs$maxP)) {
      warning(paste("marker", marker, "model", model,"NAs in maxP")) #TODO: remove debug check
      #stop(paste("marker", marker, "model", model,"error in maxP")) #TODO: remove debug check
    }
    probs$maxgeno <- probs$maxgeno - 1 # 0..ploidy instead of 1..ng
    maxP <- probs$maxP # saves code
    modelresult$meanP[model] <- mean(maxP) #the average maxP over all samples
    #calc fractions samples with maxP>0.80 .. 0.99
    modelresult$P80[model] <- length(maxP[maxP>0.80])/length(maxP)
    modelresult$P90[model] <- length(maxP[maxP>0.90])/length(maxP)
    modelresult$P95[model] <- length(maxP[maxP>0.95])/length(maxP)
    modelresult$P975[model] <- length(maxP[maxP>0.975])/length(maxP)
    modelresult$P99[model] <- length(maxP[maxP>0.99])/length(maxP)
    #the following give the parameters of the fitted model and are used
    #by fitOneMarker for plotting
    mutrans0 <- which(names(modelresult)=="mutrans0")
    modelresult[model, mutrans0:(mutrans0 + ng - 1)] <- result$psi$mu
      #was as.numeric(result$psi$mu)
    modelresult[model, (mutrans0 + ng):(mutrans0 + 2*ng - 1)] <-
      result$psi$sigma # was as.numeric(result$psi$sigma)
    #psi$p is now a matrix with one row per current population.
    if (model.ps$allcombined) {
      #psi$p has only one row and we put this in the slots for the first
      #population in modelresult
      modelresult[model, (mutrans0 + 2*ng):(mutrans0 + 3*ng - 1)] <-
        result$psi$p[1,]
    } else {
      #psi$p has rows for all original populations, and we
      #put all in one go into modelresult
      colstochange <-
        mutrans0 + ((2*ng):((2+nrow(origPopstruct$orig_pop.parents))*ng - 1))
      modelresult[model, colstochange] <- t(result$psi$p)
        #was as.numeric(t(result$psi$p))
      #the p for all original populations
    }
    #next the mu's and sigma's backtransformed to 0..1 scale:
    mu0 <- which(names(modelresult)=="mu0")
    modelresult[model, mu0:(mu0 + ng - 1)] <- result$back$mu
    modelresult[model, (mu0 + ng):(mu0 + 2*ng - 1)] <- result$back$sigma
    #after the general model data we process the data per sample,
    #and obtain the actual mus, sds, ps, and check for dips:
    probslist[[model]] <- probs
    Pact <- do.call(rbind, tapply(probs$maxgeno+1, INDEX=model.ps$population,
                                  FUN=tabulate, nbins=ng))
    #       Pact is actual fraction of samples in each peak, per population.
    #       Because seldat is ordered by SampleName in fitOneMarker, the sample
    #       order in probs is the same as in seldat and in model.ps$population
    colnames(Pact) <- 0:(ng-1)
    totPact <- colSums(Pact)
    Pact <- Pact / rowSums(Pact) # may produce NaN's, are detected with is.na()
    totPact <- totPact / sum(totPact)
    #print(paste("model=", model," totPact=",paste(totPact, collapse=" ")))
    Pactlist[[model]] <- Pact
    #check if there is a dip (smaller peak surrounded by larger peaks);
    #in any of the non-parent populations:
    dip <- length(find.dips.populations(Pact, pops=model.ps$nonparents)) > 0
    modelresult$dip[model] <- as.integer(dip)
    # calculate the actual mean, sd and P of the samples in each peak
    # Even with multiple populations we give actual means or counts OVER ALL
    # populations: in fitOneMarker we make use of the combined Pact!
    asr <- asin(sqrt(ratio))
    muact <- tapply(asr, INDEX=probs$maxgeno, FUN=mean) #mean per peak
    sdact <- tapply(asr, INDEX=probs$maxgeno, FUN=sd) #sd per peak
    nam <- as.numeric(names(muact)) #some classes may be missing
    muact0 <- which(names(modelresult) == "muact0")
    for (n in 0:(ng-1)) { #muact and sdact may have missing categories
      w<-which(nam == n)
      if (length(w) == 1) {
        modelresult[model, muact0 + n] <- muact[w]
        modelresult[model, muact0 + ng + n] <- sdact[w]
      } else {
        modelresult[[muact0+n]][model] <- modelresult[[mutrans0+n]][model]
        #                                 if no samples in peak use model mean
      }
    }
    modelresult[model, (muact0 + 2*ng):(muact0 + 3*ng - 1)] <- totPact
    if (plothist) {
      drawCompositePlot(rr=ratio, rrpop=model.ps$population, drr=numeric(0),
                        pop.parents=model.ps$pop.parents,
                        nonparents=model.ps$nonparents,
                        #geno=rep(NA, length(y)),
                        psi=result$psi,
                        maintitle=paste(marker, " ", markername, " ",
                                        model, ": ", modelname, sep=""),
                        nbreaks=40, drawGeno=FALSE)
    } #plothist
  } #CodomMarker succesful
  if (plothist) {
    savePlotall(model=model, orig_nonparents=origPopstruct$orig_nonparents,
                modelcount=modelcount)
  }
  list(stats=modelresult, probs=probslist, Pact=Pactlist)
} # MarkerResult

lengthNoNA <- function(x) { sum(!is.na(x)) }

#padded returns a string representing positive integer x left-padded with 0's
#to the length of integer maxx or of max(x)
padded <- function(x, maxx=0) {
  #paste(substr("00000",1,nchar(as.character(maxx))-nchar(as.character(x))),x,sep="")
  formatC(x, width=nchar(max(x, maxx, na.rm=TRUE)), flag="0")
}

#errorplot plots an empty plot with an error message
errorplot <- function(main, message="not fitted") {
  #print(screen())
  oldpar <- par(c("cex", "mar", "mex"))
  nwcex <- 0.66 * oldpar$cex
  nwmex <- oldpar$mex
  plot(0, main=main, xlim=c(0,1), ylim=c(0,1),
       xaxt="n", yaxt="n", xlab="", ylab="", pch=NA) #empty space with title
  text(message, x=0.5, y=0.5)
  #box("figure", col="orange")
  par(oldpar)
}

drawRejectedPlot <-  function(rr, rrpop, drr=numeric(0), ploidy, pop.parents,
                              nonparents, maintitle) {
  drawCompositePlot(rr, rrpop, drr, pop.parents, nonparents,
                    #geno=rep(NA, length(rr)),
                    psi=list(mu=rep(NA, ploidy+1)),
                    maintitle=maintitle, nbreaks=40, drawGeno=TRUE)
}

drawFittedPlot <- function(rr, rrpop, drr=numeric(0), pop.parents,
                           nonparents, geno=NULL, psi, maintitle) {
  drawCompositePlot(rr, rrpop, drr, pop.parents, nonparents,
                    geno=geno,
                    psi=psi,
                    maintitle=maintitle, nbreaks=40, drawGeno=TRUE)
}

drawCompositePlot <- function(rr, rrpop, drr=numeric(0), pop.parents, #parents,
                              nonparents,
                              geno=NULL, psi,
                              maintitle,
                              nbreaks=40, drawGeno=FALSE) {
  #This function draws a composite plot for one marker.
  #The plot consists of one histogram per F1-population or association panel,
  #plus one rug plot showing the assigned genotypes for all samples if
  #drawGeno is TRUE (even if geno is NULL!)
  #rr: vector of signal ratios (asin(sqrt(x)) transformed) of all polyploid
  #    samples for one marker
  #rrpop: integer of same length and same order as rr; values indicate rows of pop.parents
  #drr: vector of signal ratios (asin(sqrt(x)) transformed) of all diploid
  #     samples for one marker
  #pop.parents: matrix with one row per population indicated by rrpop values:
  #           nrow(pop.parents) == length(levels(rrpop))
  #           each row has 2 integers: both 0 for panels and parentals,
  #           two other row numbers for F1 populations
  #           row names are used for histogram titles
  #nonparents: integer vector with the rows in pop.parents that are not a parent
  #            (i.e. for which a histogram must be drawn)
  #   DONE?: use the popstruct info to determine the number of panels in the plot
  #   (based on orig_pop.parents)
  #   deal with psi$p being NA for populations and parents without data
  #   Note that populations in orig_pop.parents for which no data are available
  #   must show an empty plot with the name of the population!
  #geno: integer vector of same length and same order as rr, with the assigned
  #      genotypes
  #      (0..ploidy or NA)
  #psi: a list with components
  # $mu: vector of length ng with means of the mixture components
  #      (on transformed scale)
  # $sigma: vector of length ng with stdev of the mixture components
  #      (on transformed scale)
  # $p: matrix with ng columns and same number of rows as pop.parents.
  #     each row has the mixing proportions of the ng mixture components for one
  #     population
  # In case the first element of mu is NA the model is not drawn
  #maintitle: the base title; for each histogram the population name is appended
  #nbreaks: the number of breaks for the tetraploid histogram
  #drawGeno: if TRUE (default) a rug plot with the genotypes of all samples
  #          is drawn
  nhists <- length(nonparents)
  #plotheights: the histograms together take up 4/5 of the height,
  #the rug plot 1/5
  # new way, using split.screen instead of layout, as this needs to draw
  # a composite plot into one of the subscreens when called from fitTetra
  # with plotall TRUE:
  # Currently the subscreen of the whole page is assumed to be the active screen
  oldpar <- par(c("cex", "mar", "mex"))
  nwcex <- 0.66 * oldpar$cex
  nwmex <- oldpar$mex
  figs <- matrix(nrow=nhists + 2, ncol=4)
  #for each screen in the split screen one row; the 4 columns indicate the
  #boundaries (left-right-bottom-top) whre 0=bottom/left and 1=top/right of page
  figs[1,] <- c(0, 1, 0.95, 1) #title pane
  if (drawGeno) figs[nhists + 2,] <- c(0, 1, 0, 0.2) else #rug plot
    figs[nhists + 2,] <- c(0, 1, 0, 0.05) #empty, for X axis
  hh <- (figs[1, 3] - figs[nhists+2, 4]) / nhists #height of each hist plot
  for (h in 1:nhists) figs[h + 1,] <- c(0, 1, 0.95 - h*hh, 0.95 - (h-1)*hh)
  #print(paste("DrawCompositePlot: all screens were ", paste(split.screen(), collapse=" ")))
  #print(paste("DrawCompositePlot: current screen is", screen()))
  scr <- split.screen(figs)
  #print(paste("DrawCompositePlot: new screens are ",paste(scr, collapse=" ")))
  #print(paste("DrawCompositePlot: all screens are ",paste(split.screen(), collapse=" ")))

  #draw title:
  screen(scr[1])
  par(mar=c(0,4,1,1)) #margins within screen, bottom-left-top-right
  par(cex=nwcex, mex=nwmex)
  #print(paste("screens=",paste(split.screen(),collapse=" "),
  #            " cex=",par("cex"), " mex=",par("mex"),sep=""))
  plot(1, type="n", bty="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA,
       main=maintitle) #draws only maintitle above empty space for plots
  #box("figure")
  #draw histograms & models:
  suppressWarnings({
    minrr <- min(c(rr, drr), na.rm=TRUE)
    maxrr <- max(c(rr, drr), na.rm=TRUE)
  })
  if (is.na(minrr) || is.infinite(minrr) || maxrr==minrr) {
    nbars <- nbreaks + 1
  } else nbars <- ceiling(nbreaks/(maxrr - minrr))
  breaks <- seq(0, 1, by=1/nbars)
  xaxt <- "n"; xlab <- NULL
  for (his in seq_along(nonparents)) {
    pop <- nonparents[his]
    screen(scr[his+1])
    par(mar=c(0,4,0,1))
    par(cex=nwcex, mex=nwmex)
    if (length(nonparents) == 1) main <- NULL else main <- rownames(pop.parents)[pop]
    if (is.na(pop.parents[pop, 1])) {
      rrP1 <- rrP2 <- numeric(0)
    } else {
      rrP1 <- rr[rrpop == pop.parents[pop, 1]]
      rrP2 <- rr[rrpop == pop.parents[pop, 2]]
    }
    rrmain <- rr[rrpop == pop]
    if (!drawGeno && his == nhists) {
      xaxt <- "s"; xlab <- "signal ratio"
    }
    if (his > 1) drr <- numeric(0) #only plot diploids in first histogram
    PlotHistDensity(rr=rrmain,
                    drr=drr,
                    rrP1=rrP1,
                    rrP2=rrP2,
                    #geno=geno[rrpop == his],
                    mu=psi$mu,
                    sigma=psi$sigma,
                    p=psi$p[pop,],
                    breaks=breaks,
                    maintitle=main,
                    xaxis=xaxt, xlabel=xlab)
    #box("figure", col=c("red","blue")[his %% 2 + 1])
  }
  if (drawGeno) {
    #lower plot: the assigned genotypes
    #NOTE that all samples in rr are plotted, including the parentals and also
    #     any samples with an rrpop not in pop.parents
    screen(scr[length(scr)])
    par(mar=c(5,4,0,1))
    par(cex=nwcex, mex=nwmex)
    #box("figure", col="green")
    if (is.null(geno)) geno <- rep(-2, length(rr)) #all missing
    geno[is.na(geno)] <- -2 #plot well separated from the valid scores 0:4
    sampcol <- rep("blue",length(rr))
    sampcol[geno < 0] <- "red"
    #draw the empty box with axes for the rug plots:
    plot(0, type="n",
         ylab="geno", ylim=c(-3, length(psi$mu)),
         xlab="signal ratio", xlim= c(0, 1) )
    for (g in c(-2, 0:(length(psi$mu)-1))) {
      lines(c(g, g) ~ c(par("usr")[1]+0.004, par("usr")[2] - 0.004),
            col="lightgrey", lty=1) #dotted lines
    }
    points(geno ~ rr, pch="|", cex=0.5, col=sampcol)
  }
  par(oldpar) #needed if CodomMarker called directly by user
} #drawCompositePlot

checkPlotType <- function(plot.type) {
  #returns a character vector: [1] is the checked plot.type, [2] a message or ""
  result <- c("none", "")
  plot.type <- tolower(plot.type[1]) #ignore all but the first element
  if (plot.type != "none") {
    if (!requireNamespace("grDevices", quietly=TRUE)) {
      result[2] <- "package grDevices not available, no plots generated"
    } else {
      cap <- capabilities()
      #first check if we can generate a requested non-png type:
      if (plot.type == "emf") {
        if (.Platform$OS.type == "windows" &&
            requireNamespace("devEMF", quietly=TRUE) && exists("emf")) {
          result[1] <- "emf"
        } else {
            result[1] <- "png"
            result[2] <- "package devEMF (emf graphics) not available"
        }
      } else if (plot.type == "svg") {
        if (cap["cairo"] && exists("svg")) result[1] <- "svg" else {
          result[1] <- "png"
          result[2] <- "cairo (svg) graphics not available"
        }
      } else if (plot.type == "pdf") {
        if (exists("pdf")) result[1] <- "pdf" else {
          result[1] <- "png"
          result[2] <- "pdf graphics not available"
        }
      } else if (plot.type != "png") {
        result[1] <- "png"
        result[2] <- paste("plot type", plot.type, "not supported")
      }
      if (!cap["png"] || !exists("png")) {
        if (plot.type == "png") {
          result[1] <-"none"
          result[2] <- "png graphics not available, no plots produced"
        } else {
          #one of the other plot types requested
          if (result[1] == "png") {
            #the other type was not possible, was changed to png
            result[1] <- "none"
            result[2] <- paste(result[2], ", no plots produced", sep="")
          }
        }
      } else {
        #png graphics available
        if (plot.type == "png") result[1] <- "png" else {
          if (result[1] == "png") {
            #another type was not possible, was changed to png
            result[2] <- paste(result[2], ", png plots produced instead", sep="")
          }
        }
      }
    } #grDevices available
  } #plot.type != "none"
  result
} #checkPlotType

checkPlot <- function(plot, plot.type, plot.dir, filePrefix) {
  #for use in fitMarkers and fitOneMarker
  #plot, plot.type and plot.dir are parameters of one of those
  #logfile is the (valid, checked) logfile specified in fitMarkers, or ""
  #plot can be "none", "fitted", "all"
  #plot.type can be "png", "emf", "svg", "pdf"
  #plot.dir can be NULL, "" or an absolute or relative path;
  #         a final "\" is added if needed
  #result is a list of 4 items: plot, plot.type, message, plot.dir
  result <- c(tolower(plot), "none", "", "")
  if (result[1] != "none") {
    result[2:3] <- checkPlotType(plot.type)
    if (result[[2]] == "none") result[1] <- "none"
  }
  if (result[1] != "none") {
    if (missing(plot.dir) || is.null(plot.dir)) {
      result[3] <- "plot.dir must be specified"
    } else {
      if (is.na(plot.dir[1])) plot.dir <- "" else
        plotdir <- gsub("(^ +)|( +$)", "", plot.dir[1]) #strip leading and trailing blanks
      if (plot.dir == "" && filePrefix != "") {
        plot.dir <- paste(filePrefix, "_plots/", sep="")
      }
      if (plot.dir == "") {
        plot.dir <- "plots/"
      } else {
        if (substring(plot.dir, nchar(plot.dir), nchar(plot.dir)) != "/")
          plot.dir <- paste(plot.dir, "/", sep="")
        #plot.dir <- paste("plots", plot.dir, sep="_")
      }
      #check if plot.dir exists and else create it:
      testdir <- substring(plot.dir, 1, nchar(plot.dir)-1)
      if (!(dir.exists(testdir) || dir.create(testdir))) {
        result[4] <- ""
        result[3] <- "cannot create plot.dir, use working directory instead"
      } else result[4] <- plot.dir
    }
  }
  names(result) <- c("plot", "plot.type", "message", "plot.dir")
  result
} #checkPlot

# get the directories for the plots (subdirs of the current working directory)
# with name "plots_" followed with the name of the logfile if any,
# and else with name "plots"
# if this cannot be created use the current working dir (return value is "")
# getPlotsDirectory <- function(logfile) {
#   if (!is.na(logfile))
#     logfile <- gsub("(^ +)|( +$)", "", logfile) #strip leading and trailing blanks
#   if (is.na(logfile) || logfile == "")
#     plotdir <- "plots" else {
#       logfile <- sub("[.][^.]*$", "", logfile, perl=TRUE) #remove the extension
#       if (logfile == "") plotdir <- "plots" else
#         plotdir <- paste("plots_", logfile, sep="")
#     }
#   dir.create(plotdir)
#   #return value:
#   if (file.exists(plotdir)) paste(plotdir, "/", sep="") else ""
# } #getPlotsDirectory

scr.info <- function(header) {
  #only for debugging split.screen problems
  cat(paste("***", header, "\n"))
  cat(paste("dev.cur: ", dev.cur(), "dev.list:\n"))
  dl <- dev.list(); print(dl)
  cat(paste("split.screen:", paste(split.screen(), collapse=" "),
            " cur.screen:", screen(), "\n"))
  cat(paste("cex=", par("cex"), " mex=", par("mex"), "\n"))
}

startPlot <- function(plot.type, plot.filename, ncol=1, nrow=1,
                       width=6.7, height=10.1) { #A4 with 2 cm margins
  # Starts a new device with a new page and calls split.screen if
  # more than 1 plot will be drawn
  #scr.info("start startPlot")
  switch (plot.type,
          pdf = pdf(file    =paste(plot.filename, ".pdf", sep=""),
                    width=width, height=height),
          png = png(filename=paste(plot.filename, ".png", sep=""),
                    width=width, height=height, units="in", res=96),
          emf = devEMF::emf(file=paste(plot.filename, ".emf", sep=""),
                    width=width, height=height),
          svg = svg(filename=paste(plot.filename, ".svg", sep=""),
                    width=width, height=height)
  )
  #scr.info("mid1 startPlot")
  close.screen(all.screens=TRUE) #This is essential!!
  #scr.info("mid2 startPlot")
  if (ncol>1 || nrow>1)
    #layout(matrix(1:(nrow*ncol),byrow=FALSE,ncol=ncol))  #adapt to the number of models tested
    suppressWarnings(split.screen(c(nrow, ncol)))
  # suppresses warning "calling par(new=TRUE) with no plot"
  #scr.info("end startPlot")
} #startPlot

savePlot <-function() {
  # Used for plotting the fitted or rejected model, with all populations
  #scr.info("start savePlot")
  dev.off()
  #scr.info("end savePlot")
} #savePlot

startPlotall <- function(model, orig_nonparents,
                         plot.type, plot.dir, mrknr, markername) {
  #Used for plotting the individual models
  #starts a new 2-, 4- or 8-panel page for drawing the models if needed
  #and selects the screen
  #print(paste("startPlotall: model=",model," all screens were ",paste(split.screen(), collapse=" ")))
  lnp <- length(orig_nonparents)
  ppp <- ifelse(lnp==1, 8, ifelse(lnp<=3, 4, 2)) #plots per page, in 2 columns
  rowcount <- ppp / 2
  # do we have to start a new page?
  if (model %% ppp == 1) {
    page <- (model-1) %/% ppp + 1
    plot.fname <- paste(plot.dir, "plots ", mrknr, " ", markername, " ",
                        formatC(page, width=2, flag="0"), sep="")
    startPlot(plot.type=plot.type, plot.filename=plot.fname,
               ncol=2, nrow=rowcount)
  }
  # select the screen (numbered by row, filled by column)
  nop <- (model - 1) %% ppp #current nr of plot on page (0 .. ppp-1)
  #print(paste("startPlotall: all screens are ",paste(split.screen(), collapse=" ")))
  #print(paste("model, ppp, nop=",model,ppp,nop," sel.screen: ",2 * nop + 1 + ((2 * nop) >= ppp)))
  #print(paste("model, ppp, nop=",model,ppp,nop," sel.screen: ",2 * (nop %% rowcount) + 1 + (nop >= rowcount)))
  screen(2 * (nop %% rowcount) + 1 + (nop >= rowcount))
  #print(paste("startPlotall: current screen is", screen()))
  # set cex and mex for the selected screen:
  par(cex=0.66, mex=1)
} #startPlotall

savePlotall <- function(model, orig_nonparents, modelcount) {
  #saves the 2-, 4- or 8-panel page if needed
  #modelcount is the number of models in a series of tried models:
  #4 if parentalPriors specified or if try.HW FALSE, else 8
  lnp <- length(orig_nonparents)
  ppp <- ifelse(lnp==1, 8, ifelse(lnp<=3, 4, 2)) #plots per page, in 2 columns
  if (ppp > modelcount) ppp <- modelcount
  if (model %% ppp == 0) {
    close.screen(all.screens=TRUE)
    savePlot()
  }
} #savePlotall


# ******************************************************
# functions that specify the 8 models tested in fitTetra
# ******************************************************

getModelName <- function (mutype, model.ps, ng) {
  if (length(model.ps$ptype) > 1) {
    s <- " pop"
  } else if (model.ps$ptype == "p.HW") s <- " HW" else s <- ""
  paste(getMuModelName(mutype, ng), s, sep="")
}

#getPtype <- function(model) {
#  c("p.free","p.free","p.free","p.free","p.HW","p.HW","p.HW","p.HW")[((model-1)%%8)+1]
#}

getMutype <- function(model) {
  rep(c(1, 2, 5, 6), 2)[((model-1) %% 8) + 1]
}

getAllModelNames <- function(ng) {
  mumodelnames <- getMuModelName(c(1,2,5,6), ng)
  c(mumodelnames,
    paste(mumodelnames,"HW"),
    paste(mumodelnames,"pop"),
    "none")
}

getMuModelInfo <- function(ng) {
  #Note: order must match gEMMax.mu.ratiodist !
  #Note: list items 1:7 correspond to mutype 0:6 !
  #   modelinfo <- list()
  #   modelinfo[[1]] <- list(
  #     name = paste("free, ng=",ng,sep=""),
  #     npar = ng) # 0=mu.free: ng means
  #   modelinfo[[2]] <- list(
  #     name = "b2",
  #     npar = 3) # 1 : c1, c2, f
  #   modelinfo[[3]] <- list(
  #     name = "b1",
  #     npar = 2) # c, f
  #   modelinfo[[4]] <- list(
  #     name = "b2,q2",
  #     npar = 5) # c1, c2, d1, d2, f
  #   modelinfo[[5]] <- list(
  #     name = "b1,q2",
  #     npar = 4) # c, d1, d2, f
  #   modelinfo[[6]] <- list(
  #     name = "b2,q",
  #     npar = 4) # c1, c2, d, f
  #   modelinfo[[7]] <- list(
  #     name = "b1,q",
  #     npar = 3) # c, d, f
  data.frame(
    name = c(paste("free, ng=", ng, sep=""),
             "b2", "b1", "b2,q2",  "b1,q2", "b2,q", "b1,q"),
    npar = c(ng, 3, 2, 5, 4, 4, 3),
    params = c("ng means", "c1, c2, f", "c, f", "c1, c2, d1, d2, f",
               "c, d1, d2, f", "c1, c2, d, f", "c, d, f"),
    stringsAsFactors = FALSE
  )
} #getMuModelInfo

getMuModelName <- function(mutype, ng) {
  #mutype in 0..6
  getMuModelInfo(ng)$name[mutype+1]
}

getMuModelNpar <- function(mutype, ng) {
  #mutype in 0..6
  getMuModelInfo(ng)$npar[mutype+1]
}

# ************************************
# supporting functions for CodomMarker
# ************************************

EMGaussExp <- function(D,sumD) {

  # Almost all work for E-step was already done in likelihood evaluation...
  # only some scaling is needed
  Z <- t(scale(t(D), center=FALSE, scale=sumD))
  attributes(Z)$"scaled:scale" <- NULL
  Z    # Z is matrix of posterior probabilities with
  #      as many rows as expanded bins and ng columns
} #EMGaussExp

# EMGaussExp.vectorized takes a vector of ratios (on asin(sqrt)) transformed scale)
# and a list psi with the model (ng mu's, sigma's and p's; mu's and sigma's on transformed scale)
# and produces a matrix with the probabilities that each y belongs to each of the ng peaks
EMGaussExp.vectorized <- function(y, psi) {
  ng <- length(psi$mu)
  ye <- rep(y, each=ng)
  dens.norm <- matrix(dnorm(ye, psi$mu, psi$sigma), ncol=ng, byrow=T)
  Zp <- dens.norm %*% diag(psi$p)
  Zs <- dens.norm %*% psi$p
  Z <- t(scale(t(Zp), center=FALSE, scale=Zs))
  attributes(Z)$"scaled:scale" <- NULL
  Z    # Z is matrix of posterior probabilities with as many row as length(y) and ng columns
}

EMMax.p.HW <- function(Z, w) {
  p <- apply(Z, 2, weighted.mean, w)
  #old version: hard-coded for ng=5 (tetraploid):
  #phw <- (4*p[1] + 3*p[2] + 2*p[3] + 1*p[4] + 0*p[5])/4
  #c(phw^4,4*phw^3*(1-phw),6*phw^2*(1-phw)^2,4*phw*(1-phw)^3,(1-phw)^4)
  ng <- ncol(Z)
  phw <- sum(p * ((ng-1):0)) / (ng-1)
  dbinom((ng-1):0, ng-1, phw) #over 10^7 iter with ng=5 0.12 sec FASTER than hard-coded version!
}

EMMax.p <- function(Z, BPW, p, ptype, pop.parents, segrArray) {
  #p and return value are matrices with (ploidy+1) columns and one row for
  #each population (in pop.parents? or only for populations with at least
  #one non-NA sample?), with the probabilities of the dosages
  #segrArray is a 3d array as produced by getPolysomicSegr / addErrorProb
  pnw <- p

  PopOrder <- nrow(p):1   # default order from last to first; assumes
  # that in pop.parents each F1 comes before its parents.
  # This is already done in CodomMarker or in the functions that call it.

  for (popi in PopOrder) { #treats the parents before their F1

    popidx <- (BPW[,2] == popi)
    popi.n <- sum(popidx)
    Zi <- Z[popidx,]
    Zi <- matrix(Zi, nrow=popi.n) # make sure it is a matrix, even with a single row
    wi <- BPW[popidx, 3]

    #We assume now that pop.parents is sorted correctly,
    #so the next lines are commented out:
    #if (ptype[popi] == "p.F1") {
    #  if (sum(pop.parents) > 0) parentsi <- pop.parents[popi,]
    #  else parentsi <- c(popi + 1, popi + 2)  # expect parents in populations
    #                                            just following F1 populations
    #}
    pnw[popi,] <-
      switch(ptype[popi],
             p.free  = apply(Zi, 2, weighted.mean, wi),
             #         estimate free p's by taking the average per column
             #         over all samples
             p.fixed = p[popi,],                        # p's are fixed
             p.HW    = EMMax.p.HW(Zi, wi),               # Hardy-Weinberg
             #        p.F1    = GetTetraF1Segregation(which(pnw[parentsi[1],] ==
             #                                          max(pnw[parentsi[1],]))-1,
             #                                        which(pnw[parentsi[2],] ==
             #                                          max(pnw[parentsi[2],]))-1)
             # p.F1    = GetTetraF1Segregation(which(pnw[pop.parents[popi,1],] ==
             #                                         max(pnw[pop.parents[popi,1],]))-1,
             #                                 which(pnw[pop.parents[popi,2],] ==
             #                                         max(pnw[pop.parents[popi,2],]))-1),
             #                  this calls GetTetraF1segregation with the dosages of P1 and P2
             #                  that have the maximum probability after the previous iteration
             p.F1    = segrArray[which(pnw[pop.parents[popi,1],] ==
                                          max(pnw[pop.parents[popi,1],])),
                                 which(pnw[pop.parents[popi,2],] ==
                                          max(pnw[pop.parents[popi,2],])),],
             #                  this gets the F1 segregation with the dosages
             #                  of P1 and P2 that have the maximum probability
             #                  after the previous iteration
             p.nodata = rep(1/ncol(pnw), ncol(pnw)),
             p.part  = stop("p.part should not occur") # apply(Zi, 2, weighted.mean, wi)
             #         ptype partly fixed(?), TODO: adapt in case pi=0
      )
  }
  pnw
} #EMMax.p

# EMMax.mu.ratiodist returns the predicted mu values (0:(ng-1))
# for one of the ratiodisttype models (mutype) 1:6 (0 is treated in EMMax.mu)
EMMax.mu.ratiodist <- function(z, y, mutype) {
  ng  <- ncol(z)
  nr  <- nrow(z)
  yw  <- rep(y, ng)
  wgt <- as.vector(z)
  x   <- rep(0:(ng-1), each=nr)

  #start values:
  startf <- 1    # the intrinsic ratio of the two signal responses
  startc <- 0.1  # the signal background parameter
  startd <- -0.1 # the quadratic (nonlinearity of signal response) parameter;
  #startd was 0; with mutype 3 a 0 causes problems and testing shows that -0.1 is not worse than 0
  #lower values:
  lowf <- 0.04 #lower-upper f in 0.04-25 instead of 0.001-1000 gives better convergence but (very slightly) lower LL
  lowc <- 0.001 #with lowc=0, many more convergence errors
  lowd <- -0.25 #was 0; with negative lowd often infinities or convergence errors
            #with lowd=0.001 a few less convergence errors than with lowd=0, but often
            #somewhat (and in some cases much) lower likelihoods
  #upper values:
  uppf <- 25   #see lowf
  uppc <- pi/2 #gives better convergence than uppc=Inf but slightly lower LL
               #pi/4 and pi/6 go further: still better convergence but lower LL
  uppd <- 1 #was Inf, values above 1 hardly change mu values
  # mutype is not 0 (treated separately in EMMax.mu):
  switch (mutype,
          { # 1 : c1, c2, f
            formul <- yw ~ asin(sqrt((c1+x)/(c1+x + c2+f*((ng-1)-x))))
            start <- list(c1=startc, c2=startc, f=startf)
            lower <- c(lowc, lowc, lowf)
            upper <- c(uppc, uppc, uppf)
          },
          { # 2 : c, f
            formul <- yw ~ asin(sqrt((c+x)/(c+x + c+f*((ng-1)-x))))
            start <- list(c=startc, f=startf)
            lower <- c(lowc, lowf)
            upper <- c(uppc, uppf)
          },
          { # 3 : c1, c2, d1, d2, f
            formul <- yw ~ asin(sqrt((c1+x+d1*x^2) /
                                       (c1+x+d1*x^2 + c2+f*((ng-1)-x)+d2*((ng-1)-x)^2 )))
            start <- list(c1=startc, c2=startc, f=startf, d1=startd, d2=startd)
            lower <- c(lowc, lowc, lowf, lowd, lowd)
            upper <- c(uppc, uppc, uppf, uppd, uppd)
          },
          { # 4 : c, d1, d2, f
            formul <- yw ~ asin(sqrt((c+x+d1*x^2) /
                                       (c+x+d1*x^2 + c+f*((ng-1)-x)+d2*((ng-1)-x)^2 )))
            start <- list(c=startc, f=startf, d1=startd, d2=startd)
            lower <- c(lowc, lowf, lowd, lowd)
            upper <- c(uppc, uppf, uppd, uppd)
          },
          { # 5 : c1, c2, d, f
            formul <- yw ~ asin(sqrt((c1+x+d*x^2) /
                                       (c1+x+d*x^2 + c2+f*((ng-1)-x)+d*((ng-1)-x)^2 )))
            start <- list(c1=startc, c2=startc, f=startf, d=startd)
            lower <- c(lowc, lowc, lowf, lowd)
            upper <- c(uppc, uppc, uppf, uppd)
          },
          { # 6 : c, d, f
            formul <- yw ~ asin(sqrt((c+x+d*x^2) /
                                       (c+x+d*x^2 + c+f*((ng-1)-x)+d*((ng-1)-x)^2 )))
            start <- list(c=startc, f=startf, d=startd)
            lower <- c(lowc, lowf, lowd)
            upper <- c(uppc, uppf, uppd)
          }
  ) #switch

  #first we try nls with the default Gauss-Newton algorithm.
  #This usually gives a better fit than the "port" algorithm but fails more often
  suppressWarnings({
    success <- tryCatch( {
      res.nls <- nls(formul, start=start, weights=wgt)
      TRUE }, error=function(x) {FALSE} )
    if (!success) {
      #if unsuccessful we try the "port" algorithm which allows to specify lower
      #and upper boundaries
      res.nls <- nls(formul, start=start, weights=wgt, algorithm="port",
                     lower=lower, upper=upper)
    }
    mu <- predict(res.nls, data.frame(x=0:(ng-1)), type="response")
  })
  mu
} #EMMax.mu.ratiodist

wmean <- function(w, y) { weighted.mean(y, w) } #for use in apply

EMMax.mu <- function(y, z, mutype) {
  #z: probabilities of all ng dosages for each sample
  if (mutype == 0) {
    apply(z, 2, wmean, y) # estimate free mu's
  } else EMMax.mu.ratiodist(z, y, mutype=mutype) # yw arsin-sqrt transformed,
  #                               but distances between mu's on (0,1) scale
  #                               according to ratios based on dosages
}

wmean2 <- function(zmu, y) {
  l <- length(zmu)
  weighted.mean( (y-zmu[l])^2, zmu[1:(l-1)] )
}

EMMax.sd.free <- function(y, z, mu) {
  sigma <- sqrt(apply(rbind(z, mu), 2, wmean2, y))
  sigma[sigma < 0.01] <- 0.01
  sigma
}

EMMax.sd.const <- function(y, z, mu) {
  sd <- sqrt(weighted.mean((rep(y, length(mu)) - rep(mu, each=length(y)))^2,
                            as.vector(z)))
  sigma <- rep(sd, length(mu))
  sigma[sigma < 0.01] <- 0.01
  sigma
}

EMMax.sd <- function(y, z, mu, sdtype, sd.fixed) {
  switch(sdtype,
    sd.const = EMMax.sd.const(y, z, mu),   # estimate constant sigma
    sd.free  = EMMax.sd.free(y, z, mu),    # estimate free sigma's
    sd.fixed = sd.fixed )
}

orderpsi <- function(psi, o) {
  # o is the order of the mu's: psi$mu[o] will be sorted in increasing order
  psi$mu <- psi$mu[o]
  psi$sigma <- psi$sigma[o]
  psi$p <- psi$p[, o, drop=FALSE]
  psi
} #orderpsi

# loglikf calculates the total loglik of vector y, given mu, sigma and p in psi,
# and population information in matrix BPW=Bins-Populations-Weights.
loglikf <- function(y, psi, BPW) {
  D <- t(sapply(y, dnorm, psi$mu, psi$sigma))
  # normal densities; dimension: n x ng, with n the number of observations or
  # number of non-empty bins;
  # in this way the evaluation of a minimum number of normal densities is achieved
  #Next we expand to all possible combinations of bin and population
  #(corresponding to the rows of BPW, i.e. we get n*npop rows of ng probabilities
  #instead of just n):
  De <- D[BPW[,1],]  # normal densities expanded for bins that contain
  #                    observations from multiple populations
  Pe <- psi$p[BPW[,2],] # prior probabilities expanded for populations
  DePe <- De*Pe
  rSDePe <- rowSums(DePe)  # likelihood of observation / bin mid over all populations
  LL <- sum(log(rSDePe) * BPW[,3])  # log likelihood of complete vector of
  #                           (binned) observations using weights=frequencies
  list(LL=LL, D=DePe, sumD=rSDePe)
} #loglikf

EMGaussMax <- function(yw, Z, BPW, mutype, sdtype, sd.fixed,
                       ptype, p, pop.parents, segrArray) {
  # yw is a vector of (non-NA) ratios (not binned) or the bin mids (binned)
  # Z is matrix of posterior probabilities with
  #      as many rows as expanded bins and ng columns
  #BPW is a data frame (long format) with columns for bin nrs (or 1:length(yw)
  #   if not binned), populations, weights
  #ptype is a character vector with the ptype for each population (row in
  #      pop.parents)
  #p is a matrix with ng columns and as many rows as there are populations
  #pop.parents gives the parent populations of each population
  #segrArray is the 3D array with F1 segregation ratios for ploidy, or NA

  W <- diag(BPW[,3]) #weights of the bins (or all 1's if not binned)
  ZW <- W %*% Z   # combine posterior probabilities with frequency weights
  #                 each row of ZW now sums not to 1 but to BPW[,3]
  #                 equivalent: ZW <- BPW[,3] * Z

  ZW <- rowsum(ZW, BPW[,1], reorder=TRUE)    # collapse weights from different
  #    populations with values in identical bins to arrive at level of bin mids

  mu    <- EMMax.mu(yw, ZW, mutype)
  sigma <- EMMax.sd(yw, ZW, mu, sdtype, sd.fixed)
  p     <- EMMax.p(Z, BPW, p, ptype, pop.parents, segrArray)

  list(mu=mu, sigma=sigma, p=p)
  #comment from earlier version:
  #so mu is optimized first, these are used to optimize sigma
  #the new p are derived from the current p's of the samples
} #EMGaussMax

#EMGaussMix from populations version:
EMGaussMix <- function(y, ng=5, pop=rep(1, length(y)),
                       ptype="p.free", mutype=0,
                       sdtype="sd.const", sd.fixed,
                       p.start=matrix(1/ng, ncol=ng, nrow=nrow(pop)),
                       mu.start=seq(0.02, pi/2-0.02, length.out=ng),
                       sd.start=rep(0.075, ng),
                       npar, maxiter, maxn.bin, nbin,
                       pop.parents, segrArray)  {

  #y is on the transformed (arcsine-square root) scale, like mu.start and sd.start
  isna <- is.na(y)
  yw <- y[!isna] #NOT done in CodomMarker!
  n <- length(yw)
  popw <- pop[!isna] #integer vector, indexing rows of pop.parents
  popw <- as.factor(popw) #needed for table(), keeps same values but omits
  #       rows of pop.parents not indexed by the remaining pop
  binned <- FALSE
  B <- 1:length(yw)
  P <- popw
  W <- 1
  BPW <- cbind(B,P,W)  # Bins, Populations, Weights, nrow=length(yw)
  if (n > maxn.bin) {
    binned <- TRUE
    brks <- seq(0, pi/2, length.out=nbin+1)
    mids <- brks[1:nbin] + 0.5*(pi/2)/nbin
    brks[nbin+1] <- brks[nbin+1] + 0.00001 #makes sure last one above pi/2
    bin.i <- findInterval(yw, brks) #for each yw, get the bin nr
    Wt <- table(factor(bin.i), popw)  # frequencies in bins, per population;
    #                   table reports only bins with at least one observation
    #                   and populations that have at least one sample in yw
    #                   populations in columns, bin nrs in rows
    bins <- as.numeric(rownames(Wt)) #which bins have at least one obs
    #pops <- as.numeric(as.factor(colnames(Wt))) #same as 1:ncol(Wt)??
    #                          Different when more than 9 populations!!
    #Actually populations may be missing due to missing data; correct is:
    pops <- as.numeric(colnames(Wt)) #20160927 change from line above
    npops <- length(pops); nbins <- length(bins)
    yw <- mids[bins]           # bin mids, to be used as response data

    Bb <- rep(1:nbins, times=npops)
    Pb <- rep(pops, each=nbins)
    Wb <- as.vector(Wt)
    BPW <- cbind(Bb, Pb, Wb) # Bins, Populations, Weights; integer matrix
    BPW <- BPW[BPW[,3] > 0 , ] # Save only rows with at least one observation
    #                            (i.e. weight > 0 for that bin/pop combi)
  }

  result <- list (loglik=NA, npar=npar, AIC=NA, BIC=NA, psi=NA, post=NA, nobs=n,
                  iter=NA, message="")
  psi <- list(mu=mu.start, sigma=sd.start, p=p.start)
  o <- order(psi$mu); psi <- orderpsi(psi, o)
  LLnw <- loglikf(yw, psi, BPW) # returns not only LogLik but also
  #                               calculation results to be used in E-step
  iter <- 0
  tryCatch( {
    repeat { # start of EM loop
      iter <- iter+1
      LL <- LLnw
      Z <- EMGaussExp(LL$D, LL$sumD)
      result$message <- tryCatch( {
        psi <- EMGaussMax(yw=yw, Z=Z, BPW=BPW, mutype=mutype,
                          sdtype=sdtype, sd.fixed=sd.fixed,
                          ptype=ptype, p=p.start,
                          pop.parents=pop.parents, segrArray=segrArray)
        #                 p.start is used, only for the case that ptype="p.fixed"
        ""
      }, error=function(ex) {
        paste(ex$message, "in EMGaussMix")
      } )
      if (result$message != "") break
      LLnw <- loglikf(yw, psi, BPW)
      if (is.na(LLnw$LL) || (max(LLnw$LL - LL$LL) < 0.0000001) ||
          (maxiter > 0 && iter > maxiter) ) break
    } # end of EM loop

    if (result$message == "") {
      if ( is.na(LLnw$LL) ) result$message <- "LLnw==NA in EMGaussMix" else
      if (max(LLnw$LL - LL$LL) < 0.0000001) {
        o <- order(psi$mu); psi <- orderpsi(psi, o); Z <- Z[,o]
        NAmat <- matrix(rep(!isna, ng), ncol=ng, byrow=FALSE)
        z <- matrix(nrow=length(y), ncol=ng)
        BPW <- cbind(B, P, W)  # Bins, Populations, Weights, unbinned data
        LLnw <- loglikf(y[!isna], psi, BPW) # loglikf call with original data!!
        Z <- EMGaussExp(LLnw$D, LLnw$sumD)
        z[NAmat] <- Z
        result$post <- z
        result$loglik <- LLnw$LL
        result$AIC <- -2*LLnw$LL + 2*npar
        result$BIC <- -2*LLnw$LL + log(n)*npar
        result$psi <- psi
        result$iter <- iter
      } else result$message <- "iter>maxiter in EMGaussMix"
    }
  }, error=function(ex) {
    result$message <- paste(ex$message," in EMGaussMix",sep="")
  } )
  result
} #EMGaussMix

# fitOneDist implements the case that only one component distribution is present
fitOneDist <- function(y, sdtype="sd.const", sd.fixed, npar)  {
  isna <- is.na(y)
  yw <- y[!isna]
  n <- length(yw)
  mu <- mean(yw)
  if (sdtype=="sd.fixed") {
    sigma <- sd.fixed
  } else {
    sigma <- sd(yw)
  }
  p <- 1.0
  psi <- list(mu=mu, sigma=sigma, p=p)
  llnw <- sum(log(sapply(yw, dnorm, psi$mu, psi$sigma)))
  NAmat <- matrix(!isna, ncol=1,  byrow=F)
  z <- matrix(nrow=length(y), ncol=1); z[NAmat] <- 1.0
  if (is.na(llnw)) llnw <- -1.0e99
  AIC <- -2*llnw + 2*npar
  BIC <- -2*llnw + log(n)*npar
  list(loglik=llnw, npar=npar, AIC=AIC, BIC=BIC, psi=psi,
       post=z, nobs=n, iter=1, message="")
}

# *******************************************************
# Functions to re-order a pedigree data frame,
# used in fitOneMarker to re-order the pop.parents data.frame
# *******************************************************

sortPedigree <- function(ped, colInd, colPar1, colPar2,
                         parentsFirst=TRUE, semiFounders=FALSE) {
  #based on sortPedigree in Pedimap
  #function sorts a data frame containing pedigree info such
  #that parents always occur before (default) or always after their offspring
  #with minimal shuffling
  #ped: the pedigree data frame to be ordered. Must have one column each for
  #     the individual ID's, parent1 and parent2 and may contain additional
  #     columns. Column names are free (not used).
  #colInd, colPar1, colPar2: the numbers of the columns with the names of the
  #        individuals, first and second parents. Individuals may not have
  #        missing values, but for (semi-)founders one or both parents are NA.
  #        All individuals must be different.
  #parentsFirst: if TRUE (default), parents will be sorted before any progeny,
  #              if FALSE, after all progeny
  #semiFounders: if FALSE (default) no semiFounders are allowed: parent1 and
  #              parent2 must either both be NA or both be non-NA
  #Return value: character (error message) if the input is incorrect,
  #              else the ped data frame with the rows reordered if necessary,
  #              and with the colInd, colPar1 and colPar2 as factors with
  #              the same set of levels.

  #first checks of the input (except circular pedigrees):
  if (nrow(ped) == 0) return(ped)
  ped[, colInd] <- as.character(ped[, colInd])
  ped[, colPar1] <- as.character(ped[, colPar1])
  ped[, colPar2] <- as.character(ped[, colPar2])
  if (sum(is.na(ped[, colInd])) > 0)
    return ("missing individual/population names not allowed")
  if (length(unique(ped[, colInd])) < nrow(ped))
    return("duplicate individual/population names not allowed")
  if (length(setdiff(ped[, colPar1], c(NA_character_, ped[, colInd]))) > 0 ||
      length(setdiff(ped[, colPar2], c(NA, ped[, colInd]))) > 0  )
    return("all parents should also be listed as individuals/populations")
  if (!semiFounders &&
      sum(xor(is.na(ped[, colPar1]), is.na(ped[, colPar2]))) > 0)
    return("semi-founders not allowed")
  if (sum(ped[, colPar1] == ped[, colInd], na.rm=TRUE) > 0 ||
      sum(ped[, colPar2] == ped[, colInd], na.rm=TRUE) > 0  )
    return("some individuals/populations are listed as their own parents")
  if (!parentsFirst) ped <- ped[nrow(ped):1,] #reverse ped
  #double size of ped for sorting:
  pedsize <- nrow(ped)
  ped <- rbind(ped, ped) #double number of rows, for sorting
  for (i in 1:length(ped)) ped[(pedsize+1):nrow(ped), i] <-
    rep(NA, pedsize) #lower half empty

  #first step: place first founder at line 1
  founders <- which(is.na(ped[1:pedsize, colPar1]) &
                      is.na(ped[1:pedsize, colPar2]))
  if (length(founders) == 0)
    return ("pedigree contains no founders")
  if (pedsize == 1) {
    ped <- ped[1,]
    ped[, colInd] <- factor(ped[, colInd]) #both parents are NA
    return(ped)
  }
  if (founders[1] > 1) {
    #move all individuals from 1 to founders[1]-1 to end
    #and then move up all indivs to 1
    ped[(pedsize+1):(pedsize + founders[1] - 1),] <- ped[1:(founders[1] - 1),]
    ped[1:pedsize,] <- ped[founders[1]:(pedsize + founders[1] - 1),]
  }
  #second step: loop, placing in each pass the next individual whose parents
  #(if there are any)
  pass <- 2
  while(pass < pedsize) {
    #for sorting, pass<ICount-1 would be sufficient, but in that case
    #the last Indiv might have a non-existing parent which would not be noticed
    #find first indiv i that can be placed at position pass,
    #i.e. whose parents (if present) are already placed
    #For all ind from position pass to end:
    #find position of parents (0 if parent NA, NA if parent not yet placed):
    posP1 <- match(ped[pass:pedsize, colPar1], ped[1:(pass-1), colInd])
    posP1[which(is.na(ped[pass:pedsize, colPar1]))] <- 0
    posP2 <- match(ped[pass:pedsize, colPar2], ped[1:(pass-1), colInd])
    posP2[which(is.na(ped[pass:pedsize, colPar2]))] <- 0
    #find the next individuals (from position pass) that could now be placed
    nxt <- which((posP1 < pass) & (posP2 < pass))
    if (length(nxt) == 0)
      return("circular references in pedigree")
    if (nxt[1] == 1) {
      #individual ped[pass,] is already in place and possibly a number of the
      #next individuals as well, find out which:
      last.ok <- max(which(nxt == 1:length(nxt)))
      #we set pass to the next individual after that:
      pass <- pass + last.ok
    } else {
      #nxt[1] > 1: the next individual that can now be placed is
      #ped[pass + nxt[1] - 1,]
      #move all individuals from pass to pass + nxt[1] - 2 to end
      #and then move up all indivs to pass
      #ped[(pedsize+1):(2*pedsize - pass - nxt[1] + 1),] <-
      ped[(pedsize + 1):(pedsize + nxt[1] - 1),] <-  #nxt[1]-1 indiv
        ped[pass:(pass + nxt[1] - 2),]               #nxt[1]-1 indiv
      ped[pass:pedsize,] <-                          #pedsize-pass+1 indiv
        ped[(pass + nxt[1] - 1):
              (pedsize + nxt[1] - 1),]                 #pedsize-pass+1 indiv
      #the individual at pass is now ok, move to next:
      pass <- pass + 1
    }
  } #while pass
  ped <- ped[1:pedsize,]
  if (!parentsFirst) ped <- ped[nrow(ped):1,] #reverse ped back again
  ped[, colInd] <- factor(ped[, colInd])
  ped[, colPar1] <- factor(ped[, colPar1], levels=levels(ped[, colInd]))
  ped[, colPar2] <- factor(ped[, colPar2], levels=levels(ped[, colInd]))
  ped
} #sortPedigree



isPedigreeSorted <- function(ped, colInd, colPar1, colPar2,
                             parentsFirst=TRUE, semiFounders=FALSE) {
  #parameters as sortPedigree
  #return value: error message if input incorrect (see sortPedigree),
  #              else TRUE or FALSE
  ped[[colInd]] <- as.character(ped[[colInd]])
  ped[[colPar1]] <- as.character(ped[[colPar1]])
  ped[[colPar2]] <- as.character(ped[[colPar2]])
  if (sum(is.na(ped[[colInd]])) > 0)
    return ("missing individual names not allowed")
  if (length(unique(ped[[colInd]])) < nrow(ped))
    return("duplicate individual names not allowed")
  if (length(setdiff(ped[[colPar1]], ped[[colInd]])) > 0 ||
      length(setdiff(ped[[colPar2]], ped[[colInd]])) > 0  )
    return("all parents should be listed as individuals")
  if (!semiFounders &&
      sum(xor(is.na(ped[[colPar1]]), is.na(ped[[colPar2]]))) > 0)
    return("semi-founders not allowed")
  if (parentsFirst) {
    parline <- match(ped[[colPar1]], ped[[colInd]])
    if (sum(parline > (1:nrow(ped)), na.rm=TRUE) > 0) return(FALSE)
    parline <- match(ped[[colPar2]], ped[[colInd]])
    if (sum(parline > (1:nrow(ped)), na.rm=TRUE) > 0) return(FALSE)
    return(TRUE)
  } else {
    parline <- match(ped[[colPar1]], ped[[colInd]])
    if (sum(parline < (1:nrow(ped)), na.rm=TRUE) > 0) return(FALSE)
    parline <- match(ped[[colPar2]], ped[[colInd]])
    if (sum(parline < (1:nrow(ped)), na.rm=TRUE) > 0) return(FALSE)
    return(TRUE)
  }
} #isPedigreeSorted

getModelFromFile <- function(marker, modeldata) {
  #Created for debugging, may be exported in the future
  #marker: vector of marker numbers or names occurring in modeldata
  #modeldata: either the modeldata data.frame returned by fitOneMarker,
  #     or the modelfile returned by fitMarkers (RData or text file)
  #Value:
  #A list with for each marker one element (named according to marker).
  #Each of these elements is a list with 3 elements:
  # $mu: the (ploidy+1) means of the model on the transformed scale
  # $sd: the (ploidy+1) standard deviations on the transformed scale
  #      (all of these sd's are equal under the currently used models)
  # $P: a matrix with (ploidy+1) columns for dosages 0 .. ploidy and
  #     one row per population; the rownames are the population names.
  #     This matrix contains the mixing proportions of the mixture components
  if (!is.data.frame(modeldata)) {
    #modeldata is a filename
    ext <- tolower(tools::file_ext(modeldata))
    if (substr(ext, 1, 3) == "rda") {
      vars <- load(modeldata)
      if (length(vars) != 1 || vars != "modeldata")
        stop("getModelFromFile: RData file should contain one variable, called modeldata")
    } else {
      modeldata <- read.table(modeldata, header=TRUE, sep="\t")
    }
  }
  if (!is.data.frame(modeldata) ||
      names(modeldata)[1] != "marker" ||
      length(unique(modeldata$marker)) != nrow(modeldata))
    stop("getModelFromFile: modeldata should be the modeldata from fitOneMarker ",
         "or the modelfile from fitMarkers")
  mucol <- which(names(modeldata) == "mutrans0")
  sdcol <- which(names(modeldata) == "sdtrans0")
  ng <- sdcol - mucol
  pcol <- sdcol + ng
  mu0col <- which(names(modeldata) == "mu0") # first after the P columns
  popcount <- (mu0col - pcol) / ng
  if (popcount == 1) popnames <- "" else {
    popnames <- character(popcount)
    for (p in 1:popcount) {
      pn <- names(modeldata)[sdcol + p*ng] #the nulliplex col of population p
      popnames[p] <- substring(pn, 3, nchar(pn)-2) #pn is P_population_0
    }
  }
  result <- list()
  for (m in seq_along(marker)) {
    if (is.numeric(marker)) {
      mdat <- modeldata[modeldata$marker == marker[m],]
    } else mdat <- modeldata[modeldata$MarkerName ==  marker[m],]
  if (nrow(mdat) != 1)
    stop(paste("getModelFromFile: marker ", m, " (", marker[m],
               ") does not occur or occurs multiple times in modeldata",
               sep=""))
  mu <- as.numeric(mdat[mucol:(mucol+ng-1)])
  sd <- as.numeric(mdat[sdcol:(sdcol+ng-1)])
  P <- matrix(ncol=ng, nrow=popcount, dimnames=list(popnames, 0:(ng-1)))
  for (p in 1:popcount) P[p,] <- as.numeric(mdat[(pcol+(p-1)*ng):(pcol+p*ng-1)])
  result[[m]] <- list(mu=mu, sd=sd, P=P)
  }
  names(result) <- marker
  result
} #getModelFromFile

calcScores <- function(y, psi, pop=NA) {
  #Created for debugging, may be exported in the future
  #y: data frame with SampleNames in column 1 and asin(sqrt(x))- transformed
  #    ratios in column 2
  #psi: list with components mu, sd and P as returned by getModelFromFile
  #pop: data frame with SampleNames in column 1 and populationID (matching the
  #     rownames of psi$P) in column 2. Ignored and not required if P has only
  #     one row
  #Result value:
  #A data.frame with one row for each sample in y and columns:
  # $SampleName
  # $pop
  # $transf_ratio
  # $P0 .. $P<ploidy>
  # $maxP
  # $maxgeno
  if (is1NA(pop)) {
    if (nrow(psi$P) == 1) {
      pop <- data.frame(
        SampleName=y[,1],
        Population=rownames(psi$P))
    } else stop("calcScores: pop must be specified with multiple populations")
  }
  if (anyNA(y[,1]) || length(unique(y[,1])) != nrow(y))
    stop("calcscores: all SampleNames in y must be unique and not NA")
  matchY <-
  pop <- pop[match(y[,1], pop[,1]),2]
  if (anyNA(pop))
    stop("calcScores: all SampleNames from y must also occur, ",
         "with non-NA population ID, in pop")
  ng <- length(psi$mu)
  probs <- matrix(NA_real_, nrow=nrow(y), ncol=ng,
                  dimnames=list(y[,1], paste("P",0:(ng-1),sep="")))
  for (pp in 1:nrow(psi$P)) {
    popname <- rownames(psi$P)[pp]
    poprows <- pop == popname
    if (sum(poprows) > 0) {
      probs[poprows,] <-
        EMGaussExp.vectorized(y[,2], list(mu=psi$mu, sd=psi$sd, P=psi$P[pp,]))
      #matrix, per sample ng numbers: probs that sample belongs to each peak
    }
  }
  #calculate maxgeno, maxP, geno:
  result <- data.frame(
    SampleName=y[,1],
    pop=pop,
    transf_ratio=y[,2],
    probs
  )
  maxPna <- apply(probs, 1, max) #find the maxP of each row of p-values
  result$maxP <- maxPna[!is.na(maxPna)] # omit missing values
  result$maxgeno <- max.col(probs) - 1  # for every row (sample) the
  #                max peak (0..ploidy instead of 1..ng)
} #calcScores
