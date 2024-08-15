# FitPoly: Genotype Calling for Bi-Allelic Marker Assays

FitPoly is able to genotype polyploid bi-allelic markers using signal intensities from SNP arrays. 'fitPoly' assigns genotypes (allele dosages) to a collection of polyploid samples based on these signal intensities. 'fitPoly' replaces the older package 'fitTetra' that was limited to only tetraploid populations whereas 'fitPoly' accepts any ploidy level. Version 4.0.0 onwards also include all functionalities that were added in the R package fitPolyTools.

Main functionalities of fitPoly:

-   fitPoly can fit any ploidy populations following (or not) Hardy-Weinberg equilibrium, F1 populations, with parental or sample prior information.

-   Data import from Affymetrix and Illumina SNP arrays into fitPoly format.

-   Marker and sample filtering based on signal intensities.

-   Quality checks using analysis of dosage segregation of F1 populations.

-   Other utiliy functions that can be found in the vignettes.

FitTetra reference: [Genotype calling in tetraploid species from bi-allelic marker data using mixture models. Voorrips RE, Gort G, Vosman B (2011)](http://dx.doi.org/10.1186/1471-2105-12-172)

FitTetra 2.0 reference: [FitTetra 2.0 -- improved genotype calling for tetraploids with multiple population and parental data support](https://doi.org/10.1186/s12859-019-2703-y)
