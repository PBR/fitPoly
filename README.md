# FitPoly: Genotype Calling for Bi-Allelic Marker Assays

FitPoly is able to genotype polyploid bi-allelic markers using signal intensities from SNP arrays. 'fitPoly' assigns genotypes (allele dosages) to a collection of polyploid samples based on these signal intensities. 'fitPoly' replaces the older package 'fitTetra' that was limited to only tetraploid populations whereas 'fitPoly' accepts any ploidy level.

Additional functionalities have been added to fitPoly since the development of fitTetra:

-   fitPoly is now able to fit not only F1 populations, but also diversity panels (assuming Hardy-Weinberg equilibrium) and linked F1 populations.

-   Data import from Affymetrix and Illumina SNP arrays into fitPoly format.

-   Marker and sample filtering based on signal intensities.

-   Quality checks using analysis of dosage segregation of F1 populations.

-   Other utiliy functions that can be found in the vignettes.

FitTetra reference: [Genotype calling in tetraploid species from bi-allelic marker data using mixture models. Voorrips RE, Gort G, Vosman B (2011)](<http://dx.doi.org/10.1186/1471-2105-12-172>).
