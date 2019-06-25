# spartina.genetics

These are the associated files and R code for the SNP dataset that are a part of "Microgeographic genetic divergence within a coastal foundation species contributes to ecosystem function" Collaboration with Robyn Zerebecki, Randall Hughes, Torrance Hanley, Chris Nice and others.

## data

1) "allpops_individualIDs_subset.csv" and "allpops_individualIDs.csv"  
Description: set of 324 individuals and 309 individuals used in following analyses.  

2) "spartinaNov2017.called.subset.mpgl"  
Description: set of genotype likelihoods for all individuals. values range from 0-2; HOW WAS THIS GENERATED? 1st column = location; columns 2-310: 309 individuals. 2735 rows = 2735 loci. 

3) "t.s.pops.beagle.gz"  
Description: beagle-formatted file of genotyple likelihoods generated by ANGSD. read in using snpStats::read.beagle(). Each individual X SNP combination occurs across 3 columns (AA, AB, BB).

