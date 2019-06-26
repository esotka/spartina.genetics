# spartina.genetics

These are the associated files and R code for the SNP dataset that are a part of "Microgeographic genetic divergence within a coastal foundation species contributes to ecosystem function" Collaboration with Robyn Zerebecki, Randall Hughes, Torrance Hanley, Chris Nice and others.

## data

1) ***inds/***  
Description: a set of files with individuals used in the analyses.  
a. "allpops_individualIDs_subset.csv" - set of 309 individuals - all sites.  
b. "indIDS.T.S.POPS.txt" - set of six marshes with tall-vs-short comparisons.  
c. "inds.FL" "inds.RISRIT" "inds.SBITBI" "inds.SFBTFB" "inds.SWSSWT" "inds.WESWET" - six marshes with tall-vs-short comparisons  

2) "Spartina_SNP_SiteID.csv"  
Description: meta-data for each location.  

3) "spartinaNov2017.called.subset.mpgl"  
Description: set of genotype likelihoods for all individuals. values range from 0-2; We converted the phred-scale genotype likelihoods (from samtools/bcftools) per SNP-sample combination into probabilities that summed to 1, and then converted these to a single value that ranges from 0 to 2, where 0, 1 and 2 represent the highest probability of a homozygote, heterozygote, and alternative homozygote, respectively. 1st column = location; columns 2-310: 309 individuals. 2735 rows = 2735 loci. 

4) "t.s.pops.beagle.gz"  
Description: genotype likelihoods (from samtools/bcftools) in beagle-formatted file generated by angsd. Each individual X SNP combination occurs across 3 columns (AA, AB, BB).

5) "fst/"  
Description: fst output from angsd.  

6) "admix.runs.TvS.POPS/" "admix.runs.ALLPOPS"  
Description: NGSadmix runs for all populations and tall-vs-short comparisons.

### Figure: PCAs  
Based on genotype likelihoods of individuals embedded in the .mpgl file.  Uses the prcomp() to generate PCAs. there are four versions:  
* "R/map of spartina sites.R" and "output/map of spartina sites.pdf" - this includes all sites, and is coupled with a map of the east coast of the US.  
* "R/shortVStall-pca-allDATA.R" and "shortVStall-pca-allDATA.pdf" - this compares PCAs of tall-vs-short individuals within each of 6 marshes.     
* "R/shortVStall-pca_wout10.R" and "shortVStall-pca_wout10.pdf" - same as previous, but side-by-side compares when all data or the 10% outlier loci (by Fst) are removed.  
* "R/shortVStall-pca-&admix_allData.R" and "shortVStall-pca&admix_allDATA.pdf" - This does a side-by-side comparison of admixture plots and PCA for 6 marshes.    

### Figure: admixture plots  
Based on NGSadmix analyses in angsd.  
* "R/ngsAdmix.pretty.R" and "output/ngsAdmix.pretty.pdf"  - all sites and 309 individuals  
* "R/ngsAdmix.T.S.POPS.pretty.R" and "output/ngsAdmix.T.S.POPS.pretty.pdf" - only the 6 marshes in which tall-vs-short occurs.   

### Figure: Lnl and deltaK of admixture plots  
Based on 5 independent runs of NGSadmix in angsd  
* "R/deltaK.R"and "output/deltaK.pdf"  


