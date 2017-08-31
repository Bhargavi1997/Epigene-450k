### Differentially Methylated Regions

DMRs are stretches of DNA that have consistent patterns of methylation differences.  DMRs are identified by comparing the two group types of the dataset using the [bump hunting algorithm](http://ije.oxfordjournals.org/content/41/1/200.full).

The p-values are derived using $nullMethod$ with B = $permutations$ resamples and a delta beta threshold of $delta_beta_threshold$ is applied to identify significant DMRs. 

The results for each contrast can be found below:

$data$

* start: start of the genomic location of the differentially methylated regions
* end: end of genomic locations of the differentially methylated regions
* value: average difference in methylation in the bump
* fwer column: family-wise error rate(fwer) of the regions estimated by the permutation scheme
* p.value: pvalue assigned to the region
* L: number of CpGs in the differentially methylated region
* Gene: gene where CpG site is found
