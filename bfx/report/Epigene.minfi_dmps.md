### Differentially Methylated Positions (using Minfi)

Differentially methylated positions (DMPs) are identified using the Minfi package.

The p-value cut-off chosen is $pvalue$ and methylation difference threshold is $delta_beta_threshold$.  Differentially methylation is considered to be significant if p-value is lower than than p-value cut-off.

The p-values can be adjusted for multiple testing using different methods in order to reduce false positives from the result.  P-value adjustment method chosen in this run was $padjust_method$.

The results for each contrast can be found below:

$data$


* Cpg Site: site of CpG region
* Avg Case Beta: average methylation of case (sample or Mut)
* Avg Control Beta: average methylation of control (WT)
* Avg Delta Beta: average methylation difference (ie. Avg Case Beta - Avg Control Beta)
* Chromosome: chromosomal position
* Gene: gene where CpG site is found
* P-value: raw pvalue (pvalue before adjustment)
* Adj.P-value: pvalue after adjustment
