### Differentially Methylated Positions (using Limma)

Differentially methylated positions are identified using the Limma package.

The p-value cutoff chosen is $limma_pvalue$ and methylation difference threshold is $limma_delta_beta$. Differential methylation is considered to be significant if p-value is lower than the p-value cut-off.

The p-values can be adjusted for multiple testing. P-value adjustment method chosen in this run was $limma_padjust$.

$limma_data$

* CpG Site: site of CpG region
* logFC: estimate of the log2-fold-change
* AveExpr: average log of difference in methylation for the probe over all arrays and channels
* AvgDeltaBeta: average methylation difference (ie. Avg Case Beta - Avg Control Beta)
* t: moderated t-statistic
* P.Value: raw p-value
* adj.P.Val: adjusted p-value or q-value
* B: log-odds that the gene is differentially expressed
* Gene: gene where CpG site is found
* Chromosome: chromosomal position
