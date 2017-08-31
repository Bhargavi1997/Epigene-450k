#### Failed Probes

Our data is first filtered to find failed probes (probes that have failed to hybridize) using the [minfi](http://bioconductor.org/packages/release/bioc/html/minfi.html) function *detectionP* with a percent failed threshold of $percent_failed_threshold$ and a p-value cut-off of $pvalue$.  A list of the removed failed probes can be found [here](data/quality_control/failed_probes.csv).

A detection p-value is calculated for every genomic location in the sample.  Failed p-values are identified using p-value cut-off and percent of failed p-values of each probes is calculated.  Small p-values indicate a good position and large p-values indicate poor quality sample.  Probes with a percent of failed p-values above the given failed threshold are then removed.

A total of $failed_removed$ probes were removed.


