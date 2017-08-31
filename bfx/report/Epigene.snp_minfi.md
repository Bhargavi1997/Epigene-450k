#### SNP Probes

The data then is filtered to search for SNPs using the [minfi](http://bioconductor.org/packages/release/bioc/html/minfi.html) function *getSnpInfo* followed by function *dropLociWithSnp* which removes the SNPs from the data.  A list of the removed SNPs can be found [here](data/quality_control/snp_probes.csv).
A total of $snp_removed$ probes were removed.

