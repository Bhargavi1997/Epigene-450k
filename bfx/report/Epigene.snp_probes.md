#### SNP Probes

The presence of SNPs (single nucleotide polymorphisms) inside the probe body or at the nucleotide extension can have important consequences on the downstream analysis.  Presence of SNPs at interrogated CpG sites in probes can lead to incorrect methylation measurements.  A base change at the interrogated CpG can cause an otherwise methylated site to be considered unmethylated.  

Probes that contain either a SNP directly at an interrogated CpG site or within 10 bp of an interrogated CpG island are dropped.  A list of all the removed probes that contained SNPs can be found [here](data/quality_control/snp_probes_minfi.csv).

A total of $snp_removed$ probes were removed.


