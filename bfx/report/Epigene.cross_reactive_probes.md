#### Cross Reactive Probes

Cross reactivity occurs when a probe hybridizes to not only the targetted region but to other regions of the genome (off-target areas).  Since cross reactive probes map to multiple genomic region, it is unclear which genomic region gives rise to the measured methylation status. It is recommended to remove these probes from the dataset before any subsequent analysis. 

A list of [cross reactive probes](data/quality_control/cross_reactive_probes.csv) for 450K data was provided by [Weksberg Lab](http://www.sickkids.ca/Research/Weksberg-Lab/Publications/index.html).  The list is obtained from the Chen et al., 2013 paper titled *Discovery of cross reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray* and contains a total of 30969 probes.

Cross reactive probes found in the list are removed from the 450K data.


