Introduction
------------

This document contains the description of the current MUGQIC Epigene $data_type$ analysis. The information presented here reflects the current state of the analysis as of $date$.

This pipeline is a modified implementation of the open source [R Package minfi](http://bioconductor.org/packages/release/bioc/html/minfi.html) for the Illumina $data_type$ platform. The diagram below shows the structure of the pipeline.
 
![**_Epigene pipeline structre_**](images/flowchart.png)



The pipeline takes raw .IDAT files and performs a series of optional preprocessing steps on them. Several normalization methods are also available to be used before the data undergoes analysis for differentially methylated positions and regions.  The final output is a list of CpG sites and their methylation levels, a list of differentially methylated positions and a list of differentially methylated regions.


Analysis and Results
--------------------

A total of $sample_number$ samples have been analyzed.  Among these $sample_number$ samples, there are different group types:

$group_list$

Analysis can be done selectively to compare any or all of the group types specified by the contrasts. The contrasts chosen for the current analysis are described below:

| Contrast name | Number of controls | Number of cases | 
-----:|-----:|-----:
$contrast_table$
